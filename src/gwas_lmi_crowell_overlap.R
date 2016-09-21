#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

rutils::GenerateMessage("Calculate overlaps between LMI and Crowell GWAS")

# Load GWAS files
gwas.crowell.file <- "data/gwas_crowell_genes.Rds"
gwas.lmi.file <- "data/gwas_lmi_genes.Rds"

rutils::GenerateMessage("Loading files")
gwas.crowell <- readRDS(gwas.crowell.file)
gwas.lmi <- readRDS(gwas.lmi.file)

# deal with redundant lmi regions
rutils::GenerateMessage("Generating regions: LMI")
gwas.lmi <- gwas.lmi[, .(Chr, Start, End, Trait, Place, Group, bin.name)]
gwas.lmi <- gwas.lmi[!is.na(End)]
setkey(gwas.lmi, Chr, Start, End)
gwas.lmi.long <- gwas.lmi[, .(
  position = c(Start:End)
), by = .(Chr, Start, End, Trait, Place, Group, bin.name)]
gwas.lmi.long[, c("Start", "End") := NULL]

rutils::GenerateMessage("Sorting intermediate table")
setkey(gwas.lmi.long, Chr, position)

rutils::GenerateMessage("Grouping LMI regions by Chr and position")

system.time(
  gwas.long.tmp <- gwas.lmi.long[Chr == 1, .(
    trait = list(unique(Trait)),
    subpop = list(unique(Group))
  ), by = .(Chr, position)]
)

system.time(
  gwas.long.dp <- gwas.lmi.long %>%
    filter(Chr == 1) %>%
    group_by(Chr, position) %>%
    summarise(
      trait = paste0(unique(Trait, collapse = " ")),
      subpop = paste0(unique(Group, collapse = " "))
    )
)

gwas.lmi.long[Chr == 1, length(unique(position))]

lmi.regions <- gwas.lmi.long[, .(
  lmi.trait = paste0(unique(Trait), collapse = " "),
  lmi.subpop = paste0(unique(Group), collapse = " "),
  lmi.year = paste0(unique(Place), collapse = " "),
  lmi.bin = paste0(unique(bin.name), collapse = " ")
), by = .(Chr, position)]

rutils::GenerateMessage("Grouping LMI regions by lmi.bin")
lmi.regions[, lmi.qtl :=
              paste0("Chr", Chr, ":", min(position), "-", max(position),
                     "|", unique(lmi.trait), "|", unique(lmi.subpop), "|",
                     unique(lmi.year)),
            by = .(lmi.bin)]

lmi.regions[, lmi.bin := NULL]
setnames(lmi.regions, "Chr", "chr.num")

rutils::GenerateMessage("Sorting regions: LMI")
setkey(lmi.regions, chr.num, position)

# long qtls for crowell
rutils::GenerateMessage("Generating regions: Crowell")
crowell.regions.wide <- gwas.crowell[, .(
  trait.collapsed = paste0(unique(TRAIT), collapse = " "),
  subpop.collapsed = paste0(unique(SUB_POP), collapse = " ")),
  by = .(CHR, region_start, region_end, Bin_id)]

crowell.regions <- crowell.regions.wide[, .(
  position = c(region_start:region_end),
  trait.collapsed, subpop.collapsed, Bin_id),
  by = .(CHR, region_start, region_end)]

crowell.regions <- crowell.regions[, .(
  chr.num = CHR, position,
  crowell.qtl = Bin_id, crowell.trait = trait.collapsed,
  crowell.subpop = subpop.collapsed
  )]

rutils::GenerateMessage("Sorting regions: Crowell")
setkey(crowell.regions, chr.num, position)

# merge
rutils::GenerateMessage("Merging regions")
merge.both <- merge(lmi.regions, crowell.regions,
      by = c("chr.num", "position"),
      all = TRUE)
rutils::GenerateMessage("Sorting merged regions")
setkey(merge.both, chr.num, position)

# collapse regions
rutils::GenerateMessage("Crowell-specific regions")
only.crowell <- merge.both[!is.na(crowell.trait) & is.na(lmi.trait), .(
  start = min(position), end = max(position), overlap = FALSE,
  lmi.trait = NA, lmi.subpop = NA, lmi.year = NA, lmi.qtl = NA
), by = .(chr.num, crowell.trait, crowell.subpop, crowell.qtl)]
rutils::GenerateMessage("LMI-specific regions")
only.lmi <- merge.both[is.na(crowell.trait) & !is.na(lmi.trait), .(
  start = min(position), end = max(position), overlap = FALSE,
  crowell.trait = NA, crowell.subpop = NA, crowell.qtl = NA
), by = .(chr.num, lmi.trait, lmi.subpop, lmi.year, lmi.qtl)]
rutils::GenerateMessage("Overlapping regions")
overlap <- merge.both[!is.na(crowell.trait) & !is.na(lmi.trait), .(
  start = min(position), end = max(position), overlap = TRUE
), by = .(
  chr.num, lmi.trait, lmi.subpop, lmi.year, lmi.qtl, crowell.trait,
  crowell.subpop, crowell.qtl
)]

rutils::GenerateMessage("Combining regions")
regions <- rbindlist(list(overlap, only.lmi, only.crowell),
                     use.names = TRUE)

rutils::GenerateMessage("Sorting regions")
setkey(regions, chr.num, start, end)

setcolorder(regions, c(
  "chr.num", "start", "end", "overlap", "lmi.qtl", "lmi.trait", "lmi.subpop",
  "lmi.year", "crowell.qtl", "crowell.trait", "crowell.subpop"
))

# output
rutils::GenerateMessage("Writing")
write.table(regions, file = "output/overlap_regions.tsv", sep = "\t",
            quote = FALSE, na = "", row.names = FALSE)

rutils::GenerateMessage("Done")
quit(save = "no",
     status = 0)

# get unique bins
setkey(gwas.crowell, CHR, region_start, region_end)
gwas.crowell[, region_start := as.integer(region_start)]
gwas.crowell[, region_end := as.integer(region_end)]
pd.test <- gwas.crowell[, .(
  trait.collapsed = paste0(unique(TRAIT), collapse = " "),
  subpop.collapsed = paste0(unique(SUB_POP), collapse = " ")
), by = .(CHR, region_start, region_end)]

setkey(gwas.lmi, Chr, Start, End)
pd.test2 <- gwas.lmi[, .(
  trait.collapsed = paste0(unique(Trait), collapse = " "),
  subpop.collapsed = paste0(unique(Group), collapse = " ")
), by = .(Chr, Start, End)]


pd.test[CHR == 1, range(region_start)]
pd.test2[Chr == 1, range(Start)]

which.chromosome <- 4

ggplot() +
  geom_rect(data = pd.test[CHR == which.chromosome],
            mapping = aes(xmin = region_start, xmax = region_end,
                          fill = trait.collapsed),
            ymin = 0.25, ymax = 0.75) +
  geom_rect(data = pd.test2[Chr == which.chromosome],
            mapping = aes(xmin = Start, xmax = End,
                          fill = trait.collapsed),
            ymin = 1.25, ymax = 1.75) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0.5, 1.5),
                     labels = c("Crowell", "LMI"))

