library(data.table)
library(GenomicRanges)

#############
# FUNCTIONS #
#############

CollapseMetadata <- function(y) {
  # collapse rows of mcols() data (metadata) by pasting unique values
  y.n <- names(y)
  out <- lapply(y.n, function(z) {
    records <- as.character(unique(y[[z, ]]))
    records <- records[!is.na(records) & records != ""]
    records <- unique(unlist(sapply(records, strsplit, split = "[[:blank:]]")))
    paste(records, collapse = " ")
  })
  names(out) <- y.n
  return(out)
}

# x <- comb
# y <- mc.orig[c(14L, 220L, 221L, 222L, 223L, 224L, 225L, 226L, 227L), ]

ReduceWithMetadata <- function(x) {
  # original metadata
  mc.orig <- mcols(x)
  
  # reduce ranges and keep original mappings
  x.red <- GenomicRanges::reduce(x, with.revmap = TRUE)
  
  # generate mappings from x.red to x
  mc.new <- data.table(i.new = 1:length(x.red))
  mc.new <- mc.new[, .(
    i.orig = mcols(x.red)$revmap[[i.new]]
  ), by = i.new]
  
  # collapse metadata by row in x.red
  mc.new <- mc.new[, CollapseMetadata(mc.orig[i.orig,]), by = i.new]
  
  # add metadata, excluding i
  mcols(x.red) <- mc.new[, setdiff(names(mc.new), "i.new"), with= FALSE]
  
  return(x.red)
}

########
# DATA #
########

# load seqlengths
seqlengths.file <- "data/Osativa_323_v7_chrom_length.tab"
seqlengths.df <- read.table(seqlengths.file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE)

# make seqinfo
seqinf <- Seqinfo(seqnames = seqlengths.df$V1, seqlengths = seqlengths.df$V2,
                  genome = "Osativa_323_v7")

# load lmi data
gwas.lmi.file <- "data/gwas_lmi_genes.Rds"
gwas.lmi <- readRDS(gwas.lmi.file)

# load crowell data
gwas.crowell.file <- "data/gwas_crowell_genes.Rds"
gwas.crowell <- readRDS(gwas.crowell.file)

####################################
# 1. `reduce` OVERLAPS IN LMI DATA #
####################################

# format lmi data for GRanges
gwas.lmi <- gwas.lmi[, .(
  Chr, Start, End,
  lmi.trait = Trait, lmi.place = Place, lmi.group = Group,
  crowell.trait = NA, crowell.group = NA)]
gwas.lmi <- gwas.lmi[!is.na(End)]
setkey(gwas.lmi)
gwas.lmi <- unique(gwas.lmi)
gwas.lmi[, chr.name := paste0("Chr", Chr)]
gwas.lmi[, Chr := NULL]

# make GRanges object
lmi.gr <- makeGRangesFromDataFrame(
  df = gwas.lmi,
  start.field = "Start",
  end.field = "End",
  seqnames.field = "chr.name",
  ignore.strand = TRUE,
  keep.extra.columns = TRUE,
  seqinfo = seqinf
)

# find overlaps in lmi ranges
lmi.red <- ReduceWithMetadata(lmi.gr)

###########################
# 2. PREPARE CROWELL DATA #
###########################

# format for GRanges
gwas.crowell <- gwas.crowell[, .(
  region_name,
  region_start = as.numeric(region_start),
  region_end = as.numeric(region_end),
  lmi.trait = NA, lmi.place = NA, lmi.group = NA,
  crowell.trait = TRAIT, crowell.group = SUB_POP
)]
setkey(gwas.crowell)
gwas.crowell <- unique(gwas.crowell)

# make GRanges object
crowell.gr <- makeGRangesFromDataFrame(
  df = gwas.crowell,
  start.field = "region_start",
  end.field = "region_end",
  seqnames.field = "region_name",
  ignore.strand = TRUE,
  keep.extra.columns = TRUE,
  seqinfo = seqinf
)

# reduce crowell ranges
crowell.red <- ReduceWithMetadata(crowell.gr)

#####################################
# 3. OVERLAP LMI AND CROWELL RANGES #
#####################################

# combine regions
comb <- c(crowell.red, lmi.red)

# reduce again
combined.red <- ReduceWithMetadata(comb)

###################
# 4. WRITE OUTPUT #
###################

combined.dt <- as.data.table(combined.red)
combined.dt[, c("width", "strand") := NULL]
setnames(combined.dt, "seqnames", "Chr")

write.table(combined.dt, "output/condensed_qtls.csv", quote = FALSE, sep = ",",
            na = "", row.names = FALSE)

