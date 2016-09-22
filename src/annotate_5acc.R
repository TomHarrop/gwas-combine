#!/usr/bin/env Rscript

library(data.table)

# load l2fc data
rutils::GenerateMessage("Loading DE results")
stage.l2fc.dom.padj <- readRDS("data/fa/stage_l2fc_dom_padj.Rds")

# cast l2fc data (1 row per gene)
l2fc.wide <- dcast(stage.l2fc.dom.padj, gene ~ accession, value.var = c(
  "log2FoldChange", "lfcSE", "padj"
))

# join wide dts
dom.data <- unique(stage.l2fc.dom.padj[, .(
  gene, dom_all, dom_asia, dom_africa, dom_japonica, dom_indica
)])

# load tpm data
rutils::GenerateMessage("Load TPM gene expression data")
tpm.long <- readRDS("data/fa/tpm_with_calls.Rds")
ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
         "O. barthii", "O. glaberrima")
tpm.long[ , accession := factor(plyr::mapvalues(
  species,
  from = c("R", "I", "J", "B", "G"),
  to = ord),  levels = ord)
  ]
tpm.basemean <- tpm.long[, .(mean.tpm = mean(tpm)),
                         by = .(gene, accession, stage)]
tpm.wide <- dcast(tpm.basemean, gene ~ accession + stage,
                  value.var = "mean.tpm")

# merge tpm with dom.l2fc
dom.tpm.l2fc <- merge(merge(dom.data, tpm.wide, "gene"), l2fc.wide, "gene")

# subset by gwas results
disjoint.genes <- readRDS("output/regions_with_genes.Rds")
setkey(disjoint.genes, gene)
setkey(dom.tpm.l2fc, gene)

gwas.5acc <- dom.tpm.l2fc[disjoint.genes]

# add annotations
annot <- gwas.5acc[, data.table(oryzr::LocToGeneName(unique(gene)),
                               keep.rownames = TRUE, key = "rn")]
setnames(annot, "rn", "gene")
gwas.5acc.annot <- annot[gwas.5acc]

# remove genes with no expression values
sdc <- grep("[PBM|SM]$", names(gwas.5acc.annot), value = TRUE)
gwas.5acc.annot[, chuck := all(is.na(.SD)), .SDcols = sdc, by = gene]
gwas.5acc.annot <- gwas.5acc.annot[!chuck == TRUE]
gwas.5acc.annot[, chuck := NULL]

# make chromosome numeric
gwas.5acc.annot[, chr.num := as.numeric(gsub("Chr", "", Chr)), by = Chr]
gwas.5acc.annot[, Chr := NULL]
setnames(gwas.5acc.annot, "chr.num", "Chr")

# set column order
ord <- c("gene", "symbols", "names", "MsuAnnotation",
         "dom_all",  "dom_asia", "dom_africa", "dom_japonica", "dom_indica",
         sdc,
         "Chr", "start", "end", "lmi.trait", "lmi.place", "lmi.group",
         "crowell.trait", "crowell.group")
others <- setdiff(names(gwas.5acc.annot), ord)
setcolorder(gwas.5acc.annot, c(ord, others))

# choose significant genes
sdc <- grep("^dom", names(gwas.5acc.annot), value = TRUE)
sg <- gwas.5acc.annot[, any(.SD < 0.05),
                      .SDcols = sdc, by = gene][
  V1 == TRUE,unique(gene)]

# write out
setkey(gwas.5acc.annot, Chr, start, end)
write.table(gwas.5acc.annot,
            "output/gwas_5acc.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
write.table(gwas.5acc.annot[gene %in% sg],
            "output/gwas_5acc_sig.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
saveRDS(gwas.5acc.annot, "output/gwas_5acc_annot.Rds")
sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, "output/SessionInfo.annotate_5acc.txt")
