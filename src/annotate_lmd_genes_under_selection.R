#!/usr/bin/env Rscript

library(data.table)

# read window genes
annotated.window.genes <- readRDS(
  "data/annotated_window_genes.RDS")

# get tpm by expressed gene for LRT
lrt.tpm <- data.table(readRDS("data/lmd/tpm.Rds"), keep.rownames = TRUE,
                      key = "rn")
setnames(lrt.tpm, "rn", "gene")
expressed.lrt <- readRDS("data/lmd/expressedGenesAll.Rds")
lrt.tpm <- lrt.tpm[expressed.lrt]

# get LRT p-values for the LMD
dds.lrt <- readRDS("data/lmd/ddsLrt.Rds")
lrt.res <- data.table(data.frame(DESeq2::results(dds.lrt)),
                      keep.rownames = TRUE, key = "rn")
setnames(lrt.res, "rn", "gene")
lrt.p <- unique(lrt.res[, .(gene, padj)])

# combine padj and tpm
lrt.tpm.p <- lrt.p[lrt.tpm]

# combine for genes in windows
setkey(lrt.tpm.p, gene)
setkey(annotated.window.genes, gene)

selection.lmd <- lrt.tpm.p[annotated.window.genes]

# remove genes with no expression values
sdc <- c("n1r1", "n1r3", "n1r4", "n2r1", "n2r3", "n2r4", "n3r1", "n3r2",
         "n3r3", "n4r1", "n4r2", "n4r3")
keep <- selection.lmd[, !all(is.na(.SD)), .SDcols = sdc, by = gene][
  V1 == TRUE, unique(gene)]
selection.lmd <- selection.lmd[gene %in% keep]

# re-order columns
setcolorder(selection.lmd, c(
  "Chr", "window.start", "window.end", "test", "gene", "symbol", 
  "padj", "annotation", "OgroObjective", "OgroRef", "n1r1", "n1r3", 
  "n1r4", "n2r1", "n2r3", "n2r4", "n3r1", "n3r2", "n3r3", "n4r1", 
  "n4r2", "n4r3"))

# write output
outdir = "output/genes_under_selection"
if(!dir.exists(outdir)){
  dir.create(outdir)
}
saveRDS(selection.lmd, paste0(outdir, "/window_genes_x_lmd.RDS"))
write.table(selection.lmd,
            paste0(outdir, "/window_genes_x_lmd_all.tab"),
            sep = "\t",
            quote = FALSE,
            na = "",
            row.names = FALSE)
write.table(selection.lmd[padj < 0.1],
            paste0(outdir, "/window_genes_x_lmd_sig.tab"),
            sep = "\t",
            quote = FALSE,
            na = "",
            row.names = FALSE)

sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, paste0(outdir, "/SessionInfo.genes_under_selection_lmd.txt"))




