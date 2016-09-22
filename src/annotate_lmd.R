#!/usr/bin/env Rscript

library(data.table)

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

# merge tpm with gwas results
disjoint.genes <- readRDS("output/regions_with_genes.Rds")
setkey(disjoint.genes, gene)
setkey(lrt.tpm.p, gene)

gwas.lmd <- lrt.tpm.p[disjoint.genes]

# add annotations
annot <- gwas.lmd[, data.table(oryzr::LocToGeneName(unique(gene)),
                              keep.rownames = TRUE, key = "rn")]
setnames(annot, "rn", "gene")
gwas.lmd.annot <- annot[gwas.lmd]

# remove genes with no expression values
sdc <- c("n1r1", "n1r3", "n1r4", "n2r1", "n2r3", "n2r4", "n3r1", "n3r2",
         "n3r3", "n4r1", "n4r2", "n4r3")
gwas.lmd.annot[, chuck := all(is.na(.SD)), .SDcols = sdc, by = gene]
gwas.lmd.annot <- gwas.lmd.annot[!chuck == TRUE]
gwas.lmd.annot[, chuck := NULL]

# make chromosome numeric
gwas.lmd.annot[, chr.num := as.numeric(gsub("Chr", "", Chr)), by = Chr]
gwas.lmd.annot[, Chr := NULL]
setnames(gwas.lmd.annot, "chr.num", "Chr")

# re-order columns
ord <- c("gene", "symbols", "names", "padj", "MsuAnnotation", sdc, "Chr",
         "start", "end", "lmi.trait", "lmi.place", "lmi.group",
         "crowell.trait", "crowell.group")
others <- setdiff(names(gwas.lmd.annot), ord)
setcolorder(gwas.lmd.annot, c(ord, others))

# why are there tabs in the annotations???
gwas.lmd.annot[, names := gsub("\t", "", names, fixed = TRUE)]

# write out
setkey(gwas.lmd.annot, Chr, start, end)
write.table(gwas.lmd.annot,
            "output/gwas_lmd.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
write.table(gwas.lmd.annot[padj < 0.05],
            "output/gwas_lmd_sig.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
saveRDS(gwas.lmd.annot, "output/gwas_lmd_annot.Rds")
sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, "output/SessionInfo.annotate_lmd.txt")





