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
sdc <- c("n1r1", "n1r3", "n1r4", "n2r1", "n2r3", "n2r4", "n3r1", "n3r2",
         "n3r3", "n4r1", "n4r2", "n4r3")
gwas.lmd.annot[, chuck := all(is.na(.SD)), .SDcols = sdc, by = gene]
gwas.lmd.annot <- gwas.lmd.annot[!chuck == TRUE]
gwas.lmd.annot[, chuck := NULL]

------

setcolorder(lmi.lmd.annot, c("gene", "symbols", "names", "padj",
                             "MsuAnnotation", "n1r1", "n1r3", "n1r4", "n2r1",
                             "n2r3", "n2r4", "n3r1", "n3r2", "n3r3", "n4r1",
                             "n4r2", "n4r3",
                             "Chr", "Start", "End", "Trait",
                             "Place", "Group", "RapID", "OgroObjective",
                             "OgroRef"))

# why are there tabs in the annotations???
lmi.lmd.annot[, names := gsub("\t", "", names, fixed = TRUE)]

# write out
setkey(lmi.lmd.annot, Chr, Start, End)
write.table(lmi.lmd.annot,
            "output/lmi_lmd_annotated.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
write.table(lmi.lmd.annot[padj < 0.05],
            "output/lmi_lmd_annotated_sig.tab",
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)



