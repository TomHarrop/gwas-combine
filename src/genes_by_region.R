#!/usr/bin/env Rscript
library(GenomicRanges)
library(data.table)

# load disjoint gwas regions
rutils::GenerateMessage("Loading disjoint QTLs")
gwas.combined <- readRDS("output/combined_gr.Rds")

# load gtf
rutils::GenerateMessage("Loading gtf")
gtf.file <- "data/Osativa_323_v7.0.gene_exons.gffread.rRNAremoved.gtf"
gtf <- rtracklayer::import.gff(gtf.file, format = 'gtf',
                               genome = 'Osativa_323_v7',
                               feature.type="exon")

# split gtf by gene
rutils::GenerateMessage("Splitting gtf by gene_name")
gtf.split <- split(gtf, elementMetadata(gtf)$gene_name)
grl <- reduce(gtf.split)

# find genes in regions
rutils::GenerateMessage("Finding overlaps")
region.genes <- as.data.table(findOverlaps(gwas.combined, grl))

# add region info
rutils::GenerateMessage("Getting metadata for overlaps")
region.info <- region.genes[, as.list(mcols(gwas.combined[queryHits])),
                            by = queryHits]
region.info[, start := start(gwas.combined[queryHits]), by = queryHits]
region.info[, end := end(gwas.combined[queryHits]), by = queryHits]                             
region.info[, Chr := as.character(seqnames(gwas.combined[queryHits])),
            by = queryHits]

rutils::GenerateMessage("Getting gene_names for overlaps")
rutils::PrintF("Number of genes: %s\n",
               region.genes[, length(unique(subjectHits))])
region.genes[, gene := names(grl[subjectHits]), by = subjectHits]

# merge back together
rutils::GenerateMessage("Merging region info with gene names")
setkey(region.info, queryHits)
setkey(region.genes, queryHits)
gwas.genes <- merge(region.info, region.genes, by = "queryHits")

# tidy columns
gwas.genes <- gwas.genes[, .(
  gene, Chr, start, end, lmi.trait, lmi.place, lmi.group, crowell.trait,
  crowell.group
)]
setkey(gwas.genes)
gwas.genes <- unique(gwas.genes)

# write output
rutils::GenerateMessage("Saving output")
saveRDS(gwas.genes, "output/regions_with_genes.Rds")
sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, "output/genes_by_region.SessionInfo.txt")
