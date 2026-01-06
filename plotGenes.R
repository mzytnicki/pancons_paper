library(Gviz)
library(rtracklayer)
library(txdbmaker)

args <- commandArgs(trailingOnly = TRUE)
conservation_file_name <- args[[1]]
genes_file_name        <- args[[2]]

genes <- import(genes_file_name)
genes <- genes[genes$type == "gene"]
txdb <- txdbmaker::makeTxDbFromGFF(genes_file_name)

genome_track <- GenomeAxisTrack()
pancons_track <- DataTrack(range = conservation_file_name, name = "score", type = "hist", window = 100)
#pancons_track <- DataTrack(range = conservation_file_name, type = "l", name = "score")
gene_track <- GeneRegionTrack(txdb, name = "annotation", transcriptAnnotation = "symbol")

for (i in seq_along(genes)) {
    gene <- genes[i]
    chr <- seqnames(gene)
#   chromosome(genome_track) <- tolower(chr)
#   chromosome(pancons_track) <- tolower(chr)
#   chromosome(gene_track) <- tolower(chr)
    s <- start(gene) - 1000
    e <- end(gene)   + 1000
    symb <- gene$symbol
    png(paste0("plot_", symb, "_", chr, "_", s, "_", e, ".png"), width = 16, height = 10, units = "cm", res = 300)
    plotTracks(list(genome_track, gene_track, pancons_track), from = s, to = e)
    dev.off()
}
