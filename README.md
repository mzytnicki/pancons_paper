# Analysis for the `pancons` paper

## Collect data

### Pangenome graphs

*Arabidopsis thaliana* genomes were downloaded using data from [NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1033522/), and the pangenome graph was built using [Pan1c](https://forge.inrae.fr/genotoul-bioinfo/Pan1c/pan1c).

Human pangenome graph were downloaded with the following code

```
for i in $( seq 22 )
do
  wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2022_03_11_minigraph_cactus/chrom-graphs-hprc-v1.1-mc-chm13-full/chr${i}.vg
  vg convert -fW chr${i}.vg | gzip -c > chr${i}.gfa.gz
  rm chr${i}.vg
done
```

### Annotation

### *A. thaliana*

Download from [TAIR](https://www.arabidopsis.org/download/file?path=Genes%2FAraport11_genome_release%2FAraport11_GFF3_genes_transposons.current.gff.gz)

Lifton from Col-0 to Col-PEK (using sequences downloaded from the *Arabidopsis* publication).
```
lifton -g Araport11_GFF3_genes_transposons.20250813.gff.gz -o Araport11_GFF3_genes_transposons.20250813.lifton.gff -t 10 Col-0.fasta TAIR9.fasta
```

Collect features:
```
zcat Araport11_GFF3_genes_transposons.20250813.lifton.gff | awk '$3 == "CDS"' | xz -c > Arabidopsis_thaliana_genes.coding.gff.xz
zcat Araport11_GFF3_genes_transposons.20250813.lifton.gff | awk '$3 == "gene"' | xz -c > Arabidopsis_thaliana_genes.gene.gff.xz
zcat Araport11_GFF3_genes_transposons.20250813.lifton.gff | awk '$3 == "exon"' | xz -c > Arabidopsis_thaliana_genes.exon.gff.xz
xzcat Arabidopsis_thaliana_TEs.gff.xz | awk '$3 == "Gypsy_LTR_retrotransposon"' | xz -c > Arabidopsis_thaliana_TEs.gypsy.gff.xz

for i in genes.coding genes.gene genes.exon TEs.gypsy
do
    xzcat Arabidopsis_thaliana_${i}.gff.xz | awk '{ print $1 "\t" ($4 - 1) "\t" $5 "\tregion" NR "\t.\t+"}' | sort -k1,1 -k2,2n | bedtools merge -i - | gzip -c > Arabidopsis_thaliana_${i}.bed.gz
done

bedops --difference <( zcat Arabidopsis_thaliana_genes.gene.bed.gz ) <( zcat Arabidopsis_thaliana_genes.exon.bed.gz ) | gzip -c > Arabidopsis_thaliana_genes.intron.bed.gz
bedops --difference <( zcat Arabidopsis_thaliana_genes.exon.bed.gz ) <( zcat Arabidopsis_thaliana_genes.coding.bed.gz ) | gzip -c > Arabidopsis_thaliana_genes.UTR.bed.gz

python3 computeFrames.py <( xzcat Arabidopsis_thaliana_genes.coding.gff.xz ) Arabidopsis_thaliana_genes.frame
for i in Arabidopsis_thaliana_genes.frame[0-2].bed
do
    gzip $i
done
```

### Human

Human annotations: [GenCode](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz).

## Running `pancons`

### *A. thaliana*

```
for i in $(seq 5)
do
    /usr/bin/time ./pancons -i <( zcat arabidopsis.gfa.gz ) -r Col_0 -c chr${i} -p 1 > pancons_athaliana_chr${i}.bed 2> pancons_athaliana_chr${i}.log &
done
cat pancons_athaliana_chr[1-5].bed > pancons_athaliana_all.bed
```


## Collecting results

### *A. thaliana*

Plot/fit conservation score distribution.

```
Rscript plotConservation.R pancons_athaliana_all.bed 5 pancons_athaliana_all_conservation.png
Rscript fitConservation.R
```

Display cumulative conservation distribution of each feature.

```
for i in genes.UTR genes.intron genes.frame0 genes.frame1 genes.frame2 TEs.gypsy
do
    for chr in $(seq 5)
    do
        python3 computeScores.py <( zcat Arabidopsis_thaliana_${i}.bed.gz ) pancons_athaliana_chr${chr}.bed Chr${chr}
    done > pancons_A_thaliana_${i}.txt &
done

Rscript compareScores.R pancons_A_thaliana_*.txt pancons_A_thaliana_all.png 0.7
```

Plot gene-centered distribution.

```
tmpFileName=tmp.bedgraph
samtools faidx Col-0.fasta
echo "track type=bedGraph autoScale=on alwaysZero=off" > $tmpFileName
awk -v i="$i" '{printf "%s\t%d\t%d\t%.0f\n", $1, $2, $3, $5 * 1000}' pancons_athaliana_all.bed | sed 's/chr/Chr/g' >> $tmpFileName
bedGraphToBigWig $tmpFileName Col-0.fasta.fai pancons_athaliana_all.bw
rm -rf $tmpFileName

computeMatrix scale-regions -S pancons_athaliana_all.bw -R Arabidopsis_thaliana_genes.gene.chr.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrix.mat.gz

plotProfile -m matrix.mat.gz -out profile.png
```

Extract divergent regions

```
python3 getDivergent.py -i pancons_athaliana_all.bed -c 5 > pancons_athaliana_all_div.bed
```

Plot the number of features that match divergent elements.

```
library(ChIPpeakAnno)
txdb <- txdbmaker::makeTxDbFromGFF("Araport11_GFF3_genes_transposons.20250813.lifton.gff.gz")
regions <- read.delim("pancons_athaliana_all_div.bed", header = FALSE)
regions <- regions[, 1:3]
colnames(regions) <- c("chr", "start", "end")
regions$chr <- gsub("chr", "Chr", regions$chr)
regions <- makeGRangesFromDataFrame(regions)
aCR <- assignChromosomeRegion(regions, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"), TxDb=txdb)
aCR <- assignChromosomeRegion(regions, nucleotideLevel=FALSE, precedence=c("Exons", "fiveUTRs", "threeUTRs", "Introns", "Promoters", "immediateDownstream"), TxDb=txdb)
df <- data.frame(n = aCR$jaccard) |> 
  tibble::rownames_to_column(var = "annotation") |>
  dplyr::mutate(annotation = dplyr::case_match(annotation, "Exons" ~ "exon", "Introns" ~ "intron", "fiveUTRs" ~ "5' UTR", "threeUTRs" ~ "3' UTR", "immediateDownstream" ~ "downstream", "Promoters" ~ "upstream", "Intergenic.Region" ~ "intergenic")) |>
  dplyr::mutate(annotation = factor(annotation, levels = c("exon", "5' UTR", "3' UTR", "intron", "upstream", "downstream", "intergenic")))
p <- ggplot2::ggplot(data = df, ggplot2::aes(x = annotation, y = n)) + 
  ggplot2::geom_bar(stat="identity") + 
  ggplot2::scale_y_continuous(labels = scales::percent) + 
  ggplot2::xlab("") + 
  ggplot2::ylab("") + 
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
  ggplot2::theme_minimal()
ggplot2::ggsave("annotations.png", p, width = 10, height = 5, units = "cm")
```

Plot conservation around the *FLC* gene.

```
zgrep "ID=AT5G65050\|ID=AT5G10140" Araport11_GFF3_genes_transposons.20250813.lifton.gff.gz > Araport11_GFF3_selected_genes.gff
Rscript plotGenes.R pancons_athaliana_all.bedgraph Araport11_GFF3_selected_genes.gff
```

Extract the region around it

```
vg find -E -p Col_0#1#5#0:3177218-3183284 -x mc.xg | vg mod -u - | vg convert -f -W - > ATCol0-FLC.gfa
```
