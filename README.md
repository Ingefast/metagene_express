![This is an image](/images/title.png)

# INTRODUCTION

A quick approach to perform comparative analysis of genomic data on selected regions across several samples by means of metagene plots, boxplots and heatmaps. The scripts focus mainly on coverage data from DNA methylation (Bisulphite Sequencing), histone methylation (ChIP-seq, cutNrun and similar), and small RNA (sRNA) expression experiments. The scripts consider genes and transposable elements but are easily adaptable to any genomic region of interest.

It is presented here using the *Arabidopsis* TAIR10 genome as example. The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions and hints (worth to read as a complement to this instruction). 
 
# SUPPORTED PLATFORMS

Linux/Unix, Mac OS

Shell scripts (*.sh) of this software were developed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system and the Terminal utility in macOS. R scripts were tested using the R console (v4.1.1) in macOS.

# DEPENDENCIES

The following tools need to be installed and ideally available in the PATH environment. The pipeline is fully functional and tested with the following versions of the packages listed below. Other versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

bedtools (v2.26.0)

samtools (v1.3.1)


# INSTALLATION

Shell scripts can be cloned and run directly on a linux server.

```
git clone https://github.com/Ingefast/metagene_express
cd metagene_express
```

# WORKFLOW

## 1. **`metagene.frame_builder.r`** Creation of framework files for a particular genome

The input of this script consist of an genomic annotation in bed format for the relevant features. Below follows an example of how to prepare such an input file if not available.

```
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff

#selects only nuclear genes and transposable elements as the features to analyse
#grep -w 'gene\|transposable_element' TAIR10_GFF3_genes_transposons.gff | grep -wv 'Pt\|Mt'>TAIR10_GFF3_gene_TE.ONLY_chr.gff

#transforms gff3 to bed format
#gff2bed < TAIR10_GFF3_gene_TE.ONLY_chr.gff > TAIR10_genes_transposons.bed
```
It creates a six bed files for the whole genome looking like:

```
Chr1	3630	3687	AT1G01010	1	+
Chr1	3687	3744	AT1G01010	2	+
Chr1	3744	3801	AT1G01010	3	+
Chr1	3801	3858	AT1G01010	4	+
Chr1	3858	3915	AT1G01010	5	+
```

## 2. **`metagene.matrix_maker_.r`** creates a matrix with binned coverage for every feature.

Creates a matrix with genomic features as rows and coverage in bins as columns.


Input is a bed file with normalised coverage for ChiP-seq, small-RNA of RNA-seq datasets. It can easily be prepared out of sorted BAM file with:

```
	librarySize=$(samtools idxstats out1.dedup.sorted.bam | awk '{total+=$3}END{print total}');
	b=10000000;
	factor=`echo "$b / $librarySize "|bc -l`;
	bedtools genomecov -i tmp.specific.ext.bed -g ~/genome_ref/TAIR10/TAIR10.chrom.sizes -d -scale $factor > out1.coverage.txt&
```
When having DNA methylation data
**`cx.report.track_maker.r`** should be used to create the input **`track.sorted.bed`***.

The output file **`track.matrix.txt`**.


## 3. Looking at the data: plotting metageneplots, boxplots and heatmaps.

1. **`metagene.plotter.r`** creates a combine line graph representing means and a boxplot.

Transcription Start Sites (TSS) and Transcription End Site (TES) for each gene



![This is an image](/images/figure1.png)

*Figure 1*. (A) bedGraph files of a wildtype in *Capsella* for different sRNA sizes. (B) Correlogram of 24nt sRNA values over genes in three conditions with two replicates each. (C) NMDS diagram of the same dataset.

![This is an image](/images/figure2.png)

*Figure 2*. sRNA size distribution over genes and transposable elements in *Arabidopsis*.



# REFERENCES

Modified versions of this scripts have been used to process the the datasets in the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.santos at slu.se
