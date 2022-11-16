![This is an image](/images/title.png)

# INTRODUCTION

A quick approach to perform comparative analysis of genomic data on selected regions across several samples by means of metagene plots, boxplots and heatmaps. The scripts focus mainly on coverage data from DNA methylation (Bisulphite Sequencing), histone methylation (ChIP-seq, cutNrun and similar), and small RNA (sRNA) expression experiments. The scripts consider genes and transposable elements using  the *Arabidopsis* TAIR10 genome but should be easily adaptable to any genomic region of interest in any other organisms.

The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script header in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions.
 
# SUPPORTED PLATFORMS

Linux/Unix, Mac OS

Scripts were tested using GNU bash (v4.4.20) and R (v3.4.4) in a Ubuntu 18.04 linux system , and the Terminal utility and the R console (v4.1.1) in macOS.

# DEPENDENCIES

The following additional tools should ideally be available in the PATH environment. The pipeline is fully functional and tested with the indicated versions of the packages listed below. Other versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

[bedtools](https://bedtools.readthedocs.io/en/latest/#) (v2.26.0)

[samtools](http://www.htslib.org/) (v1.3.1)


# INSTALLATION

The scripts can be cloned and run directly in R.
```
git clone https://github.com/Ingefast/metagene_express
cd metagene_express
```

# WORKFLOW

## 1. Creation of framework files for a particular genome

The script **`metagene.frame_builder.r`** delineates the binning of the genomic features to be analysed. Each feature is defined by sixty windows or bins in total: a 2kb long flank made of twenty 100bp long bins upstream the Transcription Start Sites (TSS), forty equally long bins along the gene or transposon coding body, and a second 2kb long flank of twenty 100bp long bins downstream the Transcription End Site (TES).

The input to be used is a file with the annotated features in bed format. Below an example of how to prepare such a file for TAIR10 genome:
```
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff

#selects only nuclear genes and transposable elements as the features to analyse
grep -w 'gene\|transposable_element' TAIR10_GFF3_genes_transposons.gff | grep -wv 'Pt\|Mt'>TAIR10_GFF3_gene_TE.ONLY_chr.gff

#transforms gff3 to bed format
gff2bed < TAIR10_GFF3_gene_TE.ONLY_chr.gff > TAIR10_genes_transposons.bed
```

The script creates and output made of six bed files for the whole genome looking like:
```
Chr1	3630	3687	AT1G01010	1	+
Chr1	3687	3744	AT1G01010	2	+
Chr1	3744	3801	AT1G01010	3	+
Chr1	3801	3858	AT1G01010	4	+
Chr1	3858	3915	AT1G01010	5	+
```

## 2. Creation of a matrix file with binned coverage for every feature.

The script **`metagene.matrixmaker.r`** creates a matrix with genomic features as rows and mean coverage values for each bin as columns.

Input is a bed file (**`track.sorted.bed`**) with normalised coverage values for ChiP-seq, small-RNA of RNA-seq datasets. In a linux environment can easily be prepared out of sorted alignment file (**`bam`**) with:
```
	librarySize=$(samtools idxstats out1.dedup.sorted.bam | awk '{total+=$3}END{print total}');
	b=1000000;
	factor=`echo "$b / $librarySize "|bc -l`;
	bedtools genomecov -ibam out1.dedup.sorted.bam -g ~/genome_ref/TAIR10/TAIR10.chrom.sizes -bg -scale $factor |grep -v "chloroplast\|mitochondria"> track.sorted.bed&
```
When having DNA methylation data the script **`cx.report.track_maker.r`** should be used to preprocess the cytosine methylation reports from [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) in each cytosine context. 

## 3. Looking at the data: plotting metageneplots, boxplots and heatmaps.

The script **`metagene.plotter.r`** creates a combined line graph linkig mean values across bins and a boxplot for the genomic values in the feature body. In the example below four tracks representing two conditions with two replicates are presented.

![This is an image](/images/figure1.png)

*Figure 1*. (A) bedGraph files of a wildtype in *Capsella* for different sRNA sizes. (B) Correlogram of 24nt sRNA values over genes in three conditions with two replicates each. (C) NMDS diagram of the same dataset.

The **`track.matrix.txt`** files produced at step 2 can also be used to build heatmap as in figure 2.

*Figure 2*


# REFERENCES

Modified versions of this scripts have been used to process the the datasets in the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.santos at slu.se
