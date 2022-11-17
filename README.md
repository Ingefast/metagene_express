![This is an image](/images/title.png)

# INTRODUCTION

A quick approach to perform comparative analysis of high throughput sequencing data on selected regions across several samples using metagene plots, boxplots and heatmaps. The pipeline processed coverage data from DNA methylation (Bisulphite Sequencing), histone methylation (ChIP-seq, cutNrun and similar), and small RNA (sRNA) expression experiments. The scripts consider genes and transposable elements using the *Arabidopsis* TAIR10 genome but should be easily adaptable to any genomic region of interest in any other organisms.

The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script header in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions.
 
# SUPPORTED PLATFORMS

Linux/Unix, Mac OS

Scripts were tested using GNU bash (v4.4.20) and R (v3.4.4) in a Ubuntu 18.04 linux system , and in the Terminal utility and R console (v4.1.1) in macOS.

# DEPENDENCIES

The following additional tools should ideally be available in the ``PATH`` environment. The pipeline is fully functional and tested with the indicated versions of the packages listed below. Other versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

[bedtools](https://bedtools.readthedocs.io/en/latest/#) (v2.26.0)

[samtools](http://www.htslib.org/) (v1.3.1)


# INSTALLATION

The scripts can be cloned and run directly in R.
```
git clone https://github.com/Ingefast/metagene_express
cd metagene_express
```
# WORKFLOW

![This is an image](/images/flowchart.png)


## 1. Creation of framework files for a particular genome

The script **`metagene.frame_builder.r`** delineates the binning of the genomic features to be analysed. Each feature is defined by sixty windows or bins in total: a 2kb long flank made of twenty 100bp long bins upstream the Transcription Start Sites (TSS), forty equally long bins along the gene or transposon coding body, and a second 2kb long flank of twenty 100bp long bins downstream the Transcription End Site (TES).

The input to be used is a file with the annotated features in bed format (**`TAIR10_genes_transposons.bed`**). Below an example of how to prepare such a file for the *Arabidopis* TAIR10 genome (chloroplast and mitochondrion are discarded):
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

The script **`metagene.matrixmaker.r`** creates a matrix with genomic features as rows and average coverage values for each bin as columns.

Input should be a bed file (**`track.sorted.bed`**) with normalised coverage values for ChiP-seq, cutNrun, small-RNA or similar datasets. Looking like e.g.
```
Chr1	725	850	3.58513
Chr1	920	1045	3.58513
Chr1	1100	1225	3.58513
Chr1	1303	1428	3.58513
Chr1	1685	1739	3.58513
```

If not available as such, it can easily be prepared out of sorted alignment file (**`bam`**) in a linux environment with:
```
	librarySize=$(samtools idxstats out1.dedup.sorted.bam | awk '{total+=$3}END{print total}');
	b=1000000;
	factor=`echo "$b / $librarySize "|bc -l`;
	bedtools genomecov -ibam out1.dedup.sorted.bam -g TAIR10.chrom.notPt.notMt.sizes -bg -scale $factor |grep -v "chloroplast\|mitochondria"> track.sorted.bed;
```

The **`TAIR10.chrom.notPt.notMt.sizes`** referred above and also needed for **`metagene.matrixmaker.r`** is just a plain text file with chromosome sizes. It is provided under [example](https://github.com/Ingefast/metagene_express/tree/main/example):
```
Chr1    30427671
Chr2    19698289
Chr3    23459830
Chr4    18585056
Chr5    26975502
```

Read the [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) instruction website for help if problems arise.

When having DNA methylation data the script **`metagene.cx.track_maker.r`** should be used to preprocess the cytosine methylation reports from the **`bismark_methylation_extractor`** from [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) separately in each cytosine context (CG, CHG, CHH).

## 3. Looking at the data: metagene plots, boxplots and heatmaps.

The script **`metagene.plotter.r`** creates a combined line graph linking the average coverage or methylation values across each feature and its flanks, and a boxplot for the values only within the feature body. In the example below four tracks representing two conditions with two replicates are presented.


![This is an image](/images/figure1.png)
*Figure 1*. Example of metagene plots (left) and boxplots (right) of DNA methylation in CG context for genes and TEs in *Arabidopsis*. Yellow asteriscs in boxes represent mean values. It has to be beared in mind that the metagene plot line represent mean values, while the boxplots notched midline represent medians.


The same **`track.matrix.txt`** files produced at step 2 can also be used to build heatmap as in figure 2 using **`metagene.heatmap.r`**.


![This is an image](/images/figure2.png)
*Figure 2*. Heatmap of a single **`track.matrix.txt`** file representing DNA methylation in CG context for genes and TEs in *Arabidopsis* ordered after decreasing methylation values. Grey cells correspond to bins without coverage (NA data).


# REFERENCES

Modified versions of this scripts have been used to process the the datasets in the following papers:

1. Schatlowski N et al (2014). Hypomethylated pollen bypasses the interploidy hybridization barrier in *Arabidopsis*. **Plant Cell** 26 (9) 3556-68.
2. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.
3. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.santos at slu.se
