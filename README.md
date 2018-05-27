# sRNAPipe
a GALAXY-based pipeline for bioinformatic in-depth exploration of small RNA-seq data

Description
=============

The field of small RNA is one of the most investigated research areas since they were shown to regulate gene expression and play essential roles in fundamental biological processes.
sRNAPipe  a computational pipeline (sRNAPipe: small RNA pipeline) based on the Galaxy framework that takes as input a fastq file of small RNA-seq reads and performs successive steps of mapping to categories of genomic sequences: microRNAs, gene transcripts, small nuclear RNAs, ribosomal RNAs, transfer RNAs and transposable elements. It also provides individual mapping and counting for chromosomes, gene transcripts and transposable elements, normalization, small RNA length analysis and plotting of the data along genomic coordinates to build publication-quality graphs and figures. sRNAPipe evaluates 10-nucleotide 5’-overlaps of reads on opposite strands to test ping-pong amplification for putative PIWI-interacting RNAs, providing numbers of overlaps and corresponding z-scores.



Prerequisites
=============

1. Unix system with A Galaxy server (release july 2014 or later installed)

2. Some tools are used by sRNAPipe and must be installed and added to the Path.

	* Bwa aligner: you can obtain it here: https://sourceforge.net/projects/bio-bwa/files/ . Please download version  0.7.12-r1039 or higher
	* BedTools powerful toolset for genome arithmetic is also needed. It should be found here: http://bedtools.readthedocs.io/en/latest/ . We recommend to use v2.24.0 or higher.
	* Samtools : you can obtain it here :  https://sourceforge.net/projects/samtools/files/samtools/1.2/. Please download version 1.2

3. Perl version higher than 5.1 is needed with packages : "perl-statistics", "Parallel::ForkManager", "Statistics::R", "Getopt::Long" , "String::Random" installed.

4. R project version higher than 3.1 is needed with libraries "plotrix", "bioconductor-sushi", "RColorBrewer" and "ggplot2"  installed. You can find respectively these libraries here: https://cran.r-project.org/web/packages/plotrix/index.html and https://bioconductor.org/packages/release/bioc/html/Sushi.html and https://cran.r-project.org/web/packages/RColorBrewer/index.html and https://cran.r-project.org/web/packages/ggplot2/index.html

Installation
=============

The process has to be completed by an administrator of your Galaxy server to install sRNAPipe.

1. Download sRNAPipe
You can find sRNAPipe here: https://github.com/brassetjensen/sRNAPipe

2. Put the tool into Galaxy's tools directory
You need to add files into tools/ directory , where all tool-related files are stored, within your Galaxy installation.

3. Make Galaxy aware of the new tool sRNAPipe
Now that the tool and its definition file are ready, the final step is to make Galaxy aware of the new files.
Galaxy recognizes installed tools by reading the tool_conf.xml tool configuration file. Thus, letting Galaxy know about the new tool is as easy as adding a few lines to the tool_conf.xml file located in the config/ directory of the Galaxy installation. New tools can either be added to existing sections or added to new sections defined in the following way:

```xml
 <section name="NewTools" id="mTools">
    <tool file="sRNAPipe.xml" />
 </section>
```

 4. Start or Restart Galaxy to use it.
 
 User Manual
=============
 
 ![User_manual]("https://github.com/brassetjensen/sRNAPipe/blob/master/sRNAPipe_User_Manual.pdf")
