# Exploiting the power of a High Performance Computing cluster for Bioinformatics Analysis 

**Center for Research Informatics, University of Chicago**

Saturdays 04/27/2019; 9:00 AM - 12:00 PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Apply your Linux, shell scripting and HPC skills, to develop an efficient ananlysis pipeline for NGS-based RNA-Seq data analysis

## 1. Introduction

## The Dataset

For this tutorial, we will use the The RNA-seq data from the [paper](https://www.ncbi.nlm.nih.gov/pubmed/25499759), entitled: *Loss of PRDM11 promotes MYC-driven lymphomagenesis* The authors use a functional screening approach to identify novel tumor suppressors among the PR-domain (PRDM) family of genes that encodes transcriptional regulators, several of which are deregulated in cancer. They demonstrate oncogenic collaboration between depletion of the previously uncharacterized PR-domain family member Prdm11 and overexpression of MYC.

We will use the dataset from the PRDM11 knockdown and wildtype samples. The full dataset is available at the Gene Expression Omnibusrepository with ID [GSE56065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56065).

## The pipeline

In this pipeline we will use **FASTQ** tool for rawdata quality control, followed by transcription quantification using a Psuedo-alignment tool [kallisto](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html)

**kallisto** is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, it is based on the idea of **pseudoalignment** for rapidly determining the compatibility of reads with targets, without the need for alignment. 

kallisto can quantify million of RNA-Seq reads in just a few minutes, this is extremely fast, when compared with the time it takes to perform the alignment with mapping tools like [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/). Pseudoalignment of reads preserves the *key* information needed for quantification.

A kallisto manual is available [here](https://pachterlab.github.io/kallisto/manual) 

After the quantification task is completed for each sample, we will use the R/Bioconductor package: **tximport** to import and summarize transcript-level estimates for further transcript and gene-level analysis.

After summarization R/Bioconductor tools like edgeR and DEseq2 tools will be used for Diretential Gene Expression (DEG) analysis.


## 2. An inefficient bioinformatics ananlysis pipeline for RNA-Seq data

### Hands-on tutorial:

Under your home directory, create a folder called rnaseq

```
mkdir rnaseq
cd rnaseq
```

Create pipeline.pbs with the code below:

```
#!/bin/bash


#PBS -N rnaseq_tutorial
#PBS -l mem=12gb,nodes=1:ppn=8,walltime=8:00:00
#PBS -e error-log-${PBS_JOBID}.txt
#PBS -o output-log-${PBS_JOBID}.txt

#-------------------------------------------------------------
# set working directory, resource directory, and load modules 
#-------------------------------------------------------------

cd "$PBS_O_WORKDIR"

fastq_dir="/gpfs/data/mscbmi/rnaseq/fastq"
reference_dir="/gpfs/data/mscbmi/rnaseq/hg38"


module load gcc/6.2.0
module load udunits
module load java-jdk/1.8.0_92
module load samtools/1.6.0

module load fastqc/0.11.5
module load STAR/2.6.1d
module load subread/1.5.3
module load kallisto/0.44.0

module load R/3.5.0


#-----------------------
# Quaality Control Using FastQC (output html)
#-----------------------

mkdir RawReadQC/

fastqc \
	--extract \
	-o RawReadQC/ \
	-t 6 \
	--nogroup \
	${fastq_dir}/KO01.fastq.gz \
	${fastq_dir}/KO02.fastq.gz \
	${fastq_dir}/KO03.fastq.gz \
	${fastq_dir}/WT01.fastq.gz \
	${fastq_dir}/WT02.fastq.gz \
	${fastq_dir}/WT03.fastq.gz


#-----------------------
# Transcriptional Quantification Using Psuedo-alignment Tool kallisto
#-----------------------

mkdir -p \
	Quant/KO01

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/KO01 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/KO01.fastq.gz


mkdir -p \
	Quant/KO02

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/KO02 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/KO02.fastq.gz


mkdir -p \
	Quant/KO03

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/KO03 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/KO03.fastq.gz


mkdir -p \
	Quant/WT01

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/WT01 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/WT01.fastq.gz


mkdir -p \
	Quant/WT02

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/WT02 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/WT02.fastq.gz


mkdir -p \
	Quant/WT03

kallisto \
	quant \
	--index=${reference_dir}/transcripts.idx \
	--output-dir=Quant/WT03 \
	--single \
	--fragment-length=180 \
	--sd=20 \
	--threads=8 \
	${fastq_dir}/WT03.fastq.gz

#-----------------------
# Summarization of Expression from Transcripts to Genes
#-----------------------


#mkdir ${tutorial_dir}/Quant/

Rscript --vanilla \
	/gpfs/data/mscbmi/rnaseq/SRC/tximport.R


#-----------------------
# Sample Association Using PCA Plot
#-----------------------


Rscript --vanilla \
	/gpfs/data/mscbmi/rnaseq/SRC/pca.R

```

Submit your pipeline for execution on the HPC's batch mode:


```
$ qsub pipeline.pbs
13325292.cri16sc001
```

The pipeline will take aproximately 1 hour to complte. When completed, the pipeline should create  the following files:

```
$ ls -l
total 128
-rw------- 1 t.cri.biowksp01 t.cri.biowksp01 8142 Apr 26 12:05 error-log-13325292.cri16sc001.txt
-rw------- 1 t.cri.biowksp01 t.cri.biowksp01 3185 Apr 26 12:05 output-log-13325292.cri16sc001.txt
-rw-r----- 1 t.cri.biowksp01 t.cri.biowksp01 2782 Apr 26 11:42 pipeline.pbs
drwxr-x--- 8 t.cri.biowksp01 t.cri.biowksp01 4096 Apr 26 12:05 Quant
drwxr-x--- 8 t.cri.biowksp01 t.cri.biowksp01 4096 Apr 26 11:49 RawReadQC
```
The folder **RawReadQC** contains the Fastqc results, the folder **Quant** contains the kallisto results as well as the sumarized transcrptome profile for all samples, and the results of a PCA ananlysis.


Copy and explore the final result files to your local computer:

```
scp t.cri.biowksp01@gardner.cri.uchicago.edu:~/rnaseq/Quant/tutorial.count.kallisto.txt ./Desktop/
scp t.cri.biowksp01@gardner.cri.uchicago.edu:~/rnaseq/Quant/tutorial.pca.pdf ./Desktop/

```

## :mortar_board: Your turn: 

Develop, implement and test a better parallelization strategy to optimize the inefficient bioinformatics analysis pipeline above.

Start by executing the following commands 

```
cd
rm -rf rnaseq/
mkdir rnaseq
cd rnaseq
```




