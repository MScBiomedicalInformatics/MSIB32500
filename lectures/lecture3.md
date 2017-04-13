
# Using a High Performance Computing cluster for bioinformatics analysis 

**Center for Research Informatics, University of Chicago**

April - June 2017; Saturdays 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Learn how develop HPC batch scripts using the Linux shell to execute code and develop analysis pipelines to perform computationally expensive bioinformatics tasks

## 1. Introduction

In this hand-on tutorial, you will learn how to use CRI's High Performance Computing (HPC) cluster TARBELL to perfom basic Bioinformatics analysis. TARBELL is a cluster with 2640 cores running **Scientific Linux 6** Update 4 with 2.6.32-358.18.1.el6.x86_64 kernel.  CRI also provides 1.3-petabyte storage space to be provisioned as **lab-shares** for BSD researchers. 

In this tutorial you will learn how to **access TARBELL cluster**, **transfer files**, create simple **PBS** job submission scripts and submit jobs in four hand-on exercises. Upon finishing the hand-on practices you will be able to use TARBELL cluster to perform simple quality control and sequence alignment analysis on Next Generation Sequencing (NGS) data.

**Prerequisites**

We will use SSH (Secure Shell) to connect to CRI's HPC. SSH is included in all standard Linux and Mac Operating Systems. If you are using MacOS or LinuxOS you can open the **Terminal** application.  If you are using Microsoft Windows, you'll have to install a SSH client on your computer. 

:pushpin:**Microsoft Windows user** you will need to install a tool for remote computing: [MobaXterm](http://mobaxterm.mobatek.net) and/or
[PuTTY](http://www.putty.org). PuTTY user go to: http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html locate the apropiate binary file (executable putty.exe file) for your hardware (32bits or 64 bits laptop). Place this file on your destop (or other folder that you can access easyly). You will need to 'double-click' this file each time you need to access it.

Open a PuTTY session, on 'Host Name (or IP address)' type: **tarbell.cri.uchicago.edu**,  select SSH as the Connection Type, verify the port number in the'Port' Box is **22**. Press the 'Open' button, then type in the provided **username and password** when prompted. Type **yes** if you are prompted to accept a key before entering username.


## 2. Login

For MacOS or Unix/Linux users:
1. Open a terminal session.
2. Connect to the login node of TARBELL cluster:

```bash
$ssh username@tarbell.cri.uchicago.edu
```
**CRI's TARBELL HPC diagram**

![tarbell](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/tarbell.jpeg)

By default, you will be logged into to the 'log-in' nodes **(in01 or in02)**, from there you will have access to the **storage/labshares** and **scratch** spaces (1Gb and 56 Gb speeds respectively). The log-in node also have access to the **scheduler and RM** nodes, the scheduler is the interface, to the **compute nodes (cn1, cn2, .. cn500)**. The log-in nodes also have access to the **cri-syncmon** node which is a dedicated Input gate node with several protocols/ports for data transfer enabeled.


Create a working directory and four sub-directories under your **home** directory for this had-on tutorial:

```bash
$ mkdir ~/mscbmi
$ cd mscbmi
$ mkdir Ex1 Ex2 Ex3 Ex4
```
## 3. Transfering data files

You can transfer files from your local computer to your home directory on the cluster or download the files from public databases and repositories. 

To download data from a website directly to your **working directory** on TARBELL cluster, you can use either the command **wget** or **curl** 

**Exercise 1.1:** Download a microarray expression raw data file from NCBI's Gene Expression Omnibus (GEO) (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31736) to your home directory on TARBELL cluster.

```bash
$ cd Ex1
$ wget 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE31736&format=file' -O GSE31736_RAW.tar
$ ls  
$ rm GSE31736_RAW.tar
```

or
```bash
$ curl 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE31736&format=file' -o GSE31736_RAW.tar
$ ls
$ rm GSE31736_RAW.tar
```
The '-o' option redirect the content to a file, otherwise it will appear on stdout, i.e. the screen.

To **upload the files from your local computer to TARBELL cluster**, you can use **scp or rsync** commands on Mac/Unix/Linux or Winscp on Windows.

**Exercise 1.2: ** Upload a file from your laptop/computer to your home directory on TARBELL cluster:

First, download the file available at the following link, to your local computer:

https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/data/GSE31736_RAW.tar 

Using computer command line (open a new terminal), navigate (use the cd command) to the directory where file 'GSE31736_RAW.tar' is located.

For Mac/Unix/Linux users:

To upload files from your local computer using **scp** use the following command:

```bash
scp PATH_TO_GSE31736_RAW.tar username@tarbell.cri.uchicago.edu:~/mscbmi/Ex1
```
Your command should look like: 

```bash
scp ./GSE31736_RAW.tar jandrade@tarbell.cri.uchicago.edu:~/mscbmi/Ex1
```
To upload files from your local computer using **rsync** use the following command:

```bash
rsync -avz GSE31736_RAW.tar username@tarbell.cri.uchicago.edu:~/mscbmi/Ex1
```
Windows users can also use GUI tools like WinSCP (http://winscp.net/download/winscp514setup.exe) to transfer files.


## 4. Running jobs on TARBELL HPC cluster

**Note: DO NOT RUN JOBS on the login nodes of the cluster. Always submit jobs to the compute nodes (qsub), or use the interactive mode (qsub -I)**

In this section, you will learn how to execute jobs on CRI's TARBELL cluster. We will use a tool for raw data quality control of NGS data as an example. The tool we will be using is a Java based program called [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
This tool provides a modular set of QC analyses that can help you to evaluate the quality of your sequences, this is in general the first step on any NGS analysis pipeline.  

TARBELL cluster uses **Torque** as a resource manager (Provides low-level functionality to start, hold, cancel and monitor jobs) and **Moab** as Work-load Manager (job scheduler) to manage the cluster resources. Torque/Moab is based on the **Portable Batch System (PBS)** originally developed by NASA in the early 1990s. As such, **Torque/Moab uses PBS directives (commands)** to receive job requests from users.

 **What does a PBS script look like?**
 
![PBS](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/PBS.jpeg)

Torque provides user commands such as **qsub, qdel, qstat, etc.,** which are used to submit, delete and check the status of jobs on the cluster. 

TARBELL cluster supports two types of job submission: **Interactive**  and **batch** modes 

- Interactive mode (qsub -I): You will execute your code/commands in an 'interactive' command line window, you will be able to see the result/output of each action interactively, this mode is useful for testing and 'debugging' code.

- Batch mode: In a batch mode you first write a PBS script with all the instructions/code you want to execute, and then you submit that script to the scheduler. The batch job script contains all the information needed, such as the location of the input and output files, as well as run parameters. Once the batch job starts, you can log off and the job will remain running. 

** Exercise 2. Running job in interactive mode**

In this exercise, you will conduct the quality control analysis on a "good" quality sequence file and on a "bad" quality sequence file using FastQC program. 

First you will need to download two raw sequence files to your **home** directory on TARBELL cluster. 

You can do so using the command **wget or curl**:

- seqGood.fastq (https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/data/seqGood.fastq) 
- seqBad.fastq (https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/data/seqBad.fastq) 

Your command should look like: 

```bash
cd ~/mscbmi/Ex2
wget ftp://logia.cri.uchicago.edu/bioinformatics/MSIB32500/Lecture3/Ex2/seqGood.fastq
wget ftp://logia.cri.uchicago.edu/bioinformatics/MSIB32500/Lecture3/Ex2/seqBad.fastq
ls -l
```

Let us first have a look at the sequence read fastq file. The following commands show the first and last two reads, as well as the total number of reads in seqGood.fastq file.

```bash
head -8 seqGood.fastq
tail -8 seqGood.fastq
grep --count ^@ERR030881 seqGood.fastq
```
Since each read in a fastq file has 4 lines, to find out how many reads are there in the fastq file, you can also use the folloing command:

```bash
wc -l seqGood.fastq
```
And divide the result by 4.


To **launch an interactive session using one core, use the following command:

```bash
qsub -I

cd ~/mscbmi/Ex2
```
Now You can run FastQC by typing:

```bash
module load fastqc
fastqc seqGood.fastq
ls
```
This will generate a self-contained directory called **"seqGood_fastqc.html"** which contains an HTML formatted report that can be loaded into a browser and a compressed file with all the results named **"seqGood_fastqc.zip**

Explore the files created by the execution of FastQC.

Other usefull **module** commands:

- module avail: List all 'modules' (tools) available to the current user. 
  Try: 
  $ module av R
  $ module load R/3.2.0
  $ which R
  $ module load R/3.3.1  
   $ which R
- module switch: Switch to a different version of the module 
  Try:
  $ which R
  $ module switch R/3.2.0 R/3.1.0
  $ which R

- module unload: Unloada module
  Try:
  $ which R
  $ module unload R
  which R

Let's now use the **scp** command, to copy the result files from your working directory in TARBELL to your local computer, so that we can vizualize the results. 

Open a **new command line on your local** computer and then run the following commands (you will need to use your own username and password):

```bash
scp t.cri.biowksp40@tarbell.cri.uchicago.edu:~/mscbmi/Ex2/seqGood_fastqc.html ./
scp t.cri.biowksp40@tarbell.cri.uchicago.edu:~/mscbmi/Ex2/seqGood_fastqc.zip ./
```

Now, in your local computer, go to the folder were you just copied the files (in my case is ./ or my home) and open the **seqGood_fastqc.html** file to explore the FastQC results. You can also unzip the **seqGood_fastqc.zip** file and review the **summary.txt** file for a summary of **PASS** and **FAIL** tests.


Repeat the steps to run FastQC on the **seqBad.fastq** file and compare the QC reports with the previous results (note you are still in the **interactive** mode session).

```bash
cd ~/mscbmi/Ex2
fastqc seqBad.fastq
exit
```


**Note** that after we run FastQC on the **seqBad.fastq** file, we used the command: **exit** to exit the **interactive** mode session.

Open a new command line on your local computer and then run the following commands (you will need to use your own username and password):

```bash
scp t.cri.biowksp40@tarbell.cri.uchicago.edu:~/mscbmi/Ex2/seqBad_fastqc.zip ./
scp t.cri.biowksp40@tarbell.cri.uchicago.edu:~/mscbmi/Ex2/seqBad_fastqc.html ./
```
Let's compare the FastQC results from  **seqGood.fastq** with **seqBad.fastg**, we will take a look at the Per base quality test:

**seqGood.fastq**

![good](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/perBaseGood.png)

**seqBad.fastg**

![bad](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/perBaseBad.png)

**Exercise 3: Running jobs in batch mode**

In practice, you would likely want to evaluate more than one or two sequence files at the time. Instead of running FastQC sequentially on each file, you can take advantage of the power of batch job submission. In order to do so, you need to create a job submission script (a **PBS** script). 

As described before, a job submittion script will look like the following:

```bash
#!/bin/bash

###############################
# Resource Manager Directives #
###############################
### Set the job's name
#PBS -N jobname
### Select the shell you would like the script to execute
#PBS -S /bin/bash
### Set the expected runtime as walltime=HH:MM:SS
#PBS -l walltime=0:59:00
### Set the number of CPU cores for your job. This example will allocate two cores on a single node
#PBS -l nodes=1:ppn=2
### Inform the scheduler of the amount of memory you expect to use. Use units of ‘b’, ‘kb’, ‘mb’, or ‘gb’
#PBS -l mem=512mb
### Set the destination for your program’s standard output (stdout) and error (stderr)
#PBS -o $HOME/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -e $HOME/${PBS_JOBNAME}.o${PBS_JOBID}

#################
# Job Execution #
#################
### define the program/s and/or command/s to be executed

./command &> output

```

In the following exercise, you are going to perform the quality control of two RNA-seq dataset from Illumina’s Human [BodyMap 2.0 project](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/). The sequence data, generated on HiSeq 2000 instruments in 2010, consist of 16 different human tissue types. We will use a subset of the data that contains 50bp paired-end reads (PE) from 2 tissues.

First you need to download the compressed sequence read files (*.fastq.gz) to your working directory on TARBELL cluster using the command **wget or curl**. Note that FastQC program accepts both .fastq and fastq.gz file formats.

```bash
cd ~/mscbmi/Ex3
wget ftp://logia.cri.uchicago.edu/bioinformatics/MSIB32500/Lecture3/Ex3/*.gz
ls
```
You can check the compressed fastq.gz file using the command zcat or zless:

```bash
zless heart_ERR030886.sample.1.fastq.gz
```
Type **q** to exit 

Next we will create a script for the QC for each tissue (four files) based on the above template. You can use any text editor write your script. Here we will use **nano** tool that is installed on TARBELL.

```bash
nano run_fastqc_heart.pbs
```

Copy & paste the following script to the **nano** text editor:

```
#!/bin/bash
###############################
# Resource Manager Directives #
###############################
### Set the name of the job
#PBS -N run_fastqc_heart
### Select the shell you would like to script to execute within
#PBS -S /bin/bash
### Inform the scheduler of the expected runtime
#PBS -l walltime=0:59:00
### Inform the scheduler of the number of CPU cores for your job
#PBS -l nodes=1:ppn=1
### Inform the scheduler of the amount of memory you expect
#PBS -l mem=512mb
### Set the destination for your program’s output and error
#PBS -o $HOME/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -e $HOME/${PBS_JOBNAME}.o${PBS_JOBID}

#################
# Job Execution #
#################
# load the fastqc tool
module load fastqc
# set the paths
seqPath=~/mscbmi/Ex3
seqfile1=$seqPath/heart_ERR030886.sample.1.fastq.gz
seqfile2=$seqPath/heart_ERR030886.sample.2.fastq.gz
# run fastqc
fastqc -o $seqPath $seqfile1 &> $seqPath\/heart.fastqc.log
fastqc -o $seqPath $seqfile2 &>> $seqPath\/heart.fastqc.log
```
Note: To save a file in nano, you can use Ctrl-O. To close nano Ctrl-X.

Next we will submit this job in batch mode:

```bash
qsub run_fastqc_heart.pbs
```
To check the status of your job use:

```bash
qstat
```
To watch the status of your job and keep a window to do so, use:

```bash
watch qstat
```
Use **Ctrl-c** to exit the watch window.

To delete a batch job, simply type **qdel**, followed by the Job id and return.

Now, repeat the above procedure for the kidney pair-end sequence files. You will neet to create a script named **run_fastqc_kidney.pbs** using nano.

```bash
nano run_fastqc_kidney.pbs
```
Copy & paste the following script to the **nano** text editor:

```
#!/bin/bash
###############################
# Resource Manager Directives #
###############################
### Set the name of the job
#PBS -N run_fastqc_heart
### Select the shell you would like to script to execute within
#PBS -S /bin/bash
### Inform the scheduler of the expected runtime
#PBS -l walltime=0:59:00
### Inform the scheduler of the number of CPU cores for your job
#PBS -l nodes=1:ppn=1
### Inform the scheduler of the amount of memory you expect
#PBS -l mem=512mb
### Set the destination for your program’s output and error
#PBS -o $HOME/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -e $HOME/${PBS_JOBNAME}.o${PBS_JOBID}

#################
# Job Execution #
#################
# load the fastqc tool
module load fastqc
# set the paths
seqPath=~/mscbmi/Ex3
seqfile1=$seqPath/heart_ERR030886.sample.1.fastq.gz
seqfile2=$seqPath/heart_ERR030886.sample.2.fastq.gz
# run fastqc
fastqc -o $seqPath $seqfile1 &> $seqPath\/kidney.fastqc.log
fastqc -o $seqPath $seqfile2 &>> $seqPath\/kidney.fastqc.log
```

Save and close nano. Then submit the job and check the status:

```bash
qsub run_fastqc_kidney.pbs
qstat
```
You can review the QC results in the corresponding directories as you did in Exercise 2.

Note: If you need to submit many jobs to the computation nodes, you don't have to type the above commands multiple times. Instead, you can simply create a new shell script named say **submit_jobs.sh** containing a script like the following:

```bash
#!/bin/bash
qsub run_fastqc_heart.pbs
sleep 2
qsub run_fastqc_kidney.pbs
```
Then execute the script:

```bash
sh submit_jobs.sh
```
You may need to give the user (u) the authority to execute (x) the script you just created (submit_jobs.sh) before you can run/execute your shell script:

```bash
chmod u+x submit_jobs.sh
$ ./submit_jobs.sh
```

Next up we are going to align the sequence reads to the reference human genome (hg19) using two alignment tools - BWA (http://bio-bwa.sourceforge.net) and Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Both tools use the Burrows--Wheeler transformation to reduce the memory requirement for the sequence alignment. Each has its own set of limitations, for example, the lengths of reads it accepts, how it outputs read alignments, how many mismatches there can be, whether it produces gapped alignments, etc.). Please note for RNA-Seq data, people use splicing-aware alignment tools such as TopHat and STAR aligners.

** Exercise 4: Sequence alignment on HPC cluster**

First, copy the pair-end Illumina sequence files for heart from Exercise 3

```bash
cp ~/mscbmi/Ex3/*.fastq.gz ~/mscbmi/Ex4
```

To run the alignment using **bwa** on the cluster, you need to create a new job submission script named **run_bwa_heart.pbs** in nano:

```
#!/bin/bash
###############################
# Resource Manager Directives #
###############################
### Set the name of the job
#PBS -N run_bwa_heart
### Select the shell you would like to script to execute within
#PBS -S /bin/bash
### Inform the scheduler of the expected runtime
#PBS -l walltime=0:59:00
### Inform the scheduler of the number of CPU cores for your job
#PBS -l nodes=1:ppn=4
### Inform the scheduler of the amount of memory you expect
#PBS -l mem=2gb
### Set the destination for your program’s output and error
#PBS -o $HOME/mscbmi/Ex4/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -e $HOME/mscbmi/Ex4/${PBS_JOBNAME}.o${PBS_JOBID}

#################
# Job Execution #
#################
# load the programs to be executed
module load bwa
module load samtools
# set the input files and program paths
seqPath=~/mscbmi/Ex4
outPath=~/mscbmi/Ex4/bwa
seqfile1=$seqPath/heart_ERR030886.sample.1.fastq.gz
seqfile2=$seqPath/heart_ERR030886.sample.2.fastq.gz
readgroup=$outPath/heart_ERR030886.sample
referenceSeq=/group/referenceFiles/Homo_sapiens/UCSC/hg19/hg19.GATKbundle.1.5/ucsc.hg19.fasta

# create output directory if  it does not exist
if [ ! -d $outPath ]; then mkdir $outPath; fi

# align the pair-ended sequences, use 4 threads
bwa aln -t 4 $referenceSeq $seqfile1 >  $readgroup.1.sai
bwa aln -t 4 $referenceSeq $seqfile2 >  $readgroup.2.sai
bwa sampe $referenceSeq $readgroup.1.sai $readgroup.2.sai $seqfile1 $seqfile2 > $readgroup.sam
# convert sam file to bam file; sort and index bam file
samtools view -F 4 -Sb $readgroup.sam > $readgroup.bam
samtools sort $readgroup.bam -o $readgroup.sorted.bam
samtools index $readgroup.sorted.bam
```
Next submit the job to compute nodes and monitor the job status:

```bash
qsub run_bwa_heart.pbs
qstat
```
Note:
- The program takes just a few minutes to complete because the input sequence files used in this exercise are a stripped-down version of the original compressed fastQ files, The tipical size of a complete RNA-Seq fastq file is about 6.1GB. Running the BWA aligner on a real file would take much longer time, thus the utility of HPC enviroments to reduce the time to manageable.
- BWA and Bowtie2 both support multi-threading. To enable it, turn on the option -t in BWA and -p in Bowtie2. In above script, 4 threads are used for the mapper, you also need to inform the scheduler of the number of CPU cores for your job by setting ppn=4. While multi-threading can speed up the mapping, specifying too many threads could cause your job to be waiting in the queue for sufficient resources to become available. Set the number of threads to 4, 8, or 16 depending on the read file sizes and current job load on the cluster.

We are now going to run the alignment using **bowtie2** on the cluster, you need to create a new job submission script named **run_bowtie2_heart.pbs** using nano editior:


```bash
#!/bin/bash
###############################
# Resource Manager Directives #
###############################
### Set the name of the job
#PBS -N run_bowtie2_heart
### Select the shell you would like to script to execute within
#PBS -S /bin/bash
### Inform the scheduler of the expected runtime
#PBS -l walltime=0:59:00
### Inform the scheduler of the number of CPU cores for your job
#PBS -l nodes=1:ppn=4
### Inform the scheduler of the amount of memory you expect
#PBS -l mem=2gb
### Set the destination for your program’s output and error
#PBS -o $HOME/mscbmi/Ex4/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -e $HOME/mscbmi/Ex4/${PBS_JOBNAME}.o${PBS_JOBID}

#################
# Job Execution #
#################

module load samtools bowtie2

# set the file and program paths
seqPath=~/mscbmi/Ex4
outPath=~/mscbmi/Ex4/bowtie2
seqfile1=$seqPath/heart_ERR030886.sample.1.fastq.gz
seqfile2=$seqPath/heart_ERR030886.sample.2.fastq.gz
readgroup=$outPath/heart_ERR030886.sample
referenceSeq=/group/referenceFiles/Homo_sapiens/UCSC/hg19/Sequence/IlluminaBowtie2Index/genome

# create output directory if it does not exist
if [ ! -d $outPath ]; then mkdir $outPath; fi

# align the pair-ended sequences, use 4 threads
bowtie2 -p 4 -x $referenceSeq -1 $seqfile1 -2 $seqfile2 -S $readgroup.sam

# convert sam file to bam file, sort and index bam file
samtools view -F 4 -Sb $readgroup.sam > $readgroup.bam
samtools sort $readgroup.bam -o $readgroup.sorted.bam
samtools index $readgroup.sorted.bam
```
Now, use the following commands to submit and monitor the job:

```bash
qsub run_bowtie2_heart.pbs
qstat
```

Once the alignment jobs are finished, take a look at the files generated:

```bash
ls -l ~/mscbmi/Ex4/bwa
```
You should see the following 6 files:

- heart_ERR030886.sample.1.sai
- heart_ERR030886.sample.2.sai
- heart_ERR030886.sample.bam
- heart_ERR030886.sample.sam
- heart_ERR030886.sample.sorted.bam
- heart_ERR030886.sample.sorted.bam.bai

and

```bash
ls -l ~/mscbmi/Ex4/bowtie2
```
Should list the follwoing files:

- heart_ERR030886.sample.bam
- heart_ERR030886.sample.sam
- heart_ERR030886.sample.sorted.bam
- heart_ERR030886.sample.sorted.bam.bai


You can now visualize the alignments on genome visualization tools such as the Integrative Genomics Viewer [IGV](http://www.broadinstitute.org/igv). Download and install the IGV to your local computer. 
Trasfer the *bam* and *bai* files from TARBELL to your local computer. 
To transfer a file from TARBELL to your local computer, open a new comand line in your local computer and use the **scp** command: 

```bash
scp jandrade@tarbell.cri.uchicago.edu:~/mscbmi/Ex4/bwa/heart_ERR030886.sample.sorted.bam ./
scp jandrade@tarbell.cri.uchicago.edu:~/mscbmi/Ex4/bwa/heart_ERR030886.sample.sorted.bam.bai ./
```
Move the files you just downloaded to a folder in your Desktop (bwa_r)

```bash
scp jandrade@tarbell.cri.uchicago.edu:~/mscbmi/Ex4/bowtie2/heart_ERR030886.sample.sorted.bam ./
scp jandrade@tarbell.cri.uchicago.edu:~/mscbmi/Ex4/bowtie2/heart_ERR030886.sample.sorted.bam.bai ./
```
Move the files you just downloaded to a another folder in your Desktop (bowtie2_r)

In your local computer, open the the Integrative Genomics Viewer and load the **bam** files *(File -> load from file)* from each result folder.
Now search for the gene ** TOMM40L** on IGV. You should be able to see the mapping results for the two tools, as in the figure bellow:

![igv](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/IGV.png)

## Week 3 Homework: :house: 

Investigate how to use the UCSC genome browser online, to replace the IGV visualiztion of the previous results. Send me an e-mail with the image of your UCSC genome browser session showing the **TOMM40L** gene.

