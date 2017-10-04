
# Introduction to Shell scripting; Linux commands for routine Bioinformatics tasks; and Introduction to Parallel and distributed computing

**Center for Research Informatics, University of Chicago**

Saturdays 10/07/17; 9:00 AM - 12:00 PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Learn the basics of developing useful shell scripts.
- Learn how to use important Linux command on routine tasks in bioinformatics.
- Learn the basics of parallel and distributed computing environments.
- Get familiar with basic models for parallel computation, and design of parallel and distributed algorithms. 


## 1. Introduction to Shell scripting

*Why shell scripting?* if you need to repeat a process 1000 times, you could either retype the commands 1000 times. Or you can create a small 'script' (set of instructions) that tells the system to repeat it 1000 times.

**Steps to create a Shell Script**

1. Create a file using a nano editor (or any other editor you like i.e. vi).  Name script file **with extension .sh**
2. Start the script with #! /bin/bash (or #! /bin/sh depending of your operative system and preference)
3. Write some code
4. Save the script file as filename.sh
5. Set the executable permission for filename.sh for the target users
6. To excute the script 

Let's create a very simple first shell script, **on your comand line** open a **text editor** ($ nano simple.sh).

```text
#!/bin/bash
#This script display the date, who is logged on, and your working directory
echo "Today's date is:"
date
echo "You are logged in as:"
whoami
echo "You are working on the follwoing directory:"
pwd
```
Now we need to make the script excutable 

```bash
chmod u+x simple.sh         
```
We can now execute our first shell script

```bash
./simple.sh                
```

Let's see another example **on your comand line** open a **text editor** ($ nano hello.sh). On the editor type the following:

```text
#!/bin/bash
echo "Hello World!."
echo "I am going to generate 10 files: file1.txt, file2.txt, ..., file10.txt"
for i in {1..10}
do echo "this is file $i" > file$i.txt
done
echo "Done"
```

Now we need to make the script excutable 

```bash
chmod u+x hello.sh         
```

We can now execute your shell script

```bash
./hello.sh                
ls file*                   ### List the files you just created 
```

## :mortar_board: Your turn, create a shell script to remove (delete) the 10 files you just created.

**Shell Variables**

Shell variables store data in the form of characters and numbers, the following creates a shell variable and then prints it:

```text
variable ="Hello"
echo $variable
```
Let's use the the concept of variables is a small script. Create the following script using a text editor (nano greetings.sh)

```bash
#!/bin/bash
echo "what is your name?"
read name
echo "How are you, $name?"
read remark
echo "I am $remark too!"
```
Now, execute your script.

You can also access Linux environment variables ($HOME, $USER, etc.) from within your shell scripts. Try the follwoing script:

environment.sh

```bash
#!/bin/bash
#Display user information from the system
echo "User info for user: $USER"
echo "User's HOME: $HOME"
echo "Your UID is: $UID" 
```

**Shell Loops**

**for** loops iterate through a set of values until the list is exhausted. Let's see some examples:

**for1.sh**
```bash
#!/bin/bash
for i in 1 2 3 4 5
do
  echo "Looping ... number $i"
done
```

Let's try something a bit more useful:

**for2.sh**
```bash
#!/bin/bash
for i in {1..100}
do
  echo "I will create folder # $i"
  mkdir "folder$i"
done

echo "Now I will copy the sequence file: SRR001655.fastq on each of the 100 folders"
cd ~

for i in {1..100}
do
  echo "Copying sequence on folder # $i"
  cp /group/mscbmi/lecture1/SRR001655.fastq "folder$i"
done
```

Explore the files you just created and confirm that you have a copy of SRR001655.fastq on each file

## :mortar_board: Your turn, create a shell script to remove (delete) the 100 files you just created.

**While Loops**

While loops will execute the instruction(s) *while* the evaluation condition is TRUE

**while1.sh**
```bash
#!/bin/bash
INPUT_STRING=hello
while [ "$INPUT_STRING" != "bye" ]
do
  echo "Type something, type 'bye' to quit"
  read INPUT_STRING
  echo "You typed: $INPUT_STRING"
done
```

Let's convine while loops with conditional statements:

First, copy on your home directory the file greetings.txt form /group/mscbmi/lecture1/, this file contains several ways to say hello in different languages.

```bash
cd ~
cp /group/mscbmi/lecture1/greetings.txt .
```

Let's now create a script that will recognize the language for each each line in the file:

**while2.sh**
```bash
#!/bin/bash
while read f
do
  case $f in
  hello)        echo English    ;;
  bonjour)	    echo French     ;;
  'guten tag')  echo German     ;;
  hej)          echo Swedish    ;;
  holla)        echo Spanish    ;;
  ciao)         echo Italian    ;;
  namaste)	    echo Hindi	;;
  marhaba)	    echo Arabic     ;;
  'ni hau')     echo Mandarinan ;;
        *)	    echo Unknown Language: $f
                ;;
   esac
done < greetings.txt
```

**if statements** 

The syntax for if...then...else... is:

```
if [ ... ]
then
  # if-code
else
  # else-code
fi
```
Note that **fi** is **if** backwards. This concept was already used in the previous example with **case** and **esac**


Let's try to create a script to evaluate if a number is a prime. As you know, prime numbers are numbers that are bigger than 1 and **cannot** be divided evenly by any other number except 1 and itself

**isprime.sh**
```bash
#!/bin/bash
while [ "$num" != "bye" ]
do
  echo -n "Enter a number, type 'bye' to quit "
  read num
  i=2
  while [ $i -lt $num ]
  do
   if [ `expr $num % $i` -eq 0 ]
    then
      echo "$num is not a prime number"
      echo "Since it is divisible by $i"
   fi
  i=`expr $i + 1`
  done
  echo "$num is a prime number "
done
```

**Numeric comparisons**

| Comparison    | Description                                      |
| ------------- |:------------------------------------------------:|
| n1 -ep n2     | Checks if n1 is **equal** to n2                  |
| n1 -ge n2     | Checks if n1 is **greater than or equal** to n2  |
| n1 -gt n2     | Checks if n1 is **greater than** n2              |
| n1 -le n2     | Checks if n1 is **less than or equal** to n2     |
| n1 -lt n2     | Checks if n1 is **less than** n2                 |
| n1 -ne n2     | Checks if n1 is **not equal to** n2              |

**File comparisons**

In Shell, you can test the status of filrs and directories on the Linux filesystem with the follwoing file comparisons

| Comparison    | Description                              |
| ------------- |:----------------------------------------:|
| -e file       | Checks if file exist                     |
| -d file       | Checks if file exist and is a directory  |
| -f file       | Checks if file exist and is a file       |


Let's see some examples

**checkdir.sh**
```bash
#!/bin/bash

goto_directory=/group/mscbmi

if [ -d $goto_directory ] 
then
  echo "The $goto_directory directory exists" 
  cd $goto_directory
  ls
else
  echo "The $goto_directory directory does not exist"
fi
```

**checkdirorfile.sh**
```bash
#!/bin/bash
# Check if either a directory or file exists #

location=$HOME
file_name="SRR001655.fastq"

if [ -e $location ]
  then #Folder exist
    echo "OK on the $location directory."
    echo "Now checking on the file, $file_name." 
    #
    if [ -e $location/$file_name ]
    then #File exist
      echo "OK on the filename"
      echo "Updating Current Date..." 
      date >> $location/$file_name
    #
    else #File does not exist
      echo "File does not exist"
      echo "Nothing to update" 
   fi
#
else #Directory does not exist
  echo "The $location directory does not exist."
 echo "Nothing to update" 
fi
```

## :mortar_board: Your turn,create a script or use a comand to remove (delete) the line that the previous schipt inserted on *file_name*


### :book: Learn more:

- A good reference book: Linux Command Line and Shell Scripting Bible 3rd Edition by Richard Blum and Christine Bresnahan 
- On-line tutorials:
  - Shell Scripting Tutorial: https://www.shellscript.sh 
  - Linux Shell Scripting Tutorial. A Beginner's handbook: http://www.freeos.com/guides/lsst/

____________________________________________________________


## 2. Linux commands for routine Bioinformatics tasks

**a. Processing FASTQ files**

Next Generation Sequencing platfoms generate sequence data in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format)
FASTQ format has 4 lines per read, e.g.:
```text
@SRR001665.1 071112_SLXA-EAS1_s_4:1:1:672:654/1
GCTACGGAATAAAACCAGGAACAACAGACCCAGCAC
+
IIIIIIIIIIIIIIIIIIIIIIIIEII9IIIEIIII
```
Where:

```
line 1: @SRR001665.1...- a unique identifier for the sequence read
line 2: GCTACGGAA.... - the sequence read
line 3: + - a separator between the read and the quality values (this sometimes replicates the sequence ID)
line 4: IIIIIIIIIIIIIIIIIII.... - the quality values. 
```

Exploring a fastq file using 'less', use <space> or f to go to the Next Page; b to the Previous Page; q to Exit;

```bash
cd ~
less SRR001655.fastq
```
1. How many reads are there in the file SRR001655.fastq? (Hint: use cat, grep and wc) 
2. Display the first 25 reads in file SRR001655.fastq (Hint: use head) 
3. Create the file 'short_list.txt' with a copy of the sequences in file SRR001655.fastq that contains GAGAGAGC (Hint: use grep) 
4. Write the last 10000 reads to a new file called bottom_10000.fastq (Hint: use tail)

```bash
cat SRR001655.fastq | grep '@SRR'| wc -l
```
```bash
head -n100 SRR001655.fastq
```
```bash
grep 'GAGAGAGC' SRR001655.fastq > short_list.txt
cat short_list.txt
```
```bash
tail -40000 SRR001655.fastq > bottom_10000.fastq
```
**b. Using the 'paste' command to format your data**

Showing fastq formated data as a table (i.e. in columns) can be very useful to explore the data, the 'paste' command writes lines in a file **as columns separated by a the tab character**. The command take character '-' as an option to represent the standard input, e.g.: the option '- - - -', will be trasalated as 'read four lines', and write them out as four columns:

```bash
cat bottom_10000.fastq | paste - - - - | head -10
```
```bash
cat bottom_10000.fastq | paste - - - - | head -1000 > top_1000_tab.txt 
cat top_1000_tab.txt 
```
**c. Using 'awk' to work with data in columns**

Here is an example on how to identify reads with fail sequences:
```bash
###awk '{print $3}' top_1000_tab.txt 
awk '$3 ~ /N/ {print $1,"\t",$3}' top_1000_tab.txt 
```
The follwong command will convert FASTQ to (FASTA){https://en.wikipedia.org/wiki/FASTA_format} file format:
```bash
cat SRR001655.fastq | paste - - - - | awk '{print ">"$1,$2,"\n"$3}'
```
Now try saving the result of the conversion on a new file.



**d. Text manipulation with sed**

Explore the result of the following comand:

```bash
sed 's/N/0/g' top_1000_tab.txt 
```
Check that the orininal top_1000_tab.txt  file was not altered.

To logout of gardner, type:

```bash
exit
```

________________________________________________________________
## 3. Introduction to Parallel Computing

Moore's law refers to an observation made by Intel co-founder Gordon Moore in 1965. He noticed that the number of transistors per square inch on integrated circuits had doubled every year since their invention. Moore's law predicts that this trend will continue into the foreseeable future ***the number of transistors per square inch has since doubled approximately every 18 months.***
The figure bellow shows the number of trasistors on integrated circuits chips from 1971 to 2016: 

![Moore's law](https://ourworldindata.org/wp-content/uploads/2013/05/Transistor-Count-over-time.png)


Note: Exponential growth, it is unlikely to continue indefinitely, some new studie are showing that physical limitations could be reached this year - **2017**

**Factors contributing to validate Moore's law:**

- Denser integrated circuits
- Improvements on circuits architecture 

**Measures of processor performance**

Instructions per second: 

- MIPS: A million instructions per second.
- GIPS: A billion instructions per second
- TIPS: A trillion instructions per second
- PIPS: A Peta (10^5) instructions per second

Floating-point operations per second:

- MFLOPS: A million floating-point operations per second (FLOPS)
- GFLOPS: A GigaFlop is a billion FLOPS.
- TFLOPS: A TeraFlops is a trillion FLOPS, one trillion floating-point operations per second
- PFLOPS: A Peta (10^5) floating-point operations per second

Obviously, **there is a limit to the speed of a single processor*** (the speed-of-light), however the need for **FLOPS** in science and engineering **continue to grow exponetially** (wheather forecast, astrophisics, Monte Carlo simulation of nuclear reactors, etc.)

The only alternative for **efficient** analysis of BIG DATA is paralellel and/or distributed computing.

**Biologist are now joining the bid data club**

Due to available high-throughput genomics technologies (NGS), in January this year, Illumina was the first companny to anounce that their technology will be capable of sequencing a human genome for $1000 (with their HiSeq X 10 technology, a lab would be capable of sequencing 18,000 human genomes per year. Breaking down the price of an individual genome). 

How much data is produced? takea look at the troughput of the latest Illumina sequencing platforms: 

https://www.illumina.com/systems/sequencing-platforms.html

**Two examples:**
- The European Bioinformatics Institute _(EBI)_ in the UK currently stores **20 petabytes of data** and back-ups about genes, proteins and small molecules. 
- The Genomics Data Commons initiative _(GDC)_ at the University of Chicago, **currently host 5.42 petabytes of data and 87.96 Terabytes of RAM).** See the GDC's infrastructure statistics at: (https://gdc-portal.nci.nih.gov)

## 2. Hardware Paradigms for Parallel Computation: 

**CPU Design: Multiple-Core Processor**

A **multi-core processor** is a single computing component with *two or more independent actual processing units* (called "cores"), which are units that read and execute program instructions.

Although multicore chips were designed for game playing, they are becomming more and more popular in scientific computing. Now days, dual–core, quad–core, or even sixteen–core chips are commun on modern computers. These chips attain more integrated speed with less heat and more energy efficiency than single–core.

As an example: The *Sparc64-IXfx chip* has 16 Sparc cores that run at 1.85GHz and delivers 236 GFLOPS of double-precision floating point operations. Bellow the real architectural design of the chip:

![Fujitsu's Sparc64-IXfx processor](https://regmedia.co.uk/2011/11/20/fujitsu_sparc64_ixfx_chip.jpg)

**High Performance Computing** 

High-performance computing (HPC) is the use of **parallel processing** for running advanced application programs efficiently, reliably and quickly. The term applies especially to systems that function **above a TFLOP** or 10^12 floating-point operations per second. The term HPC is not a synonym for supercomputing, technically a supercomputer is a system that performs at the currently highest operational rate of more than a PFLOPS (10^15 floating-point operations per second).

![HPC](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/HPC.jpeg)

**GRID Computing**

**Aggregation** of computer resources **remotely located in multiple locations** to reach a common goal. Grid computing is distinguished from conventional HPC systems in: 

- In Grid computing, each node can be set to perform a different task/application. 
- Grid computers also tend to be more heterogeneous and geographically dispersed.
- Grid sizes can be quite large

Although a single grid can be dedicated to a particular application, commonly a grid is used for a variety of purposes. Grids are often constructed with general-purpose grid middleware software libraries. The idea of connecting heterogeneous computers together in a Grid was developed at the European Organization for Nuclear Research [CERN](https://home.cern/about).

Grid enviroments consist at least of 3 main componets: the hardware, the middleware (like the Advanced Resource Connector ARC)  and the applications (like software).

Recomened read on this topic: [The Grid: Blueprint for a New Computing Infrastructure](https://books.google.com/books?id=ONRQAAAAMAAJ) by the University of Chicago Professor [Ian Foster](https://cs.uchicago.edu/directory/ian-foster)

**CLOUD Computing**

Cloud Computing: on-demand access to a shared pool of configurable computing resources 

![CloudComputing](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/CloudComputing.jpeg)

**Cloud Services**
- Software as a Service (SaaS) applications are hosted by a vendor and made available to customers over the Internet
- Platform as a Service (PaaS) is a way to rent hardware, operating systems, storage and network capacity (and all associated service) over the internet
- Infrastructure as a Service (IaaS) an organization outsources the equipment used to support operations, including storage, hardware, servers and networking components. The service provider owns the equipment and is responsible for housing, running and maintaining it. The client typically pays on as per-use basis.

![CloudServices](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/CloudServices.jpeg)

**Cloud Computing Advantages**

- Pay as you go: Could be less expensive compared to buying installing and maintaining your own 
- Flexibility: Theoretical infinite scalability
- Easy Access: Can be used from any computer or device with an Internet connection
- Easy updates: Updates occur across the service

**Cloud Computing Disadvantages** 

- Security Concerns
- Terms of Service
- Privacy Policies
- Charge model: I/O fees 


## Programming strategies for Parallel Computation: 

- Split the data

For an **'embarrassingly parallel problems'** (also known as perfectly parallel or pleasingly parallel), the obvious strategy is to split the big data in small manageable chunks. This startegy generally apply to 'data intensive' problems. A clasic example of this kind of problems in in bioinformatics is the use of the BLAST algorithm on huge amouns of sequencing data. See [Applications of Grid Computing in Genetics and Proteomics](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.458.8893&rep=rep1&type=pdf) 

By spliting the data and distributing the execution of the BLAST algorithm over thousand of Grid 'workers', the computational runtime of BLAST over big data can be reduced to manageable and/or acepatble.  See an application of the **GRID-BLAST** algorithm at [Using Grid technology for computationally intensive applied bioinformatics analyses](http://www.bioinfo.de/isb/2006060046/main.html)

As a result [The epitope space of the human proteome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2271163/) was defined. 

- Distribute the computational load 

For **'computationally intensive'** tasks (i.e. **searching for prime numbers**, **calculating large factorials**,  quantum chromodynamics, astrophysical and cosmological simulations, weather simulation, etc.), and for tasks were there is **strong 'data interdependence'** a Message Passing Interface (MPI) strategy is needed. MPI is a communication protocol for programming parallel computers. Both point-to-point and collective communication are supported. MPI is a message-passing application programmer interface, together with protocol and semantic specifications for how its features must behave in any implementation. MPI's goals are high performance, scalability, and portability. A good read for an introduction to programming parallel systems that use the MPI: [Parallel Programming with MPI](http://www.cs.usfca.edu/~peter/ppmpi/)

Recomended reading: [rapidGSEA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1244-x), a bioinformatics application for speeding-up the gene set enrichment analysis on multi-core CPUs and CUDA-enabled GPUs

- Access more memmory

There are tasks that are both computational and data intensive, they usually need to have simultaneous access to the full data model
(i.e. multiplying two big matrices), to compute those kind of problems, it is needed to access HPC resources with nodes that share large amouns of RAM memmory. In Bioinformatics a commun problem that often needs such setting is **De novo transcriptome assembly** from RNA-Seq data this task ususallu demand ~1G of RAM per ~1M pairs Illumina reads (i.e. tools like: [SOAPdenovo-trans](http://soap.genomics.org.cn/SOAPdenovo-Trans.html), [TransABYSS](https://github.com/bcgsc/transabyss) and [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki))

See Robert Bukowski's Trinity workshop and hands-on exercises at: http://cbsu.tc.cornell.edu/lab/doc/Trinity_workshop_Part1.pdf 

## Programing languages for Distributed Computing:

[BigDataScript](https://pcingola.github.io/BigDataScript/)

Develop ONE data pipeline and run exactly the same script everywhere. No matter how big the computer. Created by: [Pablo Cingolani](https://www.linkedin.com/in/pablocingolani/) Director of IT and Bioinformatics at Kew Inc.
Your program/pipeline will runs on a 25,000 cores cluster or a single CPU laptop.

![BigDataScript](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/BigDataScript.jpeg)

It is of course also possible to develop scripts for running jobs on HPC and/or multi-core processors using other scripting languages such as PERL, python, R, C++, and others.

**Running a Job on HPC using Portable Batch System (PBS)**

To run a job on CRI's HPC cluster, you will need to set up a Portable Batch System (PBS) file. This PBS file defines the commands and cluster resources used for the job.  A PBS file is a simple text file that can be easily edited with a UNIX editor, such as vi, pico, or emacs. 

**Submitting a Job**
In order to use the HPC compute nodes, you must first log in to one of the head nodes, and submit a PBS job. The **qsub** command is used to submit a job to the PBS queue and to request additional resources. The **qstat** command is used to check on the status of a job already in the PBS queue. To simplify submitting a job, you can create a PBS script and use the qsub and qstat commands to interact with the PBS queue.

## Week 2 Homework: :house: (Graded!)

1. Download the Generic Feature Format (GFF) version of the Saccharomyces Cerevisiae yeast genome to your home directory. 

    The genome is available at: http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff 

    A description of the GFF format is available [here](http://www.sequenceontology.org/gff3.shtml)

2. Using the Linux command line, asnwer the following questions:
```
    a. How many genes are there in the data?
    
    b. How many genes on chromosome 2?
    
    c. How many mRNAs are there on chromosome 10?
    
    d. Describe what the following command do: sed '/#/d' saccharomyces_cerevisiae.gff > features.gff
    
    e. How many features (lines) are there in the file features.gff
    
    f. Describe what the following command do: cut -f 1 features.gff | sort | uniq -c | sort -k1n
    
    g. Which chromosome is the longest and which one is the shortest?
    
    h. Describe what the following command do: sed 's/chrI/chr1/g' features.gff > new_features.gff
```
:point_right: Send your homework via e-mail, your homework should have the answers and commands used for each item. 
    



