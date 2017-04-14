
# Introduction to Parallel and distributed computing

**Center for Research Informatics, University of Chicago**

April - June 2017; Saturdays 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Learn the basics of parallel and distributed computing environments.
- Get familiar with basic models for parallel computation, and design of parallel and distributed algorithms. 

## 1. Why Parallel Computing

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


## 3. Programming strategies for Parallel Computation: 

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

## 4. Programing languages for Distributed Computing:

[BigDataScript](https://pcingola.github.io/BigDataScript/)

Develop ONE data pipeline and run exactly the same script everywhere. No matter how big the computer. Created by: [Pablo Cingolani](https://www.linkedin.com/in/pablocingolani/) Director of IT and Bioinformatics at Kew Inc.
Your program/pipeline will runs on a 25,000 cores cluster or a single CPU laptop.

![BigDataScript](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/BigDataScript.jpeg)

It is of course also possible to develop scripts for running jobs on HPC and/or multi-core processors using other scripting languages such as PERL, python, R, C++, and others.

**Running a Job on HPC using Portable Batch System (PBS)**

To run a job on CRI's HPC cluster, you will need to set up a Portable Batch System (PBS) file. This PBS file defines the commands and cluster resources used for the job.  A PBS file is a simple text file that can be easily edited with a UNIX editor, such as vi, pico, or emacs. 

**Submitting a Job**
In order to use the HPC compute nodes, you must first log in to one of the head nodes, and submit a PBS job. The **qsub** command is used to submit a job to the PBS queue and to request additional resources. The **qstat** command is used to check on the status of a job already in the PBS queue. To simplify submitting a job, you can create a PBS script and use the qsub and qstat commands to interact with the PBS queue.




## Week 2 Challange: :white_check_mark:
Read the Nature technology feature article **'Biology: The big challenges of big data'** available at: http://www.nature.com/nature/journal/v498/n7453/full/498255a.html and submit via e-mail a one page critical commentary.

