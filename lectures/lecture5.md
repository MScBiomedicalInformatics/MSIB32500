
# Data Visualization in  R  

**Center for Research Informatics, University of Chicago**

Saturdays 04/28/2018; 9:00 AM - 12:00 PM

**Instructor:** Jorge Andrade, Ph.D.

## Learning Objectives

- Develop a working knowledge of how to use R and Bioconductor for basic data visualizations of genomics data
- Develop a working knowledge of how to use R and Bioconductor for basic multiple testing for gene expression data

## 1. Introduction

R Programming offers a set of builtin function and libraries (i.e. base, ggplot2, leaflet, lattice) to build data visualizations and present data analysis results. In this tutorial, you will learn the steps to create common as well as advanced visualizations in R using the base and some ggplot2 functions.

## 2. Data download

For this tutorial, we will use the  Gene expression dataset from **Golub et al.** and the R Bioconductor package for **Resampling-based multiple hypothesis testing: multtest**. 

* The **multtest** package includes non-parametric bootstrap and permutation resampling-based multiple testing procedures (including empirical Bayes methods) for controlling the family-wise error rate (FWER), generalized family-wise error rate (gFWER), the tail probability of the proportion of false positives (TPPFP), and false discovery rate (FDR). Several choices of bootstrap-based null distribution are implemented (centered, centered and scaled, quantile-transformed). Single-step and step-wise methods are available. Tests based on a variety of t- and F-statistics (including t-statistics based on regression parameters from linear and survival models as well as those based on correlation parameters) are included. 
* When probing hypotheses with t-statistics, users may also select a potentially faster null-distribution which is multivariate normal with mean zero and variance-covariance matrix derived from the vector influence function. 
* Results are reported in terms of adjusted p-values, confidence regions, and test statistic cutoffs. 
* The procedures are directly applicable to identifying differentially expressed genes in DNA microarray experiments.


#### Installing multtest Bioconductor package:

```{r}
> source("http://bioconductor.org/biocLite.R")
> biocLite("multtest")
> library(multtest)
> data(golub)
```

First, let us take a look what is available on the **golub** dataset:

```{r}
?golub
```    
As we can see, it contains gene expression data from the leukemia micro-array study of Golub et al. There are 3051 genes and 38 tumor mRNA samples. 3 datasets are available:

- **golub:** matrix of gene expression levels for the 38 tumor mRNA samples, rows correspond to genes (3051 genes) and columns to mRNA samples.
- **golub.cl:** numeric vector classifying the tumors, 27 acute lymphoblastic leukemia (ALL) cases (code 0) and 11 acute myeloid leukemia (AML) cases (code 1).
- **golub.gnames:** a matrix containing the names of the 3051 genes for the expression matrix golub. The three columns correspond to the gene *index,* *ID,* and *Name*, respectively.
        

## 3. Basic plots  

Let's start by looking into the golub data (3051 gene expression values for 38 tummor mRNA samples):

```{r}
> dim(golub)
[1] 3051   38
> golub[1:4, 1:7]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
[1,] -1.45769 -1.39420 -1.42779 -1.40715 -1.42668 -1.21719 -1.37386
[2,] -0.75161 -1.26278 -0.09052 -0.99596 -1.24245 -0.69242 -1.37386
[3,]  0.45695 -0.09654  0.90325 -0.07194  0.03232  0.09713 -0.11978
[4,]  3.13533  0.21415  2.08754  2.23467  0.93811  2.24089  3.36576
```

```{r}
> golub.cl
[1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1
```
The gene names for golub are stored in **golub.gnames**, let's list the first 10 rows:

```{r}
> golub.gnames[1:10, ]
[,1] [,2]                                               [,3]                         
 [1,] "36" "AFFX-HUMISGF3A/M97935_MA_at (endogenous control)" "AFFX-HUMISGF3A/M97935_MA_at"
 [2,] "37" "AFFX-HUMISGF3A/M97935_MB_at (endogenous control)" "AFFX-HUMISGF3A/M97935_MB_at"
 [3,] "38" "AFFX-HUMISGF3A/M97935_3_at (endogenous control)"  "AFFX-HUMISGF3A/M97935_3_at" 
 [4,] "39" "AFFX-HUMRGE/M10098_5_at (endogenous control)"     "AFFX-HUMRGE/M10098_5_at"    
 [5,] "40" "AFFX-HUMRGE/M10098_M_at (endogenous control)"     "AFFX-HUMRGE/M10098_M_at"    
 [6,] "41" "AFFX-HUMRGE/M10098_3_at (endogenous control)"     "AFFX-HUMRGE/M10098_3_at"    
 [7,] "42" "AFFX-HUMGAPDH/M33197_5_at (endogenous control)"   "AFFX-HUMGAPDH/M33197_5_at"  
 [8,] "43" "AFFX-HUMGAPDH/M33197_M_at (endogenous control)"   "AFFX-HUMGAPDH/M33197_M_at"  
 [9,] "45" "AFFX-HSAC07/X00351_5_at (endogenous control)"     "AFFX-HSAC07/X00351_5_at"    
[10,] "46" "AFFX-HSAC07/X00351_M_at (endogenous control)"     "AFFX-HSAC07/X00351_M_at"    
```
The column 2 contains the gene name, column 3 the probeID.

### Task: Investigate the expression values for a given gene in all 38 mRNA samples 

* First we will search for the position (row) of the gene of interest (say: **CCND3**)

```{r}
> grep("CCND3", golub.gnames[,2])
[1] 1042

```
The mRNA expression for the 38 samples for gene **CCND3** is in line 1042

### Scatter plot 

* Make a scatter plot with basic graph **plot()** function: 

```{}
> mygene <- golub[1042,]
> plot(mygene)
```
![scatterplot](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/scatterplot.png)

* Two import arguments in the plot function
    + `pch=`: type of symbols (point types)
    + `col=`: color of symbols  

Let's try them:
    
```{r}
> par(mar=c(6,4,2,1)) # Set margins
> plot(mygene, pch=15, col="red")
> mtext("Gene Expression values of Gene CCND3")
``` 

![RED](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/scatterplotRED.png)

Try altering the values of 'pch' and 'col'
   
```{r}
> plot(mygene, pch=22, col="green")
```

R colors palete can be found [here](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/Rcolor.pdf)

### Bar plot 

To compare and contrast the gene expression of several samples, we can us the **barplot()** command in R. For example, to plot the expression of 4 Cyclin genes expression values in 3 ALL and 3 AML samples:

```{r}
> grep("CCN", golub.gnames[,2])
[1]   85 1042 1212 2240
> golub.gnames[c(85,1042,1212,2240),2]
> mygenelist <- golub[c(85,1042,1212,2240),c(1:3, 36:38)]
> barplot(mygenelist, beside=FALSE)
```

![barplot](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/barplotplain.png)

We are now able to see the 'stacked' bar plot of expression of the 4 genes in 6 samples. It is still difficult to read, we will improve the plot: 


* There are several import arguments for us to use:
    + `beside=`: to be portrayed as 'stacked' bars or 'juxtaposed' bars
    + `horiz=`: bars to be drawn vertically or horizontally
    + `ylim=`: limits for the y axis  
    + `xlim=`: limits for the x axis
    + `col=`: color choices

Now let's try:

```{r}
> par(mar=c(6,6,1,1)) #Define margins 
> barplot(mygenelist, beside=T, horiz=F, col=c(4:7), ylim=c(-2,2))
> mtext("Expression of 4 'Cyclin' genes for 3 ALL and 3 AML samples", side=1, at=12.5, line=1)
> legend(16,2.1, golub.gnames[c(85,1042,1212,2240),2] ,cex=1.2, col=c(4:7), pch=15, bty="n")
```
![barplotB](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/newbox4c.png)


#### Task: Fing a way to add a red vertical line separating the expression of the ALL from the AML samples

### Histograms and Density Plots  


* You can create histograms with the function **hist(x)** where **x** is a numeric vector of values to be plotted. The option **freq=FALSE** plots probability densities instead of frequencies. The option **breaks=** controls the number of bins.

    + `freq=T`: represents frequencies (default)
    + `freq=F`: represents probability densities.
 
```{r}
> par(mar=c(6,4,2,2))
> hist(mygene, freq=F, main = "Histogram of gene CCND3 expression values", ylab = "Density", xlab = "CCND3")
```
![histo](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/histogram.png)

* Other import arguments
    + `col=`: color to fill the bars
    + `border=`: color of the border around the bars
    + `breaks=no.`: the number of cells 
   
```{r}
> hist(mygene, col="green", border=2, breaks=50, freq=F, 
main = "Histogram of gene CCND3 expression values", ylab = "Density", xlab = "CCND3")
```

![histoG](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/histogramGreen.png)



### Line charts 

* Line chart of the expression values of gene CCND3 in 38 mRNA samples 

```{r}
> par(mar=c(6,4,2,1))
> plot(mygene, type="l")
> mtext("CCND3 gene expression values", side=1, line=5)
```
![linechart](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/linechart.png)

* Import arguments for the plot function
    + `type="l"`: plot line chart
    + `lty=`: type of lines
    + `col=`: color of lines
    + `lwd=`: width of lines

Let's use them:

```{r}
> par(mar=c(6,4,2,1)) #margin in line
> plot(mygene,type="l", col="red", lwd=3, lty=2);
> mtext("Gene expression values of CCND3", side=1, line=5)
``` 

![plotlinered](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/plotlinered.png)

One can also add lines to an existing plot with the command **lines()** 

For example, to add a **density curve line** to a histogram we can use:

```{r}
> hist(mygene, col="green", border=1, breaks=50, freq=F,  main = "Histogram of gene CCND3 expression values", ylab = "Density", xlab = "CCND3")
> lines(density(mygene), col="red", lwd=5)
```
![histplus](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/histplusline.png)

* To add vertical/horizontal lines to the plot, we use the command: **abline(v=)** or **abline(h=)**
* Add a vertical and horizontal line to the gene expression plot

```{r}
> par(mar=c(8,6,3,3))
> plot(mygene, type="l", lwd=5, col=1)
> mtext("CCND3 gene expression values", side=1, line=4.5)
> abline(v=20, col="red", lwd=3, lty=3)
> abline(h=1.5, col="green", lwd=5, lty=2)
```

![ABline](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/ABline.png)


* Add the diagonal lines to the plot: `abline(a= , b= )` 
* E.g.3: Add a diagonal line to the gene expression plot between patient 1 (ALL) and patient 38 (AML)

```{r}
> plot(golub[,1], golub[,38], ylab = "Expression of Sample 38", xlab = "Expression of Sample 1")
> abline(a=0, b=1, col="red", lwd=3)
```
![sample1and38](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/sample1and38.png)

### Box-plots 
*  Box plots can be used to observe the distribution of all gene expression values in all patients

```{r}
> par(mar=c(8,6,3,3))
> boxplot(golub, pch=18)
> mtext("Gene expression of 38 Samples", side=1, line=3)
```
![boxplotN](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/boxplotsN.png)


* You can also use a box plot for a subset of gene expression values
    + Observe the distribution of all gene expression values in patient 1, 10, 38. Plot them in different colors

```{r}
> par(mar=c(8,6,3,3))
> boxplot(golub[,1], golub[,10], golub[,38], border=1, col= c(2,3,4), pch=15, lwd=1)
> mtext("Gene expression distribution of 3 samples", side=1, line=4)
> axis(side=1, line=1.5, cex.axis=0.8, lwd=0, labels=c("Sample 1", "Sample 10", "Sample 38"), at=c(1,2,3))
```

![boxcolor](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/boxplotcolor.png)



### Computing simple test statistics

The **mt.teststat** and **mt.teststat.num.denum** functions of the **multtest** package provide a convenient way to compute test statistics for each row of a data frame, e.g., **two-sample Welch t-statistics, Wilcoxon statistics, F-statistics, paired t-statistics, and block F-statistics.**.

Let's compute the two-sample t-statistics that compare the gene expressions for each gene in the ALL and AML cases. This can be done with the mt.teststat function. The default test is the two-sample [Welch t-test](https://en.wikipedia.org/wiki/Welch%27s_t-test), it is more reliable when the two samples have unequal variances and unequal sample sizes.

```{r}
> teststat = mt.teststat(golub, golub.cl)
```

Now let's plot the t-test statistics:

```{r}
> install.packages('ggplot2') ## Another way to install a package
> library(ggplot2) ## ggplot2 is a plotting system for R, based on the grammar of graphics
> plt = ggplot(data.frame(teststat), aes(sample = teststat)) + stat_qq() + theme_grey()
> plt

```

![qqplot](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/qqplot.png)

The quantile-quantile (stat_qq) plot is a graphical technique used for determining if two data sets come from populations with a common distribution. Those points on the plot that look like outliers correspond to genes whose expression levels are different between the ALL and AML groups. We will learn how to detect statistically significant differentially expressed genes in the next lecture.


### Adjusting p-values

* The **mt.rawp2adj** function of the **multtest** package computes adjusted p-values for simple multiple testing procedures from a vector of raw (unadjusted) p-values. The procedures include the Bonferroni, Holm (1979), Hochberg (1988), and Sidak procedures for strong control of the family-wise Type I error rate **(FWER)**, and the Benjamini and Hochberg (1995) and Benjamini and Yekutieli (2001) procedures for (strong) control of the false discovery rate **(FDR)**. 

* First we will calculate the raw nominal two-sided p-values. For this data, we will assume that a standard normal distribution for the 3,051 test statistics.

```{r}
> rawp = 2 * (1 - pnorm(abs(teststat)))

```
We can adjust these p-values using the specified procedure as follows:

```{r}
> procedures = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
> adjusted = mt.rawp2adjp(rawp, procedures)
> adjusted$adjp[1:10, ]
             rawp   Bonferroni         Holm     Hochberg      SidakSS      SidakSD
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [3,] 8.881784e-16 2.709832e-12 2.708056e-12 2.708056e-12 2.709832e-12 2.708056e-12
 [4,] 1.332268e-15 4.064749e-12 4.060752e-12 4.060752e-12 4.064749e-12 4.060752e-12
 [5,] 1.554312e-15 4.742207e-12 4.735989e-12 4.735989e-12 4.742207e-12 4.735989e-12
 [6,] 4.418688e-14 1.348142e-10 1.345932e-10 1.345932e-10 1.348142e-10 1.345932e-10
 [7,] 1.207923e-13 3.685372e-10 3.678124e-10 3.678124e-10 3.685372e-10 3.678124e-10
 [8,] 2.009504e-13 6.130996e-10 6.116929e-10 6.116929e-10 6.130996e-10 6.116929e-10
 [9,] 2.606804e-13 7.953358e-10 7.932504e-10 7.932504e-10 7.953358e-10 7.932504e-10
[10,] 3.419487e-13 1.043285e-09 1.040208e-09 1.040208e-09 1.043285e-09 1.040208e-09
                BH           BY
 [1,] 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00
 [3,] 9.032775e-13 7.768732e-12
 [4,] 9.484413e-13 8.157168e-12
 [5,] 9.484413e-13 8.157168e-12
 [6,] 2.246903e-11 1.932472e-10
 [7,] 5.264817e-11 4.528061e-10
 [8,] 7.663745e-11 6.591283e-10
 [9,] 8.837064e-11 7.600409e-10
[10,] 1.043285e-10 8.972885e-10
```

------------------------------
### Basic plots practice at home

* Read the file: "NeuralStemCellData.tab" (the file is available at the *data* folder on the GitHub repository.

```{r}
> mydata<-read.delim( "https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/NeuralStemCellData.tab", row.names=1, header=T) 
> # mydata<-read.delim( "NeuralStemCellData.tab", row.names=1, header=T) #if you downloaded the file to your working directory
> head(mydata, n=7)
```
* The first four columns are corresponding to glioblastoma-derived (GD) neural stem-cells (column name start with 'G').
* The last two columns are from non-cancerous neural stem (NS) cells (column name start with 'CB').
 
### Tasks: 

1. Plot a Histogram and density line charts of the G144 cell
    + Distribution of gene expression level in the G144 cell
    
2. Scatter plots of last two NS cells (biological replicates)
    + Are the expression levels for the biological replicates (NS) consistent?
    
3. Box plot of the gene expression distribution in all 6 cells
    + The overall gene expression level distribution in all 6 cells
    
4. Bar plot of gene "BCAR1" expression level
    + Gene ("BCAR1") expression level in all 6 cells

### Do:

1. Histogram and density line charts of the G144 cell

```{r}
> par(mfrow=c(1,2), oma=c(0,0,0,0))
> par(mar=c(5,4,1,1))
> hist(log10(mydata$G144), breaks=30, xlab ='Histogram of G144', col="red", freq=T)
> d <- density(log10(mydata$G144))
> plot(d, col="blue", type="l", lwd=5, xlab ='Density of G144')
> abline(v=2, col="green", lwd=3)
```
![task1](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/task1.png)

2. Scatter plots of last two NS cells (biological replicates)

```{r}
> par(mar=c(5,4,1,1))
> plot(log10(mydata[,5]), log10(mydata[,6]), xlab ='NS replicate 1', ylab ='NS replicate 2')  
> abline(a=0, b=1, col="red", lty=2, lwd=5)
```

![task2](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/task2.png)


3. Box plot of the gene expression distribution in all 6 cells

```{r}
> par(mar=c(5,4,1,1))
> boxplot(log10(mydata), col=c(2:7), border="black", lwd=2, pch=23)
```
![tast3](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/task3.png)


4. Bar plot of gene "BCAR1" expression level

```{r}
> par(mar=c(4,4,1,1))
> mygene <- t(mydata[rownames(mydata)=="BCAR1" ,] )
> barplot(mygene, beside=T, col=c(2:7))
```
![task4](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/task4.png)

---------------

##4. Exportting polts as an image

* R function: **pdf()** and **png()**:

Export as pdf: 

```{r}
> pdf(file="myplot.pdf", width=11, height=8)
> hist(log10(mydata$G144), main = '', xlab ='Histogram of G144', breaks=50, col="red")
> graphics.off()
```

Export as png: 

```{r}
> png(file="myplot.png", width=800, height=680, units="px")
> plot(d, col="blue", main = '', xlab ='Density of G144', type="l", lwd=5)
> graphics.off()
```
 

### Saving multiple plots on one page:  

```{r}
> par(mfrow=c(1,3), mar=c(5,3,2,1))
> hist(log10(mydata$G144), main = '', xlab ='Histogram of G144', breaks=50, col="red")
> plot(d, col="blue", main = '', xlab ='Density of G144', type="l", lwd=5)
> boxplot(log10(mydata))
> graphics.off()
```

Using the data file: NeuralStemCellData.tab

* In only one page, plot:
  -  A Histogram of sample G166
  -  A density plot of sample G166
  -  The expression of gene TP53 across all samples 
  -  The Pearson correlation of Biological Replicates for neural stem (NS) cells 
 
## Week 5 Homework: :house: (Graded!)

- bam_stat.py is a utility from the [RSeQC] (http://rseqc.sourceforge.net/#) package, used to summarize the mapping statistics of a BAM or SAM file.
- To download large files from the web to CRI's HPC, connect to CRI's **cri-syncmon** server:

```
ssh t.cri.biowksp01@cri-syncmon.cri.uchicago.edu
```

1. Download the sequencing files from run ERR030885 from the Illumina bodyMap2 project, available at https://www.ebi.ac.uk/ena/data/view/ERR030885 (Note this is a PE library)
2. Develop an analysis pipeline to perform the QC, mapping using Bowtie2 and the summary of mapping statistics using bam_stat.py for this sample.
        - Report following metrics:
                - Total records:
                - Reads map to '+':
                - Reads mapped in proper pairs:
3. Using a Genome Browser, visualize and screen capture the reads mapping for the following genes: TP53, TNF and APOE.

* Submit all your scripts and screenshots via email.





