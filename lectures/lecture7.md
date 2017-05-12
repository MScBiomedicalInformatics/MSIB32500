
# Analysis of gene expression data from RNA-seq experiemnts with edgeR
**Center for Research Informatics, University of Chicago**

April - June 2017; Saturdays 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Develop a working knowledge on how to use R and Bioconductor for analysis of RNA-seq data using **edgeR** package

## 1. Introduction

* In this lecture you will learn how to use the [edgeR package in R](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) for detecting statistically significan differentially expressed genes using data from an RNA-Seq experiemnt. 
* A user manual for the edgeR package is available [here](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)


### What is edgeR?

edgeR is an R package that performs differential tag/gene expression using count data under a negative binomial model. The methods used in edgeR *do NOT support FPKM, RPKM* or other types of data that are not reads counts.


## 2. Installing the package

* We will start by installing the edgeR and baySeq packages:

```{r}
> source("http://bioconductor.org/biocLite.R")
> biocLite("edgeR")
> library(edgeR)
```
## 3. The dataset


* For this tutorial we will use the RNA-Seq data set from the paper by Li et al. [2008] entitled: ”Determination of tag density required for digital transcriptome analysis: Application to an androgen-sensitive prostate cancer model.” The paper is available [here.](https://www.ncbi.nlm.nih.gov/pubmed/19088194) 
* This study compared *hormone treated cancer cells* to *non-treated cancer cells.* For the non-treated, there are *4 biological replicates.* For the treated, there are *3 biological replicates.*

* Download to your local computer the raw version of the dataset: 'expression.txt', from the *data* folder on the github, and load that data to your R working enviroment:

```{r}
> raw.data <- read.table( file = "expression.txt" , header = TRUE )
> head(raw.data)
		ensembl_ID lane1 lane2 lane3 lane4 lane5 lane6 lane8  len
1 ENSG00000215696     0     0     0     0     0     0     0  330
2 ENSG00000215700     0     0     0     0     0     0     0 2370
3 ENSG00000215699     0     0     0     0     0     0     0 1842
4 ENSG00000215784     0     0     0     0     0     0     0 2393
5 ENSG00000212914     0     0     0     0     0     0     0  384
6 ENSG00000212042     0     0     0     0     0     0     0   92
```

* edgeR requires the the dataset to contain only the counts with the row names as the gene ids and the column names as the sample IDs, so we will reformat the data:

```{r}
> counts <- raw.data[,2:8]
> rownames(counts) <- raw.data[ , 1 ] # gene names
> colnames(counts) <- paste(c(rep("C_R",4),rep("T_R",3)),c(1:4,1:3),sep="") # Adding sample names
> head(counts)
				 C_R1 C_R2 C_R3 C_R4 T_R1 T_R2 T_R3
ENSG00000215696    0    0    0    0    0    0    0
ENSG00000215700    0    0    0    0    0    0    0
ENSG00000215699    0    0    0    0    0    0    0
ENSG00000215784    0    0    0    0    0    0    0
ENSG00000212914    0    0    0    0    0    0    0
ENSG00000212042    0    0    0    0    0    0    0
```
Data summaries:

```{r}
> dim(counts)
[1] 37435     7

> colSums(counts)/1e06 # Library Sizes in millions of reads

C_R1     C_R2     C_R3     C_R4     T_R1     T_R2     T_R3 
0.978576 1.156844 1.442169 1.485604 1.823460 1.834335 0.681743


> table(rowSums(counts))[ 1:30 ] # Number of genes with low (<30) counts

    0     1     2     3     4     5     6     7     8     9    10    11    12    13 
15558  1826  1149   717   635   388   447   291   308   240   248   191   252   157 
   14    15    16    17    18    19    20    21    22    23    24    25    26    27 
  190   155   162   126   166    99   169   125   123   100   124   100    87   100 
   28    29 
  103    73 
```

## 4. Building the edgeR Object

* **DGEList()** is the function that coverts the count matrix into an edgeR object. First, we create a group variable that tells edgeR which samples belong to which group and supply that to DGEList in addition to the count matrix. We can then see the elements that the object contains by using the **names()** function. These elements can be accessed using the $ symbol.


```{r}
> group <- c(rep("C", 4) , rep("T", 3))
> group
[1] "C" "C" "C" "C" "T" "T" "T"

> cds <- DGEList(counts, group=group, remove.zeros=T)
> cds
An object of class "DGEList"
$counts
                C_R1 C_R2 C_R3 C_R4 T_R1 T_R2 T_R3
ENSG00000124208  478  619  628  744  483  716  240
ENSG00000182463   27   20   27   26   48   55   24
ENSG00000124201  180  218  293  275  373  301   88
ENSG00000124205    0    0    5    5    0    0    0
ENSG00000124207   76   80   85   97   80   81   37
21872 more rows ...

$samples
     group lib.size norm.factors
C_R1     C   978576            1
C_R2     C  1156844            1
C_R3     C  1442169            1
C_R4     C  1485604            1
T_R1     T  1823460            1
T_R2     T  1834335            1
T_R3     T   681743            1

> names(cds)
[1] "counts"  "samples"

> cds$counts
> cds$samples

```
* We now need to filter out low count reads since otherwise it would be impossible to detect differential expression. 
* The method used in the edgeR vignette is to keep only those genes that have *at least 1 read per million in at least 3 samples.* 
* Once this filtering is done, we can calculate the normalization factors which correct for the different compositions of the samples. * The effective library sizes are then the product of the actual library sizes and these factors.

```{r}
> cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
> dim(cds)
[1] 16494     7
```

```{r}
> cds <- calcNormFactors(cds)
> cds
An object of class "DGEList"
$counts
                C_R1 C_R2 C_R3 C_R4 T_R1 T_R2 T_R3
ENSG00000124208  478  619  628  744  483  716  240
ENSG00000182463   27   20   27   26   48   55   24
ENSG00000124201  180  218  293  275  373  301   88
ENSG00000124207   76   80   85   97   80   81   37
ENSG00000125835  132  200  200  228  280  204   52
16489 more rows ...

$samples
     group lib.size norm.factors
C_R1     C   978576    1.0293595
C_R2     C  1156844    1.0377231
C_R3     C  1442169    1.0357600
C_R4     C  1485604    1.0393389
T_R1     T  1823460    0.9534767
T_R2     T  1834335    0.9524461
T_R3     T   681743    0.9576010

> # effective library sizes
> cds$samples$lib.size * cds$samples$norm.factors

```

## 5. Multi-Dimensional Scaling (MDS) Plot

* An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions. Let's create an MDS plot for our samples:

```{r}
> plotMDS(cds , main="MDS Plot for Count Data", labels= colnames(cds$counts))

```
![mdsplot](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDSplotC.png)


## 6. Estimating Dispersions

We will now proceed to estimate the common dispersion values across all genes by the Conditional Maximum Likelihood **(CML)** method *see Robinson MD and Smyth GK 2008* [paper](http://biostatistics.oxfordjournals.org/content/9/2/321).

The **CML** method involves computing a matrix of quantile-quantile normalized counts, *called pseudo-counts.* The pseudo-counts are adjusted in such a way that the library sizes are equal for all samples, while preserving differences between groups and variability within each group. The pseudo-counts are included in the output of the function, but are intended mainly for internal use of edgeR.

```{r}
> cds <- estimateCommonDisp(cds)
> cds
An object of class "DGEList"
$counts
                C_R1 C_R2 C_R3 C_R4 T_R1 T_R2 T_R3
ENSG00000124208  478  619  628  744  483  716  240
ENSG00000182463   27   20   27   26   48   55   24
ENSG00000124201  180  218  293  275  373  301   88
ENSG00000124207   76   80   85   97   80   81   37
ENSG00000125835  132  200  200  228  280  204   52
16489 more rows ...

$samples
     group lib.size norm.factors
C_R1     C   978576    1.0293595
C_R2     C  1156844    1.0377231
C_R3     C  1442169    1.0357600
C_R4     C  1485604    1.0393389
T_R1     T  1823460    0.9534767
T_R2     T  1834335    0.9524461
T_R3     T   681743    0.9576010

$common.dispersion
[1] 0.02001527

$pseudo.counts
                     C_R1      C_R2      C_R3      C_R4      T_R1      T_R2      T_R3
ENSG00000124208 605.87002 658.20589 536.40077 615.27355 353.45120 524.28061 468.50458
ENSG00000182463  33.44225  21.34180  22.98983  21.29626  34.86219  40.20951  45.59685
ENSG00000124201 228.31842 231.87309 250.63614 227.27128 275.14468 219.80966 176.09331
ENSG00000124207  95.71876  85.06503  72.35275  80.09793  58.50832  58.97925  70.74059
ENSG00000125835 167.80525 212.52454 170.71896 188.58612 207.32975 148.80645 106.82496
16489 more rows ...

$pseudo.lib.size
[1] 1276768

$AveLogCPM
[1] 8.722936 4.689348 7.508801 5.897509 7.095842
16489 more elements ...
```
* Once the common dispersion is estimated we will proceed to estimate **the tagwise dispersions.** 
* In this scenario, each gene will get its own unique dispersion estimate, but the *common dispersion* is still used in the calculation.
*  The **tagwise dispersions are squeezed toward the common value.** The amount of squeezing is governed by the paramter **prior.n.** The higher prior.n, the closer the estimates will be to the common dispersion. 
* The recommended **prior.n.** value is the nearest integer to: 50/(#samples − #groups). For this data set that’s 50/(7 − 2) = **10.**

```{r}
> cds <- estimateTagwiseDisp(cds, prior.n = 10)
> names(cds)
> summary(cds$tagwise.dispersion)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.004014 0.016910 0.031670 0.070610 0.073990 1.281000
```

## 7. Mean-Variance Plot

* Now that we have estimated the dispersion parameters we can see how well they fit the data by plotting the **mean-variance relationship**

* Four things are shown in the plot: 
	- 	The raw variances of the counts (grey dots) 
	-  The variances using the tagwise dispersions (light blue dots)
	-  The variances using the common dispersion (solid blue line) 
	-  The variance = mean a.k.a. poisson variance (solid black line) 
	
* The plot function outputs the variances which will be stored in the data set meanVarPlot.


```{r}
> meanVarPlot <- plotMeanVar( cds , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

```

![meanvar](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/Mean-Variance.png)


## 8. Detecting Differentially Expressed Genes

* The function **exactTest()** performs pair-wise tests for differential expression between two groups of  Negative-Binomial Counts 
* The important parameters are **dispersion** which can have the following values: *common, tagwise, and auto*, and **pair** which indicates which two groups should be compared. 
* The output of **exactTest()** is a list of elements, one of which is a **table of the results.** 


```{r}

> de.common <- exactTest(cds, pair=c("C","T"), dispersion="common" )
> de.tagwise <- exactTest(cds, pair=c("C","T"), dispersion="tagwise" )
> de.auto <- exactTest(cds, pair=c("C","T"), dispersion="auto" )

> names(de.tagwise)
> de.tagwise$comparison # which groups have been compared
> head(de.tagwise$table) # results table 
> head(cds$counts)

```
* As you can see, the **table of results** from **exactTest()** contains 'raw' p-values, they are not adjusted for multiple testing. 
* To correct for multiple testing, we will use the function **topTags()**. It takes the output from exactTest(), adjusts the raw p-values using the **False Discovery Rate (FDR) correction**, and returns the **top differentially expressed genes.** * The output is similar to that of exactTest() but with **a column of adjusted p-values and sorted by increasing p-value.** 
* The *sort.by* argument allows you to sort the table by p-value, concentration or fold-change. 
* We can also use topTags() to return the original counts of the top differentially expressed genes, by setting the *n* parameter to the total number of genes, we can save the entire topTags() results table.

```{r}
# Top tags for tagwise analysis
> options(digits=3) # print only 3 digits
> topTags(de.tagwise, n = 20 , sort.by = "p.value" ) # top 20 DE genes
Comparison of groups:  T-C 
                logFC logCPM    PValue       FDR
ENSG00000151503  5.82   9.71  0.00e+00  0.00e+00
ENSG00000096060  5.01   9.95  0.00e+00  0.00e+00
ENSG00000166451  4.68   8.84 1.55e-249 8.51e-246
ENSG00000127954  8.12   7.21 1.53e-229 6.32e-226
ENSG00000162772  3.32   9.74 2.74e-204 9.02e-201
ENSG00000115648  2.60  11.47 5.03e-180 1.38e-176
ENSG00000116133  3.25   8.79 7.57e-175 1.78e-171
ENSG00000113594  4.11   8.05 1.11e-160 2.28e-157
ENSG00000130066  2.61   9.99 2.73e-154 5.01e-151
ENSG00000116285  4.20   7.36 7.52e-150 1.24e-146
ENSG00000140526  2.24  10.35 5.82e-143 8.73e-140
ENSG00000123983  3.66   8.59 1.36e-137 1.87e-134
ENSG00000140263  2.82   8.20 2.81e-132 3.57e-129
ENSG00000167751  2.95   8.75 6.16e-119 7.26e-116
ENSG00000067113  2.32   9.93 5.35e-115 5.89e-112
ENSG00000104419  3.25   7.42 1.46e-114 1.50e-111
ENSG00000117143  2.11   9.59 4.65e-112 4.51e-109
ENSG00000166086  5.48   6.28 2.40e-111 2.20e-108
ENSG00000175928 -4.56   6.92 1.30e-104 1.13e-101
ENSG00000064042  1.67  11.40 1.54e-100  1.27e-97
```

Store full topTags results table

```{r}
> resultsTbl.tgw <- topTags(de.tagwise, n=nrow(de.tagwise$table))$table
> head(resultsTbl.tgw)
```

* Now that we have the full results table with the adjusted p-values we can compare those p-values to a defined *significance level* and see how many differentially expressed genes we find. 
* Here we use a *significance level of 0.05.* In addition, we can use the function decideTestsDGE() to determine how many DE genes are *up* or *down* regulated compared to control.

```{r}
> de.genes.tgw <- rownames(resultsTbl.tgw)[resultsTbl.tgw$FDR <= 0.05 ]

# How many DEGs are significant?
> length(de.genes.tgw)
[1] 4431

# Percentage of total genes
> length(de.genes.tgw)/nrow(resultsTbl.tgw) * 100
[1] 26.9

# Up/Down regulated summary for tagwise results

> summary(decideTestsDGE(de.tagwise, p.value=0.05))
  [,1] 
-1  2086
0  12063
1   2345
```

## 9. Visualizing the DEG Results

* We will start by vizualing the spread of expression levels for the top DEGs. Below, we plot a histogram of log concentrations for the top 100 genes:

```{r}
> hist(resultsTbl.tgw[de.genes.tgw[1:100],"logCPM"], breaks=25, xlab="Log Concentration", col="blue",  freq=FALSE, main="Tagwise: Top 100" )

```
![tagwise100](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/tagwisetop100.png)

* The next plot is an MA plot that shows the relationship between concentration and fold-change across the genes
* The differentially expressed genes are colored in red and the non-differentially expressed are colored black. T
* The orange dots represent genes in which the counts *were zero in all samples of one of the groups.* 
* The blue line is added at a log-FC of 2 to represent a level for biological significance.

```{r}

> plotSmear(cds, de.tags=de.genes.poi, main="Poisson",
           pair = c("C","T"),
           cex = .35,
           xlab="Log Concentration", ylab="Log Fold-Change")
> abline(h = c(-2, 2) , col = "dodgerblue")
> plotSmear(cds, de.tags=de.genes.tgw, main="Tagwise",
           pair = c("C","T"),
           cex = .35,
           xlab="Log Concentration", ylab="Log Fold-Change")
> abline( h=c(-2,2), col="dodgerblue")

```

![maPlottag](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MAplotTag.png)

## 10. Saving the results


Now that the analysis is completed, we want to create an output table (or a csv file) containing all the results. We will pull together the concentrations, fold-changes, p-values, adjusted p-values, up/down regulated variable, dispersions, and the count matrix

```{r}
# Change column names to be specific to the analysis
> colnames(resultsTbl.tgw) <- c("logConc", "logFC", "pVal.Cmn","adj.pVal.Cmn")

# Re-order the count matrix to be in line with the order of the results.

> wh.rows.tgw <- match(rownames(resultsTbl.tgw), rownames(cds$counts))

> head(wh.rows.tgw)

> combResults.tgw <- cbind(resultsTbl.tgw,
                          "Tgw.Disp"=cds$tagwise.dispersion[wh.rows.tgw],
                          "UpDown.Tgw"=decideTestsDGE(de.tagwise, p.value = 0.05)[wh.rows.tgw],
                          cds$counts[wh.rows.tgw,])

> head(combResults.tgw)
                logConc logFC  pVal.Cmn adj.pVal.Cmn Tgw.Disp UpDown.Tgw C_R1 C_R2 C_R3 C_R4 T_R1 T_R2 T_R3
ENSG00000151503    5.82  9.71  0.00e+00     0.00e+00  0.00735          1   35   35   49   59 3307 3439 1224
ENSG00000096060    5.01  9.95  0.00e+00     0.00e+00  0.00682          1   65   79  105  113 3975 3727 1451
ENSG00000166451    4.68  8.84 1.55e-249    8.51e-246  0.00899          1   41   52   57   57 1750 1654  728
ENSG00000127954    8.12  7.21 1.53e-229    6.32e-226  0.00911          1    0    0    3    3  607  602  220
ENSG00000162772    3.32  9.74 2.74e-204    9.02e-201  0.00744          1  172  204  250  304 2972 3269 1112
ENSG00000115648    2.60 11.47 5.03e-180    1.38e-176  0.00623          1  940 1084 1317 1345 9730 9942 3272
```

Finally, we will save the complete results table on a file:

```{r}
> write.csv(combResults.tgw, file='DEG_list_edgeR.csv')
```


## Week 7 Homework: :house: 

* For this tutorial, we used the **tagwise** dispersion to estimate the DEGs (see section 8. Detecting Differentially Expressed Genes). Using the other available methods *common*, and *auto*, repeat the detection of the **top 100** statistically significant genes (DEGs). One list ofr each method.
* Create a list of the overlaping DEGs
* Send your homework via e-mail  


