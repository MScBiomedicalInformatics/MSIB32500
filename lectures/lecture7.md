# Differential expression ananlysis of RNA-seq data with edgeR
**Center for Research Informatics, University of Chicago**

Saturdays 05/18/2019; 9:00 AM - 12:00 PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Develop a working knowledge of how to use R and Bioconductor for DEGs analysis of RNA-seq data using **edgeR** package

## 1. Introduction

* RNA sequencing (RNA-seq) technology has become a widely used technology for profiling transcriptional activity in biological systems. One of the most common aims of RNA-seq profiling is to identify genes or molecular
pathways that are differentially expressed (DE) between two or more biological conditions. Changes in expression can then be associated with differences in biology.
* In this lecture you will learn how to use the [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) package for detecting statistically significant differentially expressed genes using data from a RNA-Seq experiment. 
* The edgeR package, implements a range of statistical tests based on the **negative binomial distributions** expected on data from RNA-seq, including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests. The edgeR package can be also used to estimate differential **signal analysis** of other types of genomic data that produce **counts**, including ChIP-seq, Bisulfite-seq, Serial analysis of gene expression (SAGE) and Cap analysis gene expression (CAGE).
* The methods used in edgeR **do NOT support FPKM, RPKM normalizations**. 
* A user manual for the edgeR package is available [here](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

## 2. The data

* For this tutorial we will use the dataset available at the GEO
repository as series [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450). 
The goal of the experiemnt was to identify genes that are specifically expressed in lactating mammary glands, to that aim, the gene expression profiles of luminal and basal cells from different developmental stages (virgin, 18.5 day pregnant and 2 day lactating mice) were compared. Two biological replicates were collected for
each group.

* Sequencing was performed in a Illumina Hiseq sequencer, 30 million 100bp single-end reads for each sample.
* **Subread** version 1.4.4 was used to align the reads to the mouse mm10 genome and **featureCounts** was used to assign reads to *Entrez Genes* using RefSeq gene annotation. The raw FASTQ files are available at the Sequence Read Archive (SRA) repository.

### Sample description table:

The table below describes the experimental design, basal and luminal cell types are abbreviated with
B and L respectively. The GEO and SRA IDs for each RNA sample is also shown.

Define a working directory:

```
> setwd("/Users/jandrade/Desktop/w7")
> getwd()
```

For your convenience, you can read the *targets.txt* file directly from the GitHub repository:

```r
> targets <- read.delim(file = "https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/targets.txt", header = TRUE)

> targets
               GEO        SRA CellType    Status
MCL1.DG GSM1480297 SRR1552450        B    virgin
MCL1.DH GSM1480298 SRR1552451        B    virgin
MCL1.DI GSM1480299 SRR1552452        B  pregnant
MCL1.DJ GSM1480300 SRR1552453        B  pregnant
MCL1.DK GSM1480301 SRR1552454        B lactating
MCL1.DL GSM1480302 SRR1552455        B lactating
MCL1.LA GSM1480291 SRR1552444        L    virgin
MCL1.LB GSM1480292 SRR1552445        L    virgin
MCL1.LC GSM1480293 SRR1552446        L  pregnant
MCL1.LD GSM1480294 SRR1552447        L  pregnant
MCL1.LE GSM1480295 SRR1552448        L lactating
MCL1.LF GSM1480296 SRR1552449        L lactating
```

For later use, we combine the treatment factors into a single grouping factor:

```r
> group <- paste(targets$CellType, targets$Status, sep=".")
> group <- factor(group)
> table(group)
group
B.lactating  B.pregnant    B.virgin L.lactating  L.pregnant    L.virgin 
          2           2           2           2           2           2 
```
## 3. Preliminary analysis

* Download the genewise read counts for the GEO series GSE60450. 
* The file **GSE60450_Lactation-GenewiseCounts.txt.gz** will be downloaded to your working directory

```r
> FileURL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450","format=file", "file=GSE60450_Lactation-GenewiseCounts.txt.gz", sep="&")

> download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")

> GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz", row.names="EntrezGeneID")
> head(GenewiseCounts)
```
We will use the *substring* function to fromat the name of the columns

```r
> colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)
> dim(GenewiseCounts)
[1] 27179    13
```

```r
> head(GenewiseCounts)
          Length MCL1.DG MCL1.DH MCL1.DI MCL1.DJ MCL1.DK MCL1.DL MCL1.LA MCL1.LB MCL1.LC MCL1.LD MCL1.LE MCL1.LF
497097      3634     438     300      65     237     354     287       0       0       0       0       0       0
100503874   3259       1       0       1       1       0       4       0       0       0       0       0       0
100038431   1634       0       0       0       0       0       0       0       0       0       0       0       0
19888       9747       1       1       0       0       0       0      10       3      10       2       0       0
20671       3130     106     182      82     105      43      82      16      25      18       8       3      10
27395       4203     309     234     337     300     290     270     560     464     489     328     307     342

```
The row names of *GenewiseCounts *are the *Entrez Gene* Identifiers. The first column contains the length of each gene, the remaining 12 columns contain read counts for each sample.
 
Loading the The edgeR package:

```{r}
> library(edgeR)
```

The edgeR package stores data in a simple list-based data object called a *DGEList.*

```{r}
> y <- DGEList(GenewiseCounts[,-1], group=group, genes=GenewiseCounts[,1,drop=FALSE])
> y
An object of class "DGEList"
$counts
          MCL1.DG MCL1.DH MCL1.DI MCL1.DJ MCL1.DK MCL1.DL MCL1.LA MCL1.LB MCL1.LC MCL1.LD MCL1.LE MCL1.LF
497097        438     300      65     237     354     287       0       0       0       0       0       0
100503874       1       0       1       1       0       4       0       0       0       0       0       0
100038431       0       0       0       0       0       0       0       0       0       0       0       0
19888           1       1       0       0       0       0      10       3      10       2       0       0
20671         106     182      82     105      43      82      16      25      18       8       3      10
27174 more rows ...

$samples
              group lib.size norm.factors
MCL1.DG    B.virgin 23227641            1
MCL1.DH    B.virgin 21777891            1
MCL1.DI  B.pregnant 24100765            1
MCL1.DJ  B.pregnant 22665371            1
MCL1.DK B.lactating 21529331            1
7 more rows ...

$genes
          Length
497097      3634
100503874   3259
100038431   1634
19888       9747
20671       3130
27174 more rows ...
```
We can access each commponent of the DGEList object by using the *$* as follow.

```{r}
> options(digits=3)
> y$samples
              group lib.size norm.factors
MCL1.DG    B.virgin 23227641            1
MCL1.DH    B.virgin 21777891            1
MCL1.DI  B.pregnant 24100765            1
MCL1.DJ  B.pregnant 22665371            1
MCL1.DK B.lactating 21529331            1
MCL1.DL B.lactating 20015386            1
MCL1.LA    L.virgin 20392113            1
MCL1.LB    L.virgin 21708152            1
MCL1.LC  L.pregnant 22241607            1
MCL1.LD  L.pregnant 21988240            1
MCL1.LE L.lactating 24723827            1
MCL1.LF L.lactating 24657293            1
```

## 4. Adding gene annotation

We will use the **org.Mm.eg.db** package to link the *Entrez Gene Ids* to gene information in the NCBI database. The code below will add a column with the corresponding *gene symbols* to y$genes:

```{r}
> library(org.Mm.eg.db)
> y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
> head(y$genes)
          Length  Symbol
497097      3634    Xkr4
100503874   3259 Gm19938
100038431   1634 Gm10568
19888       9747     Rp1
20671       3130   Sox17
27395       4203  Mrpl15
```
To see a list of all available annotations on the org.Mm.eg.db package use:

```{r}
> columns(org.Mm.eg.db)
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
[11] "GO"           "GOALL"        "IPI"          "MGI"          "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"     
[21] "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"

```
## 5. Filtering to remove low counts

We will now remove Genes that have very low counts across all the libraries, such genes can be removed from the analysis without any loss of information, this is justified on both biological and statistical grounds.

* As a rule of thumb, we require that a gene have a count of at least 10 in at least some libraries.
* For the current analysis, we keep genes that have count-per-million (CPM) values **above 0.5 in at least two libraries.** (the method used in the edgeR vignette is to keep only those genes that have at least 1 CPM in at least 3 samples).

```{r}
> keep <- rowSums(cpm(y) > 0.5) >= 2
> table(keep)
FALSE  TRUE 
11375 15804 
```
* Other sensible filtering criteria are also possible. For example *keep <-rowSums(y$counts) > 50* is a very simple criterion that would keep genes with a total read count of more than 50. The exact
value is not important because the downstream differential expression analysis is not sensitive to the small changes in this parameter.

Try:

```{r}
> keep1 <- <-rowSums(y$counts) > 50
> table(keep1)
keep1
FALSE  TRUE 
11090 16089 
```

We will now subsetted the *DGEList object* to retain only the non-filtered genes:

```{r}
> y <- y[keep, , keep.lib.sizes=FALSE]
```
The option *keep.lib.sizes=FALSE* causes the library sizes to be *recomputed after the filtering*. This is generally recommended, although the effect on the downstream analysis is usually small.

## 6. Normalization for composition bias

Normalization by [trimmed mean of M values (TMM)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) is performed by using the calcNormFactors function, which returns the DGEList argument with only the norm.factors changed. It calculates a set of normalization factors, one for each sample, **to eliminate composition biases between libraries.** The product of these factors and the library sizes defines the *effective library size*, which replaces the original library size in all downstream analyses.

```{r}
> y <- calcNormFactors(y)
> y$samples
              group lib.size norm.factors
MCL1.DG    B.virgin 23218026        1.237
MCL1.DH    B.virgin 21768136        1.214
MCL1.DI  B.pregnant 24091588        1.126
MCL1.DJ  B.pregnant 22656713        1.070
MCL1.DK B.lactating 21522033        1.036
MCL1.DL B.lactating 20008326        1.087
MCL1.LA    L.virgin 20384562        1.368
MCL1.LB    L.virgin 21698793        1.365
MCL1.LC  L.pregnant 22235847        1.005
MCL1.LD  L.pregnant 21982745        0.923
MCL1.LE L.lactating 24719697        0.529
MCL1.LF L.lactating 24652963        0.535
```

* A normalization factor *below one* indicates that *a small number of high count genes are monopolizing the sequencing*, causing the counts for other genes to be lower than would be usual given the library size. As a result, *the effective library size will be scaled down for that sample*, (see Robinson and Oshlack (2010) on the the weighted trimmed mean of M-values  from the delta method on Binomial data) 

* Here we see that the luminal-lactating samples have low normalization factors. This is a sign that these samples contain a number of very highly upregulated genes.

## 7. Multi-Dimensional Scaling (MDS) plots

We will now explore the overall differences between the expression profiles of the different
samples, and weather samples are grouping together as expected. 

```{r}
> pch <- c(0,1,2,15,16,17)
> colors <- rep(c("darkgreen", "red", "blue"), 2)
> plotMDS(y, col=colors[group], pch=pch[group])
> legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)
```
![MDSweek7](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDSplotW7.png)

* As we can see, replicate samples from the same group cluster together while samples from different groups are well separated. 
* Samples are separated by the *cell type* in the first dimension, and by the *mouse status* in the second dimension.

## 8 Exploring the expression profiles of individual samples

The expression profiles of individual samples can be explored more closely with **mean-difference (MD) plots.** 

An MD plot visualizes the *library size-adjusted log-fold change* between two libraries *(the difference)* against *the average* log-expression across those libraries *(the mean)*.

*  To compares **sample 1** to an artificial reference library constructed from **the average of all the other samples** we can use:

```{r}
> plotMD(y, column=1)
> abline(h=0, col="red", lty=2, lwd=2)
```

![MDplotS1](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDplotS1.png)

The bulk of the genes are centered around the line of zero log-fold change. The diagonal lines in the lower left of the plot correspond to genes with counts of 0, 1, 2 and so on in the first sample.

Let's now look at one of the **luminal lactating** samples that were observed have low normalization factors:

```{r}
> plotMD(y, column=11)
> abline(h=0, col="red", lty=2, lwd=2)
```

![MDplotS11](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDplotS11.png)


For this sample, the log-ratios show noticeable positive skew, with a number of very highly upregulated genes. In particular, there are a number of points in the upper right of the plot, corresponding to genes that are both **highly expressed and highly up-regulated** in this sample compared to others.

## 9. Differential Expression Analysis

The differential expression analysis in **edgeR** requires a **design matrix** to be specified. The design matrix describes which samples belong to which groups, and it also defines how the experimental effects are parametrized in the linear models. 

```{r}
> design <- model.matrix(~0+group)
> colnames(design) <- levels(group)
> design
   B.lactating B.pregnant B.virgin L.lactating L.pregnant L.virgin
1            0          0        1           0          0        0
2            0          0        1           0          0        0
3            0          1        0           0          0        0
4            0          1        0           0          0        0
5            1          0        0           0          0        0
6            1          0        0           0          0        0
7            0          0        0           0          0        1
8            0          0        0           0          0        1
9            0          0        0           0          1        0
10           0          0        0           0          1        0
11           0          0        0           1          0        0
12           0          0        0           1          0        0
attr(,"assign")
[1] 1 1 1 1 1 1
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"
```

### Dispersion estimation

As we know, **edgeR** uses the **negative binomial (NB) distribution** to model the read counts for each gene in each sample. The **dispersion** parameter of the NB distribution **accounts for variability between biological replicates** [see McCarthy DJ, et al, Nucleic Acids Res. 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378882/)

Dispersion estimates are easily obtained from the **estimateDisp** function:

```{r}
> install.packages("statmod")
> library(statmod)

> y <- estimateDisp(y, design, robust=TRUE)
```
This returns a DGEList object with additional components (common.dispersion, trended.dispersion
and tagwise.dispersion) added to hold the estimated dispersions.

We can now plot the dispersion estimates with the **plotBCV** function:

```{r}
> plotBCV(y)
```
![poltdisp](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/plotDisp.png)

The NB dispersions tend to be higher for genes with very low counts. The *biological coefficient
of variation (BCV)* tends to be in range from **0.05 to 0.2** for genetically identical mice or
cell lines, whereas somewhat larger values (> 0.3) are observed for human subjects.

* The NB model can be extended with **quasi-likelihood (QL)** methods to account for **gene-specific variability** from both biological and technical sources. 
* Under the QL framework, the NB dispersion trend is used to describe the overall biological variability across all genes, and gene-specific variability above and below the overall level is picked up by the QL dispersion. In the QL approach, the individual (tagwise) NB dispersions are not used.

The estimation of QL dispersions is performed using the **glmQLFit** function:

```{r}
> fit <- glmQLFit(y, design, robust=TRUE)
> head(fit$coefficients)
```
The QL dispersions can be visualized by **plotQLDisp**

```{r}
> plotQLDisp(fit)
```

![QLdisp](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/QLdisp.png)

The QL functions moderate the genewise QL dispersion estimates in the same way that the limma package moderates variances. **The raw QL dispersion** estimates are **squeezed towards a global trend**, and this moderation **reduces the uncertainty of the estimates and improves testing power**.



### Testing for Differential Expression

We will now test for Differential Expression between the experimental groups. For this tutorial we will detect the differences between the **basal pregnant** and **basal lactating** groups. Other contrast can be easily constructed with the **makeContrasts** function:

```{r}
> B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
```
We will use QLF-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty in dispersion estimation.

```{r}
> res <- glmQLFTest(fit, contrast=B.LvsP)
```

The top DE genes can be viewed with the **topTags** function, the (FDR), multiple testing correction is performed using the Benjamini-Hochberg method.

```{r}
> topTags(res)
Coefficient:  1*B.lactating -1*B.pregnant 
       Length  Symbol logFC logCPM   F   PValue      FDR
12992     765 Csn1s2b  6.09  10.18 421 4.86e-11 7.68e-07
211577   2006  Mrgprf  5.15   2.74 345 1.30e-10 8.05e-07
226101   7094    Myof  2.32   6.44 322 1.97e-10 8.05e-07
381290   8292  Atp2b4  2.14   6.14 320 2.04e-10 8.05e-07
140474  11281    Muc4 -7.17   6.05 308 2.64e-10 8.34e-07
231830   3346 Micall2 -2.25   5.18 282 4.47e-10 1.18e-06
24117    2242    Wif1 -1.82   6.76 260 7.30e-10 1.65e-06
12740    1812   Cldn4 -5.32   9.87 298 8.88e-10 1.71e-06
21953     667   Tnni2  5.75   3.86 313 9.76e-10 1.71e-06
231991   2873   Creb5  2.57   4.87 241 1.16e-09 1.83e-06
```

* The top DE gene **Csn1s2b** has a large positive logFC, showing that it is far more highly expressed in the **basal cells** of **lactating** than **pregnant** mice.

* The total number of DE genes identified at an **FDR of 5%** can be shown with **decideTestsDGE**. There are in fact more than **5000** DE genes in this comparison:

```{r}
> is.de <- decideTestsDGE(res)
> summary(is.de)
       1*B.lactating -1*B.pregnant
Down                          2770
NotSig                       10529
Up                            2505
```

The magnitude of the differential expression changes can be visualized with a fitted model MD plot:

```{r}
> plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
```
![MDdes](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDdegs.png)

* The logFC for each gene is plotted against the average abundance, genes that are significantly DE are highlighted in red and blue.

To narrow down the list to genes that are more biologically meaningful, modify the statistical test so as to detect expression changes greater than a specified threshold, using the **glmTreat** function. 

```{r}
> tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
> topTags(tr)
Coefficient:  1*B.lactating -1*B.pregnant 
       Length  Symbol logFC unshrunk.logFC logCPM   PValue      FDR
12992     765 Csn1s2b  6.09           6.09  10.18 6.36e-11 1.00e-06
211577   2006  Mrgprf  5.15           5.15   2.74 1.86e-10 1.47e-06
140474  11281    Muc4 -7.17          -7.34   6.05 4.26e-10 2.25e-06
226101   7094    Myof  2.32           2.32   6.44 9.85e-10 2.94e-06
12740    1812   Cldn4 -5.32          -5.32   9.87 1.18e-09 2.94e-06
21953     667   Tnni2  5.75           5.76   3.86 1.23e-09 2.94e-06
381290   8292  Atp2b4  2.14           2.14   6.14 1.30e-09 2.94e-06
231830   3346 Micall2 -2.25          -2.25   5.18 2.45e-09 4.85e-06
231991   2873   Creb5  2.57           2.57   4.87 4.27e-09 7.50e-06
16012    1289  Igfbp6  2.87           2.87   3.67 5.38e-09 8.50e-06
```
To see how many genes passed this threshold we can use:

```{r}
> is.de <- decideTestsDGE(tr)
> summary(is.de)
       1*B.lactating -1*B.pregnant
Down                           677
NotSig                       14761
Up                             366
```

About **1100** genes are detected as DE with a FC significantly above 1.5 at an FDR cut-off of 5%.

We can now visualized these genes in an MD plot:

```{r}
> plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
```

![MDgesth](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MDdegsth.png)


## 10. Heat map clustering

Heatmaps are a popular way to display differential expression results for publication purposes. To create a heatmap, we first need to convert the read counts into **log2-counts-per-million (logCPM) values**. This can be done with the **cpm** function:

```{r}
> logCPM <- cpm(y, prior.count=2, log=TRUE)
> rownames(logCPM) <- y$genes$Symbol
> colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
```

We will create a heatmap to visualize the top 30 DE genes between **B.lactating** and
**B.pregnant**.  First we select the logCPM values for the 30 top genes:

```{r}
> o <- order(tr$table$PValue)
> logCPM <- logCPM[o[1:30],]
```
Then we scale each row (each gene) to have mean zero and standard deviation one. This scaling is commonly done for heatmaps and ensures that the heatmap displays relative changes for each gene:

```{r}
> logCPM <- t(scale(t(logCPM)))
```
Then we can create the heat map using the **heatmap.2** function in the gplots package:

```{r}
> library(gplots)
> col.pan <- colorpanel(100, "blue", "white", "red")
> heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
```
![heatmapDEGs](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/heatmapDEGs.png)

By default, **heatmap.2** clusters genes and samples based on the *Euclidean distance* between the expressionvalues. Because we have pre-standardized the rows of the logCPM matrix, the Euclidean distance between each pair of genes is proportional to (1 − r)^2, where r is the *Pearson correlation coefficient* between the two genes (the heatmap will cluster together genes that have positively correlated logCPM values).

## 11. Pathway analysis

### Gene ontology (GO) analysis

* A simple and effective way to interpret the list of DEGs is to count the number of DEGs that are annotated with each possible GO term. GO terms that occur frequently in the list of DE genes are said to be **over-represented** or **enriched**.

* In edgeR, GO analyses can be conveniently conducted using the **goana** function:

```{r}
> go <- goana(tr, species="Mm")
> topGO(go, n=15)
                                           Term Ont    N Up Down  P.Up   P.Down
GO:0007059               chromosome segregation  BP  281  1   55 0.999 1.62e-21
GO:0000070 mitotic sister chromatid segregation  BP  142  0   39 1.000 4.36e-21
GO:0000819         sister chromatid segregation  BP  165  0   41 1.000 2.51e-20
GO:1903047           mitotic cell cycle process  BP  599  7   80 0.988 9.07e-20
GO:0000280                     nuclear division  BP  352  6   58 0.834 7.67e-19
GO:0098813       nuclear chromosome segregation  BP  219  1   45 0.995 1.15e-18
GO:0000776                          kinetochore  CC  126  1   34 0.950 2.80e-18
GO:0140014             mitotic nuclear division  BP  245  2   47 0.980 3.57e-18
GO:0000775       chromosome, centromeric region  CC  181  1   40 0.986 6.68e-18
GO:0051301                        cell division  BP  538  5   71 0.996 2.44e-17
GO:0022402                   cell cycle process  BP  957 15  100 0.965 7.60e-17
GO:0000278                   mitotic cell cycle  BP  769  8   87 0.998 8.56e-17
GO:0048285                    organelle fission  BP  400  9   58 0.592 3.43e-16
GO:0042254                  ribosome biogenesis  BP  267  0   46 1.000 6.29e-16
GO:0007049                           cell cycle  BP 1451 19  129 0.999 7.08e-16
```
* The row names of the output are the universal identifiers of the **GO terms** and the **Term** column gives the human-readable names of the terms. 
* The **Ont** column shows the ontology domain that each GO term belongs to. The three domains are: **biological process (BP)**, **cellular component (CC)** and **molecular function (MF)**. 
* The **N** column represents the total number of genes annotated with each GO term. 
* The **Up** and **Down** columns indicate the number of genes within the GO term that are significantly up- and down-regulated in this differential expression comparison,respectively. 
* The **P.Up** and **P.Down** columns contain the p-values for over-representation of the GO term in the up and down-regulated genes, respectively.

The goana function uses the NCBI RefSeq annotation and requires the use of Entrez Gene Identifiers.


### KEGG pathway analysis

Another popular annotation database is the **Kyoto Encyclopedia of Genes and Genomes (KEGG)**. Much smaller than GO, this is a curated database of molecular pathways and disease signatures. A KEGG analysis can be done exactly as for GO, but using the **kegga** function:

```{r}
> keg <- kegga(tr, species="Mm")
> topKEGG(keg, n=15, truncate=34)
                                         Pathway   N Up Down     P.Up   P.Down
path:mmu03008 Ribosome biogenesis in eukaryot...  76  0   19 1.00e+00 3.04e-10
path:mmu04110                         Cell cycle 118  0   20 1.00e+00 1.29e-07
path:mmu05150    Staphylococcus aureus infection  31  0   10 1.00e+00 3.81e-07
path:mmu00240              Pyrimidine metabolism  93  0   15 1.00e+00 8.98e-06
path:mmu00100               Steroid biosynthesis  18  5    0 4.33e-05 1.00e+00
path:mmu00900    Terpenoid backbone biosynthesis  21  5    0 9.71e-05 1.00e+00
path:mmu04970                 Salivary secretion  63  8    7 9.77e-05 1.78e-02
path:mmu04972               Pancreatic secretion  66  8    4 1.36e-04 3.13e-01
path:mmu04610 Complement and coagulation casc...  48  2    9 3.06e-01 1.72e-04
path:mmu04925 Aldosterone synthesis and secre...  75  8    6 3.34e-04 1.02e-01
path:mmu00230                  Purine metabolism 152  2   16 8.71e-01 8.19e-04
path:mmu04640         Hematopoietic cell lineage  60  2    9 4.06e-01 9.61e-04
path:mmu04514     Cell adhesion molecules (CAMs) 112  2   13 7.36e-01 1.00e-03
path:mmu05144                            Malaria  38  2    7 2.20e-01 1.01e-03
path:mmu00983    Drug metabolism - other enzymes  51  2    8 3.31e-01 1.36e-03
```
* By default, the kegga function automatically reads the latest KEGG annotation from the Internet each time it is run.



