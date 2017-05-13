
# Analysis of Illumina Microarray data with R and Bioconductor 

**Center for Research Informatics, University of Chicago**

April - June 2017; Saturdays 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

This hands-on tutorial is focused on the analysis of Illumina microarray data using R and Bioconductor, this tutorial assumes that you have previous experience using R for data analysis.

## 1. The Data:

- Expression profiling: Illumina MouseRef-8 v2.0 expression beadchip array was used to identify genes, which upon **8** or **24** hr of treatment with all-trans retinoic acid, display differential expression in **RARγ knockout (RARγKO)** murine embryonic stem cells relative to CCE WT cells. The study demonstrated that following **RA** treatment the majority of inducible transcripts are present at lower levels in RARγKO ES cells compared to WT ES cells.
- Overall experiment design: Murine embryonic stem cells **(WT and RARγKO)** were treated with either **all-trans retinoic acid** (up to 24 hr) or with **vehicle control (EtOH)**
- Data: The data set: GSE43221 used in this tutorial is available in NCBI’s gene expression omnibus.
- The paper is available in [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/23264745)

### Getting the data 


Please download and unzip the dataset: **“GSE43221_non-normalized_data.txt.gz”** from the following [link](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43221)

* Open the file **GSE43221_non-normalized_data.txt** file using Excel or a text editor and explore Illumina BeadStudio version 3.1.7 format.
* Observe that the first 6 lines are meatadata, the expression profile of the samples starts in line 7 with the header. Illumina Probe_ids are on rows, and samples on the columns.

## 2. Data analysis workflow:

For any bioinformatics analysis, it is always advisable to start with an analysis plan, this can be represented in the analysis workflow. For this tutorial, we will follow the analysis workflow bellow:

![workflow](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/workflow.png)


## 3. Installing requiered packages:
The first thing we will need to start working with R and Bioconductor, is to download and install the packages we will need. For this analysis we will use the Bioconductor **lumi** package. The **lumi** package provides an integrated solution for Illumina microarray data analysis, it includes functions for Illumina BeadStudio (GenomeStudio) data input, quality control, BeadArray-specific variance stabilization, normalization and gene annotation at the probe level. Available functions include: **lumiR(), lumiB(), lumiT(), lumiN() and lumiQ()** designed for data input, preprocessing and quality control. 

Downloading and installing **lumi** package:

```{r}
> source("http://bioconductor.org/biocLite.R")
> biocLite("lumi")
> library(lumi)
```

We will also need the corresponding **annotation libraries** **lumiMouseAll.db** and **AnnotationDbi**, we can download and install these libraries with the code below:

```{r}
> biocLite("lumiMouseAll.db")
> library(lumiMouseAll.db)

> biocLite("annotate")
> library(annotate)

> biocLite("AnnotationDbi")
> library(AnnotationDbi)
```
## 4. Reading the data

The function **lumiR()** supports directly reading of the Illumina Bead Studio toolkit output from version 1 to version 3. It can automatically detect the BeadStudio output version/format and create a new *LumiBatch* object for it. The **lumiR** function will automatically determine the starting line of the data; columns with header including AVG_Signal and BEAD_STD are required for the LumiBatch object. The sample IDs and sample labels are extracted from the column names of the data file. After reading the data, **lumiR** will *automatically* initialize the **QC slot** of the LumiBatch object by calling **lumiQ.**

To read the data, we will start by defining our working directory:

```{r}
> setwd("/Users/jorgeandrade/Desktop/GSE43221")
> getwd()

[1] "/Users/jorgeandrade/Desktop/GSE43221"
```
The folder **GSE43221** contains the **GSE43221_non-normalized_data.txt** file; we will now read/load the file into R using **lumiR**. The annotation file that we just installed **lumiMouseAll.db** is also loaded into the **featureData** of the *LumiBatch* object **x.lumi**

```{r}
> x.lumi <- lumiR('GSE43221_non-normalized_data.txt', lib='lumiMouseAll')
```
Next we will need to create a **.txt** file with the *phenotype information*, please create a tab-delimited file that list all the samples and associated phenotype as follow:

![pheno](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/pheno.png)

For your convenience, I have already created that file and it is available [here](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/Phenod.txt) Download this file to your working directory.

###Importing the phenotype data: 

```{r}
> phenod <- read.table("Phenod.txt", header=TRUE, row.names=1, sep="\t")
> phenod
```

Output:

```
                  Species Type Treatment Time Replicate
CCE_WT_Ctrl_1         CCE   WT      Ctrl Ctrl         1
CCE_WT_Ctrl_2         CCE   WT      Ctrl Ctrl         2
CCE_WT_Ctrl_3         CCE   WT      Ctrl Ctrl         3
CCE_WT_RA_1hr_1       CCE   WT        RA  1hr         1
CCE_WT_RA_1hr_2       CCE   WT        RA  1hr         2
CCE_WT_RA_1hr_3       CCE   WT        RA  1hr         3
CCE_WT_RA_8hr_1       CCE   WT        RA  8hr         1
CCE_WT_RA_8hr_2       CCE   WT        RA  8hr         2
CCE_WT_RA_8hr_3       CCE   WT        RA  8hr         3
CCE_WT_RA_24hr_1      CCE   WT        RA 24hr         1
CCE_WT_RA_24hr_2      CCE   WT        RA 24hr         2
CCE_WT_RA_24hr_3      CCE   WT        RA 24hr         3
RARE_KO_Ctrl_1       RARE   KO      Ctrl Ctrl         1
RARE_KO_Ctrl_2       RARE   KO      Ctrl Ctrl         2
RARE_KO_Ctrl_3       RARE   KO      Ctrl Ctrl         3
RARE_KO_RA_1hr_1     RARE   KO        RA  1hr         1
RARE_KO_RA_1hr_2     RARE   KO        RA  1hr         2
RARE_KO_RA_1hr_3     RARE   KO        RA  1hr         3
RARE_KO_RA_8hr_1     RARE   KO        RA  8hr         1
RARE_KO_RA_8hr_2     RARE   KO        RA  8hr         2
RARE_KO_RA_8hr_3     RARE   KO        RA  8hr         3
RARE_KO_RA_24hr_1    RARE   KO        RA 24hr         1
RARE_KO_RA_24hr_2    RARE   KO        RA 24hr         2
RARE_KO_RA_24hr_3    RARE   KO        RA 24hr         3
```
To get a data summary report we simply use:

```{r}
> x.lumi
```

Output:

```
Summary of data information:
	 Data File Information:
		Illumina Inc. BeadStudio version 3.1.7
		Normalization = none
		Array Content = MouseRef-8_V2_0_R0_11278551_A.bgx.xml
		Error Model = none
		DateTime = 8/8/2008 10:01 AM
		Local Settings = en-US
		ID_REF 	SAMPLE 1	CCE_WT_Ctrl_1	Detection Pval	SAMPLE 2	CCE_WT_Ctrl_2	Detection Pval	SAMPLE 3	CCE_WT_Ctrl_3	Detection Pval	SAMPLE 4	CCE_WT_RA_1hr_1	Detection Pval	SAMPLE 5	CCE_WT_RA_1hr_2	Detection Pval	SAMPLE 6	CCE_WT_RA_1hr_3	Detection Pval	SAMPLE 7	CCE_WT_RA_8hr_1	Detection Pval	SAMPLE 8	CCE_WT_RA_8hr_2	Detection Pval	SAMPLE 9	CCE_WT_RA_8hr_3	Detection Pval	SAMPLE 10	CCE_WT_RA_24hr_1	Detection Pval	SAMPLE 11	CCE_WT_RA_24hr_2	Detection Pval	SAMPLE 12	CCE_WT_RA_24hr_3	Detection Pval	SAMPLE 13	RARE_KO_Ctrl_1	Detection Pval	SAMPLE 14	RARE_KO_Ctrl_2	Detection Pval	SAMPLE 15	RARE_KO_Ctrl_3	Detection Pval	SAMPLE 16	RARE_KO_RA_1hr_1	Detection Pval	SAMPLE 17	RARE_KO_RA_1hr_2	Detection Pval	SAMPLE 18	RARE_KO_RA_1hr_3	Detection Pval	SAMPLE 19	RARE_KO_RA_8hr_1	Detection Pval	SAMPLE 20	RARE_KO_RA_8hr_2	Detection Pval	SAMPLE 21	RARE_KO_RA_8hr_3	Detection Pval	SAMPLE 22	RARE_KO_RA_24hr_1	Detection Pval	SAMPLE 23	RARE_KO_RA_24hr_2	Detection Pval	SAMPLE 24	RARE_KO_RA_24hr_3	Detection Pval

Major Operation History:
            submitted            finished
1 2013-03-19 09:31:37 2013-03-19 09:31:54
2 2013-03-19 09:31:54 2013-03-19 09:31:55
3 2013-03-19 09:31:55 2013-03-19 09:31:57
                                                                      command
1     lumiR("GSE43221_non-normalized_data.txt", lib.mapping = "lumiMouseAll")
2        lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)
3 addNuID2lumi(x.lumi = x.lumi, lib.mapping = lib.mapping, verbose = verbose)
  lumiVersion
1      2.10.0
2      2.10.0
3      2.10.0

Object Information:
LumiBatch (storageMode: lockedEnvironment)
assayData: 25697 features, 24 samples 
  element names: beadNum, detection, exprs, se.exprs 
protocolData: none
phenoData
  sampleNames: CCE_WT_Ctrl_1 CCE_WT_Ctrl_2 ... RARE_KO_RA_24hr_3 (24
    total)
  varLabels: sampleID
  varMetadata: labelDescription
featureData
  featureNames: 9lekf9_d_rteewJTpQ KqvaKDqYkgAXyjiImk ...
    Zl4J4B0PsjTJ1N19IM (25697 total)
  fvarLabels: Probe_Id Accession ... Definition (9 total)
  fvarMetadata: labelDescription
experimentData: use 'experimentData(object)'
Annotation:  
Control Data: N/A
QC information: Please run summary(x, 'QC') for details!
```
## 5. Preliminary exploratory analysis

We will now start an exploratory analysis of the data using different methods; this analysis will allow us to easily identify outliers and/or other problems with the samples.

We will start by producing a density/intensity histogram, and a boxplot of the distribution of intensities:

```{r}
> hist(x.lumi)
```
![hist](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/hist.png)

```{r}
> boxplot(x.lumi)
```
![box](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/boxplot.png)

####These two graphics are showing us: 

1. There is no evident outlier in the data set
2. The distribution of the intensities in each array (boxplot) illustrates the need for a normalization step

We will now explore the correlation of samples before normalization:

Pairwise correlation of TR control group: 

```{r}
> pairs(x.lumi[,c(1,2,3)], smoothScatter=TRUE)
> ## plot(x.lumi[,c(1,2,3)], what='pair')
```
![correlation](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/corr.png)

You can continue exploring the correlation of technical replicates (TR) for other groups by changing the values at: *pairs(x.lumi[,c(1,2,3)]*

Next we will create the [MA plot](https://en.wikipedia.org/wiki/MA_plot) for some samples of the LumiBatch object by using plot or **MAplot** functions. MA plots are used to determine if the data needs normalization and to test if the normalization worked. **M** is defined as the log intensity ratio (or difference between log intensities) and **A** is the average log intensity for a dot in the plot.
	
MA plots for samples in the control group before normalization:

```{r}
> MAplot(x.lumi[,c(1,2,3)])
> ## with smoothing > MAplot(x.lumi[,c(1,2,3)], smoothScatter=TRUE)
```
![MAplot](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MAplot.png)

## 6. Backgroud correction

The **lumi** package provides the function **lumiB()** for background correction. By default, the **BeadStudio output data** (imported as input for our analysis) is already background corrected.

### Variance-stabilizing transformation

Variance stabilization step is critical for subsequent statistical inference of differentially expressed genes.  The Variance-stabilizing transformation **(VST)** method, takes advantages of larger number of technical replicates available on the Illumina microarray.  The function **lumiT()** performs the variance-stabilizing transformation; both input and output are LumiBatch objects.

```{r}
> lumi.T <- lumiT(x.lumi)
```
Output:

```
Perform vst transformation ...
2013-03-21 14:39:15 , processing array  1 
2013-03-21 14:39:15 , processing array  2 
2013-03-21 14:39:15 , processing array  3 
2013-03-21 14:39:15 , processing array  4 
2013-03-21 14:39:15 , processing array  5 
2013-03-21 14:39:15 , processing array  6 
2013-03-21 14:39:15 , processing array  7 
2013-03-21 14:39:15 , processing array  8 
2013-03-21 14:39:15 , processing array  9 
2013-03-21 14:39:15 , processing array  10 
2013-03-21 14:39:15 , processing array  11 
2013-03-21 14:39:15 , processing array  12 
2013-03-21 14:39:15 , processing array  13 
2013-03-21 14:39:15 , processing array  14 
2013-03-21 14:39:15 , processing array  15 
2013-03-21 14:39:15 , processing array  16 
2013-03-21 14:39:16 , processing array  17 
2013-03-21 14:39:16 , processing array  18 
2013-03-21 14:39:16 , processing array  19 
2013-03-21 14:39:16 , processing array  20 
2013-03-21 14:39:16 , processing array  21 
2013-03-21 14:39:16 , processing array  22 
2013-03-21 14:39:16 , processing array  23 
2013-03-21 14:39:16 , processing array  24
```
```{r}
trans <- plotVST(lumi.T)
```
![trans](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/trans.png)

Now we will compare the VST and log2 transformations

```{r}
> matplot(log2(trans$untransformed), trans$transformed, type='l', main='compare VST and log2 transform', xlab='log2 transformed', ylab='VST transformed')
> abline(a=0, b=1, col=2)
```
![VST](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/VST.png)


## 7. Data Normalization

**lumi** package provides several normalization method options, which include *quantile normalization*, *SSN* (Simple Scaling Normalization), *RSN* (Robust Spline Normalization), loess normalization and Rank Invariant Normalization. Robust spline normalization *(RSN)* algorithm (which combines the features of quanitle and loess nor-malization), is designed to normalize the variance-stabilized data. We will use the function **lumiN()** to perform **RSN**, both input and output are LumiBatch objects. 

```{r}
> lumi.N <- lumiN(lumi.T, method='rsn') ## Normalization
> lumi.N.Q <- lumiQ(lumi.N)  ## QC of LumiBatch object
> summary(lumi.N.Q, 'QC') ## QC summary
```

We can now plot the density and box plot of our samples, after **RSN** normalization.

```{r}
> plot(lumi.N.Q, what='density')
> boxplot(lumi.N.Q)
```
Density After Normalization:
![densityAN](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/DensityAN.png)
Box Plot After Normalization:
![boxplotAN](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/boxplotAN.png)

Next we will plot the **MAplots** of samples in the control group after normalization:

```{r}
> MAplot(lumi.N.Q[,c(1,2,3)], smoothScatter=TRUE)
```
![MAplotAN](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/MAplotAN.png)

##8. Principal Component Analysis
Principal components analysis [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis) is a statistical technique for determining the key variables in a multidimensional data set that explain the differences in the observations, and can be used to simplify the analysis and visualization of multidimensional data sets. R has several  built-in functions for PCA analysis: **prcomp()** and **princomp()**.

```{r}
> ## note: the expression matrix needs to be transposed
> pca <- prcomp(t(exprs(x.lumi)), scale=TRUE)
> summary(pca)
```

Output:

```
Importance of components:
                            PC1      PC2      PC3     PC4      PC5     PC6      PC7      PC8     PC9     PC10
Standard deviation     143.1359 36.41579 24.26316 19.5683 19.04346 17.2648 16.35583 15.73560 14.9494 13.63279
Proportion of Variance   0.7973  0.05161  0.02291  0.0149  0.01411  0.0116  0.01041  0.00964  0.0087  0.00723
Cumulative Proportion    0.7973  0.84889  0.87180  0.8867  0.90082  0.9124  0.92283  0.93246  0.9412  0.94839
                           PC11     PC12     PC13     PC14     PC15     PC16    PC17    PC18    PC19    PC20
Standard deviation     13.28207 12.81103 12.27342 11.52524 10.92886 10.04001 9.86903 9.52206 9.41028 8.26338
Proportion of Variance  0.00687  0.00639  0.00586  0.00517  0.00465  0.00392 0.00379 0.00353 0.00345 0.00266
Cumulative Proportion   0.95526  0.96164  0.96750  0.97267  0.97732  0.98124 0.98503 0.98856 0.99201 0.99467
                          PC21    PC22   PC23      PC24
Standard deviation     7.53919 7.03579 5.5426 6.326e-14
Proportion of Variance 0.00221 0.00193 0.0012 0.000e+00
Cumulative Proportion  0.99688 0.99880 1.0000 1.000e+00

```

### Plotting the PCA

For this tutorial, we will use the library **rgl**, a 3D Visualization library that uses OpenGL.

* Note: MacOS user will need to download and install: [www.xquartz.org](https://www.xquartz.org)

```{r}
> install.packages("rgl")
> library("rgl")
```
Now we will create an array of different colors to distinguish samples that belong to the same group:

```{r}
> timeseries.colors <- c(rep("Green", 3), 
rep("DimGrey", 3), 
rep("IndianRed",3), 
rep("Red", 3), #WT samples
rep("BurlyWood", 3), 
rep("grey", 3), 
rep("pink", 3), 
rep("SkyBlue", 3))#KO samples

```

And finally we plot the PCA:

```{r}
> plot3d(pca$x[, 1:3], col=timeseries.colors, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
```
You should be able to see and interactive 3D PCA plot as bellow:

![PCA](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/PCA.png)

The PCA plot shows samples that belong to the same group in similar color. Note WT samples are represented in Red and KO samples in SkyBlue.


##9. Encapsulating the processing steps

* The function **lumiExpresso()** is designed to encapsulate the major functions of Illumina preprocessing. It is organized in a similar way as the function **expresso()** in the affy package. 

* The following code will perform the same processing as the previous multi-steps, and produced the same results as **lumi.N.Q()**.

```{r}
> lumi.N.Q1 <- lumiExpresso(x.lumi, normalize.param=list(method='rsn'))
```
At this point we can save the **“expression matrix”** into a tab delimited text file:

```{r}
> dataMatrix <- exprs(lumi.N.Q)
> dim(dataMatrix) ## Checking the dimension of my matrix
[1] 25697    24

> write.table(dataMatrix, file="Data-Matrix-NQ.txt", quote=F, sep = "\t")
```

##10.Data filtering

* We will now filter our LumiBatch object for probes with only **“good”** detectable signal (detection p-value smaller that 0.01). 
* Our filtering criteria will be to require detectable signal in at least 80% of our 24 samples:

```{r}
> presentCount <- detectionCall(x.lumi)
> selDataMatrix <- dataMatrix[presentCount > 19,]
> dim(selDataMatrix)
[1] 10243    24
```
Note we are filtering for signal present (presentCount) in at least 19 samples (80% of the total 24 samples). This filter reduced the number of rows from **25697**  before filtering to **10243** after filtering.

Let's now save the Expression Matrix after filtering

```{r}
> write.table(selDataMatrix,file="Data-Matrix-NQ-filtered.txt", quote=FALSE, sep = "\t")
```
We will then create a list of probeIDs:

```{r}
> probeList <- rownames(selDataMatrix)
```

## 11. Identifying Differentially Expressed Genes (DEGs)

For this tutorial, we will proceed to identify differentially expresed genes (DEG) only between control (WT, 3 samples) and 24 hr. of treatment (WT 24hr, 3 samples); subsequent pairwise DEG analysis can be repeated for different groups:

```{r}
> M6data <- selDataMatrix[, c(1,2,3,10,11,12)]
> dim(M6data)
 [1] 10243     6
```

```{r}
> head(M6data)

```

```

                   CCE_WT_Ctrl_1 CCE_WT_Ctrl_2 CCE_WT_Ctrl_3 CCE_WT_RA_24hr_1
WpaZ9x9f_hAnoR.VBE     11.785177     11.510133     11.537132        11.688847
ZhdXp75JftSF3iWLF4      9.845989      9.252357      9.430707         9.119103
rUCmU113ENlNRdxLSs      9.298117      8.836202      8.874029         8.892999
xnx7ijE2CcOSlY55eA      9.838060      9.430353      9.549075         9.581293
xPkUXaq37oR_vSi6HY     11.497402     11.320781     11.250264        11.626142
HSHd1YIsbjZ5TS3cTo     11.747510     11.563647     11.764771        12.125247
                   CCE_WT_RA_24hr_2 CCE_WT_RA_24hr_3
WpaZ9x9f_hAnoR.VBE        11.608809        11.429682
ZhdXp75JftSF3iWLF4         9.249300         9.118371
rUCmU113ENlNRdxLSs         8.902206         8.887582
xnx7ijE2CcOSlY55eA         9.539387         9.490084
xPkUXaq37oR_vSi6HY        11.647679        11.508672
HSHd1YIsbjZ5TS3cTo        12.033338        11.962688

```
Let us review our phenotype data:

```{r}
> phenod
```


```
Species Type Treatment Time Replicate
CCE_WT_Ctrl_1         CCE   WT      Ctrl Ctrl         1
CCE_WT_Ctrl_2         CCE   WT      Ctrl Ctrl         2
CCE_WT_Ctrl_3         CCE   WT      Ctrl Ctrl         3
CCE_WT_RA_1hr_1       CCE   WT        RA  1hr         1
CCE_WT_RA_1hr_2       CCE   WT        RA  1hr         2
CCE_WT_RA_1hr_3       CCE   WT        RA  1hr         3
CCE_WT_RA_8hr_1       CCE   WT        RA  8hr         1
CCE_WT_RA_8hr_2       CCE   WT        RA  8hr         2
CCE_WT_RA_8hr_3       CCE   WT        RA  8hr         3
CCE_WT_RA_24hr_1      CCE   WT        RA 24hr         1
CCE_WT_RA_24hr_2      CCE   WT        RA 24hr         2
CCE_WT_RA_24hr_3      CCE   WT        RA 24hr         3
RARE_KO_Ctrl_1       RARE   KO      Ctrl Ctrl         1
RARE_KO_Ctrl_2       RARE   KO      Ctrl Ctrl         2
RARE_KO_Ctrl_3       RARE   KO      Ctrl Ctrl         3
RARE_KO_RA_1hr_1     RARE   KO        RA  1hr         1
RARE_KO_RA_1hr_2     RARE   KO        RA  1hr         2
RARE_KO_RA_1hr_3     RARE   KO        RA  1hr         3
RARE_KO_RA_8hr_1     RARE   KO        RA  8hr         1
RARE_KO_RA_8hr_2     RARE   KO        RA  8hr         2
RARE_KO_RA_8hr_3     RARE   KO        RA  8hr         3
RARE_KO_RA_24hr_1    RARE   KO        RA 24hr         1
RARE_KO_RA_24hr_2    RARE   KO        RA 24hr         2
RARE_KO_RA_24hr_3    RARE   KO        RA 24hr         3
```
We will now create a subset of phenod with only the phenotype for our 6 selected samples: 

```{r}
> M6pheno <- phenod[c(1,2,3,10,11,12), 4]
> M6pheno
[1] "Ctrl" "Ctrl" "Ctrl" "24hr" "24hr" "24hr"
```

We will use the limma  (Linear Models for Microarray Data) package from Bioconductor to perform the analysis of diferential gene expression:

```{r}
> #source("http://bioconductor.org/biocLite.R")
> #biocLite("limma")
> library(limma) 
```

We will create a two groups comparison design matrix:

```{r}
> design <- model.matrix(~ factor(M6pheno, levels=c("Ctrl", "24hr")))
> colnames(design) <- c("Ctrl", "24hrvsCtrl")
```

And fit a linear model for each gene in the expression set M6data given the design matrix:

```{r}
> fit <- lmFit(M6data, design)
```

Now we are going to calculate the **differential expression** by **empirical Bayes shrinkage** of the standard errors towards a common value, by computing the moderated t-statistics, moderated F-statistic, and log-odds:

```{r}
> fit <- eBayes(fit)
```

Now, we will create (and save) a table with the calculated statistics for the selected group of samples sorted by corrected p-value:

```{r}
> tab <- topTable(fit, coef = 2, adjust = "fdr", n = dim(M6data)[1])  
> write.table(tab ,file="Allprobes-pq.xls",row.names=T, quote=FALSE, sep="\t")
```
Next, we are going to filter genes that have adjusted **p-values** less that 0.001, to create a smaller list of significant genes:

```{r}
> tab.sig <- tab[tab$adj.P.Val < 0.001,]
> dim(tab.sig)
[1] 367   7
```
And save those highly significant DEGs in a Excel like format:

```{r}
> write.table(tab.sig,file="DEG-FDR0.001.xls", quote=FALSE, row.names=T, sep="\t")
```

We are now going to sort the table of significant genes by **FDR adjusted p-value** and **fold change**

```{r}
attach(tab.sig)
sortedDEG <- tab.sig[order(adj.P.Val, logFC),]
```
After sorting, we can select and save a list of the top 20 significant DGEs.

```{r}
top20DEGs<-head(sortedDEG,20)
write.table(top20DEGs,file="top20DEGs.xls", quote=FALSE, row.names=T, sep="\t")
```


##12. Visualizing the expression of DEGs in a heatmap

Next we are going to create a **heatmap** of the expression of highly significant genes:

```{r}
> selected  <- tab.sig[,1] ## select the IDs of significant genes 
> M6dataSel <- M6data[rownames(M6data) %in% selected, ] ## subsetting
> ## Ploting a heatmap of the 367 Significant DEGs
> heatmap(M6dataSel, labRow=c(""), col=topo.colors(16), cexCol=0.6,  xlab = "Samples", ylab = "Features", main = "DEGs")
```
![heatmap](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/heatmap.png)

##13. Annotation 

The final step in this analysis will be the annotation of the DEGs. First we are going to download and install the libraries we need for the annoation:

```{r}
> ##biocLite("lumiMouseAll.db")
> ##library(lumiMouseAll.db)
> ##biocLite("AnnotationDbi")
> ##library(AnnotationDbi)

> install.packages("R2HTML")
> library(R2HTML)
```
To list the objects available in this annotation package we can use:

```{r}
> ls("package:lumiMouseAll.db")
```
Output:

```
[1] "lumiMouseAll"              "lumiMouseAll.db"          
 [3] "lumiMouseAll_dbconn"       "lumiMouseAll_dbfile"      
 [5] "lumiMouseAll_dbInfo"       "lumiMouseAll_dbschema"    
 [7] "lumiMouseAllACCNUM"        "lumiMouseAllALIAS2PROBE"  
 [9] "lumiMouseAllCHR"           "lumiMouseAllCHRLENGTHS"   
[11] "lumiMouseAllCHRLOC"        "lumiMouseAllCHRLOCEND"    
[13] "lumiMouseAllENSEMBL"       "lumiMouseAllENSEMBL2PROBE"
[15] "lumiMouseAllENTREZID"      "lumiMouseAllENZYME"       
[17] "lumiMouseAllENZYME2PROBE"  "lumiMouseAllGENENAME"     
[19] "lumiMouseAllGO"            "lumiMouseAllGO2ALLPROBES" 
[21] "lumiMouseAllGO2PROBE"      "lumiMouseAllMAPCOUNTS"    
[23] "lumiMouseAllMGI"           "lumiMouseAllMGI2PROBE"    
[25] "lumiMouseAllORGANISM"      "lumiMouseAllORGPKG"       
[27] "lumiMouseAllPATH"          "lumiMouseAllPATH2PROBE"   
[29] "lumiMouseAllPFAM"          "lumiMouseAllPMID"         
[31] "lumiMouseAllPMID2PROBE"    "lumiMouseAllPROSITE"      
[33] "lumiMouseAllREFSEQ"        "lumiMouseAllSYMBOL"       
[35] "lumiMouseAllUNIGENE"       "lumiMouseAllUNIPROT" 
```

In general *ProbeIDs* are used to map microarray measurements and probes. However, Illumina’s ProbeIDs can change with different versions, and even between different batches of Illumina microarrays. To solve this problem, **lumi** designed a nucleotide universal identifier **(nuID)**, which encodes the **50mer oligonucleotide sequence** and contains error checking and self-identification code. We will now annotate **nuIDs** with the corresponding **gene** information:

```{r}
> ID <- row.names(tab.sig)
> head(ID)
```

Output:

```
[1] "ZUsnIvMeX7f0IdPiAo" "ELiaujnucbuA3t.lKU" "iohWMQGvTB5.TFycp4"
[4] "o1Px7ENOBJ3iNei47o" "03lwTqhYBkxUVSn3lI" "Nd66OX794B5zsYVNe4"
```

```{r}
> Symbol <- getSYMBOL(ID, 'lumiMouseAll')
> Name <- as.character(lookUp(ID, 'lumiMouseAll', "GENENAME"))
> Ensembl <- as.character(lookUp(ID, 'lumiMouseAll', "ENSEMBL"))
```

For each Ensembl ID (if we have it), we will now create a hyperlink that goes to the Ensembl genome browser:

```{r}
> Ensembl <- ifelse(Ensembl=="NA", NA,
                  paste("<a href='http://useast.ensembl.org/Mus_musculus/Gene/Summary?g=",
                        Ensembl, "'>", Ensembl, "</a>", sep=""))
                        
```

And make a temporary data frame with all those identifiers:

```{r}
> tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, logFC=tab.sig$logFC, pValue=tab.sig$P.Value, FDR=tab.sig$adj.P.Val, Ensembl=Ensembl)

> row.names(tmp)<-NULL
```
Finally, we are now going to create one **html** file and one **text** file with the list of DEGs. The html file contains clickable links to the Ensembl Genome Browser.

```{r}
> HTML(tmp, "out.html", append=F)
> write.table(tmp,file="target.txt",quote=F,row.names=F, sep="\t")
```

Browse the created out.html file on your working directory and try the links.

## Week 6 Homework: :house: 

In this tutorial, we detected the DEGs between WT Control (3 samples) and KO Treatment at 24 hours (3 samples). Using similar analysis workflow, develop an R analysis pipleine (script) to detect the DEGs between WT Control and KO Treatment at 8 hour. Compare the list of top 20 DEGs between KO Treatment at 24 hours and KO Treatment at 8 hour.

Submit your homework via e-mail.











