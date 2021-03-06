
# Developing R programming skills for Bioinformatics  

**Center for Research Informatics, University of Chicago**

Saturdays 05/04/2019; 9:00 AM - 12:00 PM

**Instructor:** Jorge Andrade, Ph.D.


## Learning Objectives

- Develop basic R programming skills and learn how to use R for baisc statistics
- Develop a working knowledge of how to use R and Bioconductor packages for the analysis of genomics data 

## 1. Introduction
R is a free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms, Windows, and MacOS. To download R, please choose your preferred CRAN mirror.

### Where to find R:

* [CRAN](http://cran.r-project.org/) The Comprehensive R Archive Network
	* [Manuals](http://cran.r-project.org/manuals.html)  
	* [FAQs](http://cran.r-project.org/faqs.html)
	* [Contributed Guides](http://cran.r-project.org/other-docs.html)
* [R Home](http://www.r-project.org/) The R Project for Statistical Computing
	* [R Wiki](http://wiki.r-project.org/) 
	* [R Journal] (http://journal.r-project.org/)
	* [Mailing Lists] (http://www.r-project.org/mail.html)
	* [Bioconductor] (http://www.bioconductor.org/)

#### Advantages of using R

  * R is a computational environment - somewhere between a program and a programming language
  * No buttons, no wizards, only a command line interface
  * Is a professional statistics toolset, likely the only analyses tool you will ever need 
  * Is also a programming language 
  * Can handle large datasets (like NGS data)
  * Very powerful graphics
  * State-of-the-art statistics and tool of choice for bioinformatics 
  * Lots of available bioinformatics analysis tools/packages (i.e. Bioconductor)
  * Free, and open source!


#### Disadvantages of using spreadsheets (like Excel)

* Hard to handle large dataset (>1000 data points)
* Inflexible, few analyses available
* Hard to repeat analyses systematically with new data
* Not'reproducible research' friendly -  it is hard to document process and actions


#### Getting started:

* Download and install the latest version of R: https://www.r-project.org 
* Download and install the Open Source version of [Rstudio](https://www.rstudio.com). 

  * Start Rstudio, in a File/New File/R Script type:
  * Review and get familiar with Rstudio user interface.

  
  ```{r}
  > demo(graphics)
  ```

  * Hit enter a few times


#### Getting help:

Most R functions have online documentation.

```{r}
> help(lm) 
> ?plot
> help.search("aov")
> apropos("plot")
```
  
## 2. Basic R objects:
#### Vectors
  The basic Object in R is **a vector**. In statistics, we are almost always operating with several "data points" that can be stored/organized in a vector.
  
 A vector is a collection of numbers and/or strings:
 
- ("jorge", "wenjun", "ron") 
- (10, 5.2, 1,  7,  2, 21) 
- (3) 

The last example is a vector of length 1

In R, we make a vector by the **c()** command  (c stands for concatenate)

```{r}
> c(1,5,10, 7, 2, 1)

> c("jorge", "wenjun", "ron") 
```
When creating vectors of **strings** or **characters**, we have to use  **"** or **'**.
  
If we are making vectors of size 1, we can skip **c()**:

```{r}
> 3
```
When you start a session in R, it is a good idea to first clean all R objects, see below:


```{r}
> ls()     		# List the contents of the workspace. 
> rm(list=ls())         # remove all objects in your workspace. 
> ls()                  # character(0) means "nothing to see here"
```

#### Task1:
-------------------------
* Make the following vector: 45,5,12,10

* What happens with the following commands?
  - c(1:100)   
  - c(50:2) 
<details> 
  <summary>Generate a vector containing the 50 first odd numbers </summary>
 
```{r}
> (1:50)*2 - 1

> seq(1, 99, 2)

> seq(1, by = 2, len = 50)
```
</details>

------------------------



A vector is a **data structure**, and the most fundamental in R. Almost everything in R is some kind of vector, although sometimes in several dimensions - vectors within vectors (matrix). 

#### A matrix:

We can create a memory representation of a matrix in R with the **matrix** function (the data elements must be of the same basic type). 

See below:

```{r}
> A <- matrix(c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3, byrow = TRUE)
> A                     
     [,1] [,2] [,3] 
[1,]    2    4    3 
[2,]    1    5    7

### Try replacing byrow = FALSE and removing byrow = TRUE all together
```

#### Use a Reference sheet:

You may get overwhelmed by different command names fast, use the reference sheet available [here](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/ShortRrefcard.pdf) as a guide while you learn.

#### Assignment to memory

The **c()** command is almost useless by itself - we want to keep/save/use the contents of a vector for other analyses, we can do that using the assignment concept:

```{r}
> 4+5           # add 4 and 5
> a <-  4       # store 4 as "a"
> b <-  5       # store 5 as "b"
> a             # just checking
> b
> a+b           # add a+b (4+5)
```
In R the symbol **<-** and **=** have the effect of assigning **value** to a **variable**. Many R traditionalists prefer **<-** for clarity (to distinguish **'is assigned to'** from **'is equal to'**). 

#### Expanding assignment to a whole vector:

```{r}
my_vector <- c(1,5,10, 7, 2)
```
The commands return no value to the screen, a variable of type **vector** called "my_vector" is created. Variable names are totally arbitrary.

The anatomy of the vector:

Name   | my_vector |        |        |       |       | 
------ | --------- | ------ | ------ | ----- | ----- | 
Values |    1      |   5    |   10   |   7   |   2   | 
Index  |   [1]	   |  [2]   |	[3]  |  [4]  |	[5]  | 


We can access part of the vector using the corresponding index:

```{r}
 > my_vector[5]   # will give you the 5th item in the vector
```

Explore the effect of the following commands:
  
```{r}
> my_vector<- c(1,5,10,7,2) 
> my_vector [c(1,3,5)] 
> my_vector[1:4] 
> my_vector[4:1]
```

In a matrix, an element at the *m* row and *n* column, can be accessed by the expression A[m,n]:

```{r}
> A[2, 3]      # element at 2nd row, 3rd column 
[1] 7
```
The entire *m* row of matrix A can be extracted as A[m, ]:

```{r}
> A[2, ]       # the 2nd row 
[1] 1 5 7
```

#### Task2


Given the following vector:
```{r}
> my_vector<- c(1,5,10, 7, 2)   
```
<details> 
  <summary> T2: Using the reference sheet, figure out at least three ways of making R print your vector in the reverse direction </summary>
   A2: 

```{r}
> my_vector[5:1]
> my_vector[c(5,4,3,2,1)]
> c<- c(my_vector[5],my_vector[4],my_vector[3],my_vector[2], my_vector[1])
> rev(my_vector)
```
</details>

----------------

#### Naming rules and the danger of over-writing 

A good practice, in any programming language, is to assign meaningful names to variables/vectors. In R, never start a vector name with a number:

```{r}
> a<- c(1,5,4,2)    #OK
> 1a<- c(1,5,4,2)		# NOT OK  Error: syntax error
> a1<- c(1,5,4,2)   # OK
```

#### Over-writing

As with any variable assignment, you must be aware of overwriting. Try: 

```{r}
> my_vector<- c(1,5,10, 7, 2)
> my_vector<- c(10,5,2, 3, 1)
```
What does my_vector contain now?

----------------
### Vector operations

#### Analyzing vectors

There are many prebuild R functions to work with vectors. Most have logical/meaninfull names. For instance, **length(my_vector)** gives the number of items in the vector (= 5).

#### Task3: 

<details> 
  <summary> T3: Create a vector called big_vector with values 1 to 10000, then find: Length of the vector; Sum of all items in vector;
  Mean(average) of all items in the vector; Median value of the vector.
</summary>
   A3: 
	
```{r}
> big_vector<-(1:10000); length(big_vector)
> sum(big_vector)
> mean(big_vector)
> median(big_vector)
```
</details>
	
	
#### Task4: 

<details> 
  <summary> T4: Explore at the help for the command sample(), sort(), plot(), and hist(). Then use them on big_vector
</summary>
   A4: 

```{r}
> x <-sample(big_vector,100)
> x[1:20]
> plot(x)
```

```{r}
> sort(x, decreasing = FALSE)[1:20]
> hist(x)
```
</details>

------------------------------

#### Adding and multiplying a number to a vector

Sometimes we want to add a constant (like a number), to each element in the vector, test the following code: 

```{r}
> big_vector2 <- big_vector+10 
> min(big_vector) 
> max(big_vector) 
> min(big_vector2) 
> max(big_vector2) 
```
#### Adding vectors

We can also add one vector to another vector, create the following vectors: 
  
```{r}
> A <- c(10, 20, 30, 50) 
> B <- c(1,4,2,3) 
> C <- c(2.5, 3.5) 
```
Test what happens and explain the outcome: 

```{r}
> A+B 
> A+C
```

#### Adding vectors
```{r}
> A<-c(10 ,  20, 30, 50)
> B<-c(1  ,   4,  2,  3)
> C<-c(2.5, 3.5        )
> A+B
> A+C
```
A+B is easy to understand : A[1]+B[1] , etc. 
A+C is trickier - the C vector is just of length 2. It is re-used. This is what is happening also with the earlier A+10 example. The 10 is used many times.


### Plotting vectors

Let's make up some semi-random data:

```{r}
> dat <- rnorm(100)   ### Draw 100 random normal distributed data points 
```
Test the following:

```{r}
> plot(dat)
> plot(dat,type='l')
> barplot(dat)
> hist(dat)
```

Now observe the effect o using the command: par( mfrow=c(2,2) )


```{r}
> par( mfrow=c(2,2) )
> dat<-rnorm (100) 
> plot(dat)
> plot(dat,type="l")
> barplot(dat)
> hist(dat)
> dev.off()
```

- Why are your three first plots different from mine? 
- Why is your last plot more similar to mine?
- What is the effect of the command: *> dev.off()*

#### Graph options

You can give extra 'options' to graphical commands, for example:

```{r}
> plot(dat, col='blue', type='l')
```

In plot try the effect of the following options: - note that you can use several at once  

- type='b'  
- col='hotpink'  
- main='plot'  
- type='h'  
- type='S'  

These options are 'arguments' to the *plot()* function


```{r}
> plot(dat, col='blue', type='l')
```
```{r}
> plot(dat, col='hotpink', type='h', main='Plot')
```
We will learn more on how to use R to create data vizualizations on the next lesson.

#### More about functions

In most cases, a function needs some **input**, like in *plot(dat)* 'dat' here is an *unnamed argument,* and this works because **plot()** assumes we mean: x values = dat.

We could also say: **plot(x=dat)**  - a *named argument.* If you have many arguments, most of them are named - such as in the following example:

- plot (some_vector, col="blue", type="s"). 

The help pages will tell you what type of arguments you can use for a particular command.

#### The danger of unnamed arguments.

The order of unnamed arguments, could make a big difference. 

Try the followoing commands. Observe the differences 

```{r}
> par(mfrow=c(1,3))
> a<-rnorm(100)
> b<-rnorm(100)*2
> plot(a,b, main="plot(a,b)")
> plot(b,a, main="plot(b,a)")
> plot(x=b, y=a, main="plot(x=b,y=a)")
```

```{r}
> dev.off()
```

#### Some generic R arguments to plots - the par() function

  * The **par()** function (we have usedit before), is used to set general plot properties. It has hundreds of possible arguments - see: **?par**
  * Two very handy **par()** arguments are **mfrow()** and **mfcol()** - these will allow many plots in one page
  * You give these functions a vector of length 2 - this gives the number of cells in the page, see examples bellow:


#### A:

To print in 4 windows in the same page:
```{r}
> par( mfrow=c(2,2) )
> plot(a,b); plot(b,a); plot(x=b, y=a)
```
```{r}
> dev.off()
```

#### B:

To print the three plots in one row:

```{r}
> par( mfrow=c(1,3) )
> plot(a,b); plot(b,a); plot(x=b, y=a)
```
```{r}
> dev.off()
```

#### Overlaying plots

  * Sometimes we need to add additional information to a plot. This is often made by the **lines()** or **points()** commands:


```{r}
> plot(b, type="l", col="blue")
> lines(a, col="red")
```

- Why did we start with plotting b? 
- What would have happed if we use **points()** instead of **lines()**?


#### Sizing graphs

  * Change X scale: **xlim=c(start_value, end_value)**
  * Change Y scale: **ylim=c(start_value, end_value)**


Let's try the code bellow:

```{r}
> par(mfrow=c(1,2))
> plot(a, type="l", col="blue")
> plot(a, type="l", col="blue", ylim=c(-5,5))
```

#### Saving graphs

  * All systems can use the **device()** function - see: **?device**

``` {r}
> # Saving to a .pdf file 
> setwd('/Users/jandrade/Desktop')
> pdf('plot.pdf') 
> plot(a,b) 
> dev.off() 
> getwd()
```

``` {r}
> # Saving to a .jpg file 
> setwd('/Users/jandrade/Desktop')
> jpeg('rplot.jpg') 
> plot(a,b) 
> dev.off()
> getwd()
```

## 2. Some Basic Statistics with R

### Mean, median, histogram and boxplot


#### Summary statistics

  - **hist()** (= Histogram) We used the hist() function before, as showed, hist() is a graphical way of summarizing distributions - it creates a number of "bins" and calculates how many of the data points fall into each bin. 
  - We can also summarize by the center point in the data:
  - **mean():**

![mean](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/mean.png)

  - **median():**
Sort the data, pick the number in the center. If the number of data points is even, take the mean of the two center points.

#### Task:
 
  * Create a vector of 100 elements
  * Add a few extra points
  * Test: **mean()** and **median()** 
  * Are they the same? 
  * Explain the differences by plotting a histogram 
  * What is the advantage/disadvantage of each measure?


```{r}
> dat<-rnorm(100)
> dat2<-c(dat, 102, 105, 130 )
> median(dat2)
> mean(dat2)
> par(mfrow=c(1,2))
> hist(dat)
> hist(dat2)
```

```{r}
> dev.off()
```
**Means are sensitive to outliers! Very common situation in genomics.** 

### Boxplots are a better way to visualize outliers

```{r}
> par(mfrow=c(1,2))
> boxplot(dat, main='dat') 
> boxplot(dat2, main='dat2') 
```

```{r}
> dev.off()
```

### Setting margins

As we can see, it is important to be able to set the margins for plots, we can do that with the **mar** function.
**mar** – A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. 

The default is: c(5.1, 4.1, 4.1, 2.1).

```{r}
> par(mar=c(4, 4, 4, .1))
> dat <- rnorm(100, mean = 50, sd = 2)
> dat2<-c(dat, 105, 110, 200 )
> par( mfrow=c(1,2) )
> boxplot(dat, main='dat'); boxplot(dat2, main='dat2')
```

Note the parameters of **rnorm(n, mean = , sd = )**


```{r}
> dev.off()
```

#### Percentiles

  * An extension of the median concept
  * Best explained by example: 
  * The 20th percentile is the value (or score) below which 20 percent of the observations may be found. 
  * The median is the same as the 50th percentile 
  * The first quartile is the 25th percentile, the third is the 75th

#### Try: summary(dat) and summary(dat2)

```{r}
> summary(dat)
> summary(dat2)
```

The command **ecdf()** (Empirical Cumulative Distribution Function) calculates "all" percentiles in your data - and also understands the **plot()** function. 

#### Try: 

```{r}
> dev.off()
> plot (ecdf(dat2))
> plot (ecdf(dat1))
```

The figure shows what fraction of the data has been covered at point X. 

#### Boxplots

An "easier" way to represent the distribution of data is a Boxplot. A Boxplot show you both the center point and the "spread" of the data:

![BOXPLOT](https://upload.wikimedia.org/wikipedia/commons/1/1a/Boxplot_vs_PDF.svg)

Boxplot (with an interquartile range) and a probability density function of a Normal N(0,σ2) Population


 	 - The first quartile (Q1) is defined as the middle number between the smallest number and the median of the data set (splits off the lowest 25% of data from the highest 75%).
	 - The second quartile (Q2) is the median of the data (cuts data set in half)
	 - The third quartile (Q3) is the middle value between the median and the highest value of the data set (splits off the highest 25% of data from the lowest 75%)
	 - The interquartile range (IQR), also called the midspread or middle 50%, is a measure of statistical dispersion, being equal to the difference between 75th and 25th percentiles, or between upper and lower quartiles. In other words, the IQR = 3rd quartile - 1st quartile.
	-  Outliers here are defined as observations that fall below Q1 − 1.5 IQR, or above Q3 + 1.5 IQR
 
 
 * R can calculate all these values and plot them with a simple command: **boxplot()**

```{r}
> boxplot(dat) 
> rug(dat, side=2)
```

Outliers are common observations in Genomics data.
  

![outliers](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/outliers.png)

______________________________

#### Variance, standard deviation and  data spread

Try the follwoing code:

```{r}
> par(mfrow=c(1,3))
> hist(rnorm(100,sd=1), xlim=c(-100, 100), main="histogram of rnorm(100,sd=1)")
> hist(rnorm(100,sd=10),xlim=c(-100, 100), main="histogram of rnorm(100,sd=10)")
> hist(rnorm(100,sd=100),xlim=c(-100, 100), main="histogram of rnorm(100,sd=100)")
```


![distributions](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/distributions.png)


What is the difference between these distributions?

  * Same mean and median, but different spread over the x-axis.
  * This can be measured by the variance of the data. 
  * The variance (σ2) is a measure of how far each value in the data set is from the mean.

![variance](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/variance.png)
  
  * It is basically the difference between each point and the mean, squared.  

  * The Standard Deviation is simply the variance squared:

  ![Standard](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/StandardDeviation.png)
  
#### Task:

Observe the ooutput of the two code blocks bellow:

```{r}
> x<-rnorm(100, sd=1)
> y<-rnorm(100, sd=10)
> z<-rnorm(100, sd=100)
> par( mfrow=c(1,3) )
> hist(x, xlim=c(-100, 100))
> hist(y, xlim=c(-100, 100))
> hist(z, xlim=c(-100, 100))
```

```{r}
> x<-rnorm(100, sd=1, xlim=c(-100, 100)) 
> y<-rnorm(100, sd=10, xlim=c(-100, 100))
> z<-rnorm(100, sd=100, xlim=c(-100, 100))
> par( mfrow=c(1,3) )
> hist(x); hist(y); hist(z)
```

#### Why is variance and standard deviation important?

  * Variance inform us about the quality of the measurements. It measures how far a set of (random) numbers are spread out from their mean
  * The higher the variance, the harder it is to say with certainty that two measurements are different
  * Standard Deviation (SD): is a measure that quantifies the amount of variation or dispersion of a set of data values

  * R functions for variance and standard deviation are: **var()** and **sd()**. 
   
Let us try them with some random data:

```{r}
> smallset<-rnorm(100) 
> largeset<-rnorm(10000) 
```
  * What is the variance and standard deviation for these?
  * Is the standard deviation really the square root of the variance (what is the R function for square root?)    


```{r}
> var(smallset)
> var(largeset)
> sd(largeset)
> sd(smallset)
> sqrt(var(smallset)) ##Calculating SD manually
```
Why do we get about the same variance?  

```{r}
> ?rnorm
```

## 4. Basic R programming


### Control Flow

The flow control in R is defined by the **if** reserved word, The syntax of if statement is:

```
if (test_expression) {
   statement
}
```
If the *test_expression* is TRUE, the **statement** gets executed. But if it's FALSE, nothing happens.

Here, test_expression can be a **logical** or **numeric vector**, but only the first element is taken into consideration.

In the case of numeric vector, zero is taken as FALSE, rest as TRUE.

#### Example: if statement


```{r}
x <- 5
if(x > 0){
   print("Positive number")
}
```
```
[1] "Positive number"
```

#### if...else statement

The syntax of **if...else** statement is:

```
if (test_expression) {
   statement1
} else {
   statement2
}
```

The else part is optional and is only evaluated if test_expression is FALSE.



#### Example of if...else statement

```{r}
x <- -5
if(x > 0){
   print("Non-negative number")
} else {
   print("Negative number")
}
```
Output
```
[1] "Negative number"
```

The above conditional can also be written in a single line as follows:

```{r}
> if(x > 0) print("Non-negative number") else print("Negative number")
```

#### Nested if...else statement

We can nest as many if...else statement as we want as follows:

#### The syntax of nested if...else statement is:

```
if (test_expression1) {
   statement1
} else if (test_expression2) {
   statement2
} else if (test_expression3) {
   statement3
} else
   statement4
```


Only **one statement** will get executed depending upon the *test_expressions*.

####  Example of nested if...else

```{r}
x <- 0
if (x < 0) {
   print("Negative number")
} else if (x > 0) {
   print("Positive number")
} else
   print("Zero")
   
```
Output

```
[1] "Zero" 
```

### For, while, and repeat loops

A loop is a way to repeat a sequence of instructions under certain conditions. They allow you to automate parts of your code that are in need of repetition


#### For loop

```{r}
for (i in 1:10) {
  print(i)
}
```
####  Example of a for loop

```{r}
for (year in c(2012,2013,2014,2015,2016,2017,2018))
{
  print(paste("The year is", year))
}
```
Output:

```
[1] "The year is 2012"
[1] "The year is 2013"
[1] "The year is 2014"
[1] "The year is 2015"
[1] "The year is 2016"
[1] "The year is 2017"
[1] "The year is 2018"
```
Another example to count the number of **even** numbers in a vector.

```{r}
x <- c(2,5,3,9,8,11,6)
count <- 0
for (val in x) {
    if(val %% 2 == 0)  count = count+1
}
print(count)
```

Output:

```
[1] 3
```

#### While loop

Syntax of while loop:

```
while (test_expression)
{
   statement
}
```
Here, test_expression is evaluated and the body of the loop is entered **if** the result is TRUE.

The statements inside the loop are executed and the flow returns to evaluate the *test_expression* again.

This is repeated each time until *test_expression* evaluates to FALSE, in which case, the loop exits. This means than the statement in a **while** loop, can get executed 0 or more times.

#### Example of a while loop

```{r}
i <- 1

while (i < 6) {
   print(i)
   i = i+1
}
```

Output:

```
[1] 1
[1] 2
[1] 3
[1] 4
[1] 5
```


#### Repeat loop

A repeat loop is used to iterate over a block of code multiple numbers of times. There is no condition check in repeat loop **to exit** the loop.

We must ourselves **put a condition explicitly inside the body** of the loop and use the break statement to exit the loop. Failing to do so will result in an infinite loop.

#### Syntax of repeat loop

```
repeat {
   statement
}
```

In the statement block, we must use the **break** statement to exit the loop.

#### Example of a repeat loop

```{r}
x <- 1

repeat {
   print(x)
   x = x+1
   if (x == 6){
       break
   }
}
```

Output

```
[1] 1
[1] 2
[1] 3
[1] 4
[1] 5
```


**Note:** If you make the like with 'break' a comment (as in # break), you would have created an infinite loop.



#### User defined functions

In programming, we use functions to incorporate sets of instructions that we want to use repeatedly, or that because of their complexity are better self-contained in a sub program and called when needed. A function is a piece of code written to carry out a specified task.

In R you define a function with the construct:
```
function ( arglist )  {body}
```
Where the code in between the curly braces is the body of the function.

Once the definition of the function is done, somewhere else in the code, we call the function (we use it). The following code defines a function that computes the square of the argument and then calls it after assigning a value for its argument.

```{r}
# define a simple function

myFirstFun<-function(n)
{
  n*n  #  compute the square of integer n
}

# define a value
> k<-10

# call the function with that value
> m<-myFirstFun(k)
```

We can define a function that receive several 'arguments':

```{r}
findSumSquare <- function(a,b) 
{
  return(a^2+b^2)
}
# Using the findSumSquare function:
> a=3
> b=5
> findSumSquare(a,b)
> findSumSquare(1,2)
```
Bellow a simple example of how to use a user defined function in R to compute the first *n* numbers of the [fibonacci](https://math.temple.edu/~reich/Fib/fibo.html) series:

```{r}
fibonacci <- function(n)
{
len <- n
fibvals <- numeric(len)
fibvals[1] <- 1
fibvals[2] <- 1
for (i in 3:len) { 
  fibvals[i] <- fibvals[i-1]+fibvals[i-2]
                 } 
return (fibvals) 
}

## Calling our function

> fibonacci(5)
[1] 1 1 2 3 5
```


## 5. Basic data exploration

### How to load data files

R can process input data sets can be in various formats (.XLS, .TXT, .CSV, JSON ). In R, it is easy to load data from any source, due to its simple syntax and availability of predefined libraries. Here, we will take examples of reading a CSV file and a tab separated file.



```{r}
> mydata <- read.csv(file="https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/binary.csv", header=T) 
> head(mydata, 4)
> summary(mydata[,2:3])
```

* Reading a comma delimitted text file (csv) from your compueter: (files are available in the 'data' folder on the GitHub)

```{r}
# To read a file from your working directory use:
# > data1<-read.csv("NeuralStemCellData.csv", row.names=1, header=T)
# To read a file from a web link:
> data1<-read.csv("https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/NeuralStemCellData.csv", row.names=1, header=T)
> head(data1, n=7)
```

* Loading data from a tab delimited text file 
  
```{r}
# To read a file from your working directory use:
# > data2<-read.delim( "NeuralStemCellData.tab", row.names=1, header=T)
# To read a file from a web link:
data2<- read.delim("https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/NeuralStemCellData.tab",row.names=1, header=T)
> head(data2, n=7)
```

* One way to read several formats: **read.table()**

    + comma-separated values data (csv)
    
        - `read.table("filePath/fileName", sep=",", ...)` 
        
    + tab-delimited data
    
        - `read.table("filePath/fileName", sep="\t", ...)` 
        
    + space-delimited data
    
        - `read.table("filePath/fileName", sep=" ", ...)`
		

### Loading data with missing values:

```{r}
> mydata.missing <- read.csv( "https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/data/NeuralStemCellDataMissing.csv", row.names=1, na.strings = c(" ", "NA", "NaN"))
> head(mydata.missing , n=15)
```

## Knowning the data 

**head()** and **tail()** - returns the first or last lines of a vector, matrix, table, data frame or function

```{r}
> head(mydata)
> head(mydata, n = 10)
> tail(mydata)
> tail(mydata, n = 10)
```
### Names of the variable in the dataset
```{r}
> names(mydata)
```

### Number of rows and columns in the dataset
```{r}
> dim(mydata)
```

We will learn more on data exploration on the following lesson.

## 6. Basic data transformation and model fitting

#### Add a new column: id
```{r}
> mydata$id=1:400
> head(mydata, 5)
```


#### Selecting a column: gre, rank
```{r}
> mydata$gre
> mydata$rank
```


#### Subset the dataset
 
 1. By index

```{r}
> mydata[1:10, 2:3]
```  
Output

```
  gre  gpa
1  380 3.61
2  660 3.67
3  800 4.00
4  640 3.19
5  520 2.93
6  760 3.00
7  560 2.98
8  400 3.08
9  540 3.39
10 700 3.92
```

2. By conditions

```{r}
> mydata[mydata$rank <= 1,]  
> mydata[mydata$id>=10 & mydata$id<=20,]   

``` 

Output

```
admit gre  gpa rank id
10     0 700 3.92    2 10
11     0 800 4.00    4 11
12     0 440 3.22    1 12
13     1 760 4.00    1 13
14     0 700 3.08    2 14
15     1 700 4.00    1 15
16     0 480 3.44    3 16
17     0 780 3.87    4 17
18     0 360 2.56    3 18
19     0 800 3.75    2 19
20     1 540 3.81    1 20
```


#### Boxplot for GRE and admit
```{r}
> boxplot(gre~admit, data=mydata, xlab="Admit", 
        ylab="GRE", main="Boxplot GRE and admit")
```

### Data management - sorting 
* Sort or order a data vector into ascending or descending order

```{r}
> head(sort(rownames(mydata), decreasing=F))
> head(order(rownames(mydata), decreasing=T))
```    

#### Data management - sorting 

```{r}
# sort by gpa
> head(sortdata <- mydata[order(mydata$gpa),])

# sort by gre and gpa
> head(sortdata <- mydata[order(mydata$gre, mydata$gpa),])

```


## 7. Writing datasets to a text file

```{r}
> write.table(mydata, file="test.txt", sep="\t", row.names=FALSE, quote=FALSE)

> list.files(path=getwd(), pattern="test.txt", full.names=T) ## Check file exixt; return name

```

## 8.  Statistical hypothesis testing

The **t.test()** function produces a variety of t-tests. Unlike most statistical packages, the default assumes unequal variance and applies the [Welsh df modification](https://en.wikipedia.org/wiki/Welch%27s_t-test) 

#### T-test for gre scores for admit=0 vs admit=1

```{r}
> t.test(gre~admit,data=mydata)
```

#### Simple regression between gre and gpa

**lm()** is used to fit linear models. It can be used to carry out regression, single stratum analysis of variance and analysis of covariance

```{r}
> fit1 <- lm(gre~gpa, data=mydata) 
> plot(gre~gpa, data=mydata)
> abline(fit1,col="red")
```

#### Logistic regression

Generalized linear models are fit using the **glm()** function:

```{r}
> fit <- glm(admit~gpa+gre+factor(rank), data=mydata, family=binomial)
> print(summary(fit)$coef, digits=2)
```

## Week 5 Challange (5 Extra points)


- Copy all .fasta available at: **/gpfs/data/mscbmi/hw4/** into your R working directory:

- The follwoing R code uses the [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) package to perform pairwise alignment of all the sequences in the input file. The R code will generate an output file named **alignment.txt** containing the sequence alaignmets and alignment scores.

```{r}

library(seqinr)
library(Biostrings)

rm(list=ls()) 

seqsim <- function(seq)
{
  sink("alignment.txt") ## direct R output to a file
  for(i in 1:(length(seq)-1))
  {
    n1 <- toupper(c2s(seq[[i]]))
    n2 <- toupper(c2s(seq[[i+1]]))
    globalAlign<- pairwiseAlignment(n1, n2)
    print(globalAlign)
  }
  sink()
}

seq <- read.fasta(file = 'reference.fasta')
start.time <- Sys.time()
seqsim(seq)
end.time <- Sys.time()
time.taken <- end.time - start.time
write.table(time.taken, "timestamp.txt")
```

* A. Study and understand the code and all the R functions used
* B. Develop an analysis pipeline that uses CRI's GARDNER cluster to run the code
* C. Sort the **'score'** value for each aligned sequence and save the sorted scores in a file named: 'scores.txt' 
* D. What is the maximun and minimum alignment scores for reference.fasta?
* E. Record, report and explain the total runtime (execution time) for your script running with the following hardware configurations:
	- 1 nodes; 4 cpus; 4gb of RAM
	- 1 nodes; 8 cpus; 4gb of RAM
	- 1 nodes; 16 cpus; 8gb of RAM
	
* F. Are there opportunities for reducing the runtime using Parallelization?  If so, what could be the strategy? Algorithm?

Submit all developed scripts and results files via e-mail.


