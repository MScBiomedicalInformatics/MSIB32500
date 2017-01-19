# Introduction to Unix/Linux: 

**Center for Research Informatics, University of Chicago**

April - June 2017; Saturdays 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.

**Teaching Assistant:** Wenjun Kang, Ms.C.

## Learning Objectives

- This training will teach you all you need to know to start using Linux, be it on the CRI’s High Performance Infrastructure (HPC), or on your own computer (even if you are already running Windows or Mac OS).
- After this training you will: 
  - Feel comfortable using Linux
  - Know how software works on Linux  and how to use it
  - Use bash to execute commands in Linux and know your way around the file system
  - Have an overview of Linux tools that are extremely useful for bioinformatics

## Log on to CRI's High Performance Computing (HPC) system

In MacOS open the Terminal

:pushpin: If you are using a Windows machine (other than Windows 10 Anniversary Edition), you need install [MobaXter](http://mobaxterm.mobatek.net/)

```bash
ssh username@tarbell.cri.uchicago.edu
```
Enter your **_password_** when prompted. Type yes if you are prompted to accept a key.

## How to get help

- Use the manual ($ man) command; to exit the manual type 'q'
- Ask for help ($ your_command --help)
- Use the comand apropos ($ apropos text)

```bash
man ls
```
```bash
mkdir --help
```
```bash
apropos secure copy
```

## Navigation
a. List files in your home directory
```bash
ls            ### List the files in your current directory
ls -a         ### List the files in your current directory; do not ignore entries starting with .
```
b. Change current directory to the root of the file system and explore the directory structure
```bash
cd /          ### / represents the root of the file system
ls -l         ### List files in long format
pwd           ### Show the current directory
cd /tmp       ### Change current directory to /tmp
pwd
ls -l
```
c. Change to parent directory; one 'step' UP on the file tree
```bash
cd ..        
pwd
ls -l
```
d. Go back to your home directory
```bash
cd ~
pwd
```
## Operations with directories and files 

a. Create a new directory
```bash
mkdir mydir  ### Make a new directory called 'mydir'
mkdir newfolder1
ls -l 
```
b. Create a new file in a directory
```bash
cd mydir
nano file1.txt          ### Use nano editor to create a new file, write some text and use Control-O to save and Control-X to exit.
pwd
ls -l
```
c. Copy files and directories 
```bash
cp file1.txt copy.txt
ls -l
more copy.txt         ### Quick look to the contents of the file copy.txt
man more              ### Review the manual of the comand more; type q to exit
cp ~/mydir/*.txt ~/newfolder1/          ### copy all .txt files from mydir to newfolder1 
cd ~/newfolder1/
ls -l                 ### List the files you just copied
cd ..                 ### Change to parent directory
pwd
mkdir newfolder2
cd newfolder2
cp ~/newfolder1/file1.txt .     ### Copy file1.txt on current (.) directory
ls -l
cd ~
pwd
```
d. Remove a file or directory

```bash
cd newfolder1        
ls -l
rm file1.txt          ### Remove file1.txt
ls -l
cd                    ### Back to your home directory (this is equivalent to cd ~)
```

```bash
rm -r mydir           ### You will delete the whole mydir folder and all files on it, there is no way to recover those files
```
e. Rename a file or folder

```bash
cd newfolder1
mv copy.txt newname.txt         
ls -l
cd                              
mv newfolder1 newnamefolder                ### Rename newfolder1 with newnamefolder
```
f. Move files from one folder to another

```bash
mv newfolder2/file1.txt newnamefolder
```
g. Compress files

```bash
cd /group/bioinformatics/shared/MSIB32500/Lecture1/linux
ls -l SRR*                     ### List all files that start with SRR
gzip SRR001655.fastq           ### Compress a file, this command should create a compressed file named SRR001655.fastq.gz
ls -l                        
gunzip SRR001655.fastq.gz      ### Decompress a file
ls -l
```

## File transfer between computers

a. Windows user download and install [WinSCP](http://winscp.net/eng/index.php). MacOS users open the Terminal
The *scp* (secure copy command) allows you to copy/move files between compters on the command line
```bash
scp example1.txt username@tarbell.cri.uchicago.edu:.
```
:bulb: The [hypexr.org](http://www.hypexr.org/linux_scp_help.php) website has nice list of examples on how to use secure copy (also some [computer comics and cartoons](http://www.hypexr.org/comics.php))

b. You can use *wget* to get a file from the internet directly to your working directory in Linux
```bash
wget http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
```
## Input/Output redirect and pipe

a. Use the symbol '>' to redirect the output of a command to a file

```bash
nano text1.txt                       ### Create a new file and write some text on it
 
cat text1.txt                        ### Print file1.txt to screen
 
cat text1.txt > text2.txt            ### Print file1.txt to file2.txt
 
nano text2.txt                       ### Use nano to view file2.txt

```
b. Use the symbol '>>' to redirect the output of a command to a file using 'append'

```bash
cat file1.txt >> file2.txt           ### Append file1.txt to file2.txt
nano file2.txt                       ### View file2.txt
```
c. Use the symbol '<' to redirect the input
```bash
cat < file1.txt                      ### Print file1.txt to screen
```
d. Use the symbol '|' to pipe the output of one command as the input of another 

```bash
ls | wc -l              ### List the files and folders in the current directory and then count them
```
## Text extraction and manipulation

A standard Unix/Linux installation will have available several text editors (like: vi, vim, nano, emacs, and others) and text viewers (like: less, more, head and tail)

a. Use the symbol '>' to redirect the output of a command to a file

```bash
cd /group/bioinformatics/shared/MSIB32500/Lecture1/linux               
less SRR001655.fastq                 ### View a fastq file: use q to exit, space or f to the next page,
                                     ### b to the previous page, and / to search a word
head -20 SRR001655.fastq             ### Show the first 20 line of the file
tail -20 SRR001655.fastq             ### Show the last 20 lines of the file

```
b. Pattern Search with 'grep'. The grep command searches specified files or other input(stdin) for patterns matching a given expression(s).

```bash
$ cat list1.txt                      ### See the contents of file list1.txt
apples
bananas
plums
carrots
 
$ cat list2.txt                      ### See the contents of file list2.txt
Apple Sauce
wild rice
black beans
kidney beans
dry apples
 
$ grep apple list1.txt list2.txt     ### Search for "apple" in list1.txt and list2.txt
list1.txt:apples
list2.txt:dry apples
```
c. Finding a file
```bash
find . -name "*.txt"            ### Find all the files whose names ended with ".txt" under the current directory
find . -mtime 0                 ### Find all the files modified today under the current directory
find . -mtime +3                ### Find all the files modified more than 3 days ago under the current directory
find . -mmin -30                ### Find all the files modified in the last 30 minutes
find . -type d -name "ex*"      ### Find all the folders whose name starts with ex
```
d. Text replacement and text operation: cat, sed, tr, and rev
```bash
$ cat list2.txt                                ### See the contents of file list2.txt
Apple Sauce
wild rice
black beans
kidney beans
dry apples
 
 
$ cat list2.txt | sed 's/bean/button/g'       ### Replace bean with button
Apple Sauce
wild rice
black buttons
kidney buttons
dry apples
 
$ cat list2.txt | tr /a-z/ /A-Z/              ### Change lower case to upper case
APPLE SAUCE
WILD RICE
BLACK BEANS
KIDNEY BEANS
DRY APPLES
 
 
$ echo 'attgctcgat' | tr /atgc/ /tacg/ | rev  ### Find a complement DNA sequence and then reverse it
atcgagcaat
```

e. Table manipulation: sort, uniq, cut, awk, and paste

```bash
$ cat table.txt                                      ### show the contents on table.txt
 CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
 
  ### Note that the header line is not showing properly aligned with the contents. Now try:

$ cat table.txt | column -t
CHR  SNP         BP         A1  C_A  C_U   A2  CHISQ    P       OR
19   rs10401969  19268718   C   222  890   T   0.03462  0.8524  0.9857
1    rs10873883  76734548   G   934  3811  A   0.5325   0.4656  0.9691
1    rs11589256  214196749  C   271  1084  T   0.01928  0.8896  0.9902
15   rs10401369  19268718   C   232  890   T   0.03232  0.2524  0.1157
11   rs10873487  767334548  G   964  3811  A   0.5525   0.2356  0.2391
 
 
$ sort -k1n table.txt  > table_sorted1.txt         ### sort the table by the first column in numerical order
                                                   ### and write the results to a new file
 
$ cat cat table_sorted1.txt
CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
 
$ sort -k2 table.txt                               ### sort the table by the second column
 
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
 CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
 
#### Note that the header is now at the end of the table, we can fix it with 'awk':

$ awk 'NR==1;NR>1{print $0 | "sort -k2"}' table.txt | column -t
 CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
 
$ cut -f1,2,3,5 table.txt                ####Extract columns 1, 2, 3, and 5 from table.txt
CHR     SNP             BP              C_A
19      rs10401969      19268718        222
1       rs10873883      76734548        934
1       rs11589256      214196749       271
15      rs10401369      19268718        232
11      rs10873487      767334548       964
 
### The paste command is used to align files side-by-side, merging records from each file respectively.
$ paste table.txt list1.txt | column -t
CHR  SNP         BP         A1  C_A  C_U   A2  CHISQ    P       OR      apples
19   rs10401969  19268718   C   222  890   T   0.03462  0.8524  0.9857  bananas
1    rs10873883  76734548   G   934  3811  A   0.5325   0.4656  0.9691  plums
1    rs11589256  214196749  C   271  1084  T   0.01928  0.8896  0.9902  carrots
15   rs10401369  19268718   C   232  890   T   0.03232  0.2524  0.1157
11   rs10873487  767334548  G   964  3811  A   0.5525   0.2356  0.2391
 
$ cat p1.txt           ### See the contents of file p1.txt
IBM
MSFT
Apple
SAP
Yahoo
 
$ cat p2.txt          ### See the contents of file p2.txt
25.23
234.02
23.03
11.22
15.8
 
$ paste p1.txt p2.txt   ### Note what happens when we paste p1.txt and p2.txt
IBM     25.23
MSFT    234.02
Apple   23.03
SAP     11.22
Yahoo   15.8
 
> paste -d ":" p1.txt p2.txt  ###P aste the files delimited by ':'  
IBM:25.23
MSFT:234.02
Apple:23.03
SAP:11.22
Yahoo:15.8
```

f. Count the number of words, lines and bytes with 'wc'
```bash
$ wc -l list1.txt list2.txt              ###Count the number of lines in list1.txt and list2.txt
 4 list1.txt
 5 list2.txt
 9 total
 
$ wc -w list2.txt                        ###Count the number of words in list2.txt
10 list2.txt
 
$ wc list2.txt                           ###There are 5 lines, 10 words and 58 bytes in list2.txt
 5 10 58 list2.txt
```
:bulb: Download and Review the [LinuxReference.pdf](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/LinuxReference.pdf) file, a compilation of basic and most useful Linux comand for bioinformatics 

## Shell scripting

*Why shell scripting?* if you need to repeat a process 1000 times, you could either retype the commands 1000 times or you can create a small 'script' (set of instructions) that tells the system to repeat it 1000 times.

Let's create a simple shell script, **on your comand line** open a text editor ($ nano hello.sh), type the following text on the editor:

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
chmod u+x hello.sh         ### Make your shell script executable
 
./hello.sh                 ### Execute your first shell script
ls file*                   ### List the files you just created 
```

## Using Linux command line in Bioinformatics

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
cp /group/bioinformatics/shared/MSIB32500/Lecture1/linux/SRR001655.fastq . 
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
```
```bash
tail -40000 SRR001655.fastq > bottom_10000.fastq
```
**b. Using the 'paste' command to format your data**

Showing fastq formated data as a table (i.e. in columns) can be very useful to explore the data, the 'paste' command writes lines in a file as columns separated by a the tab character. The command take character '-' as an option to represent the standard input, e.g.: the option '- - - -', will be trasalated as 'read four lines', and write them out as four columns:

```bash
cat bottom_10000.fastq | paste - - - - | head -10
```
```bash
cat bottom_10000.fastq | paste - - - - | head -1000 > top_1000_tab.txt 
cat top_1000_tab.txt 
```
**c. Using 'awk' to work with data in columns**

The linux command 'awk' is very useful and practical for text manipulation in bioinformatics, 'awk' works with data in tabular format (like the result files on the previous excersice). The name stands for Aho, Weinberger and Kernighan (yes, [Brian Kernighan](https://www.cs.princeton.edu/~bwk/)), the authors of the language, which started in 1977.

What is it that awk does?

awk is a utility/language designed for data extraction awk is often used with 'sed' to perform useful and practical text manipulation tasks in bioinformatics. One of the most simple and popular uses of 'awk' is selecting a column from a text file, or other command's output. 

The general syntax of awk is:

```
awk '/search pattern/{Actions}' filename
```

'awk' goes through each line in *filename* and if the line matches the *search pattern*, then action(s) is performed

```bash
awk '/N/{print}' top_1000_tab.txt
```
This will prints out all lines in the file that contain an N (i.e. when you want to know if your sequence failed!!)

Now explore the result of the following command:

```bash
awk '/N/{print $1,"\t",$2,"\t",$3,"\t",$4}' top_1000_tab.txt
```
Note the effect of $1, $2, $3 and $4 (change the order in the command). The '\t' is a 'scape' sequence the ‘\t’ is an escape sequence, used to represent a tab delimited. Other scape sequences are: '\n' = new line; '\r' = carriage return; '\\' = a literal backslash.

Now explore the result of the following command:

```bash
awk $3 '/N/{print $1,"\t",$3}' top_1000_tab.txt
```
For each line, search for 'N' in the 3th. column, if N is found, print the first and thirth column separed by a tab. Now try:

```bash
awk $2 '/389/{print $1,"\t",$2}' top_1000_tab.txt
```
Yes, you guessed right, search for 389 in the second column and print columns 1 and 2.

Here is an example on how to count the number of fail reads in our sequence file:
```bash
awk '$3 ~ /N/ {print $1}' top_1000_tab.txt | wc -l
```
The follwong command will convert FASTQ to FASTA file format:
```bash
cat SRR001655.fastq | paste - - - - | awk '{print ">"$1,$2,"\n"$3}'
```
Now try saving the result of the conversion on a new file.

**d. Text manipulation with 'sed' **
 
sed stands for **s**tream **ed**itor is a stream oriented editor which was created exclusively for executing scripts. Thus all the input you feed into 'sed' passes through and goes to the screen (STDOUT). In other words,'sed' does not change the input file.

The general syntax for sed is:

```
/pattern/action
```





