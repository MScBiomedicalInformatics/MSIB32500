# Introduction to Unix/Linux: 

**Center for Research Informatics, University of Chicago**

Saturdays 10/06/18; 9:00AM - 12:00PM

**Instructor:** Jorge Andrade, Ph.D.

## Learning Objectives

- This training will teach you the basics on how to start using Linux, be it on the CRI’s High-Performance Infrastructure (HPC), or on your own computer.
  - After this training you will: 
  - Feel comfortable using Linux
  - Know how software works on Linux and how to use it
  - Use bash to execute commands in Linux and know your way around the file system
  - Have an overview of Linux tools that are extremely useful for bioinformatics

## Log on to CRI's High-Performance Computing (HPC) system

:pushpin:**Microsoft Windows user** you will need to install a tool for remote computing: [MobaXterm](http://mobaxterm.mobatek.net) and/or
[PuTTY](http://www.putty.org). PuTTY user go to http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html locate the appropriate binary file (executable putty.exe file) for your hardware (32bits or 64 bits laptop). Place this file on your desktop (or another folder that you can access easily). You will need to 'double-click' this file each time you need to access it.

Open a PuTTY session, on 'Host Name (or IP address)' type: **gardner.cri.uchicago.edu**,  select SSH as the Connection Type, verify the port number in the 'Port' Box is **22**. Press the 'Open' button, then type in the provided **username and password** when prompted. Type **yes** if you are prompted to accept a key before entering a username.

**Login**

For MacOS or Unix/Linux users:
1. Open a terminal session.
2. Connect to the login node of GARDNER cluster:

```bash
ssh t.cri.biowksp01@gardner.cri.uchicago.edu
```
Enter your **_password_** when prompted. Type yes if you are prompted to accept a key.


## Accessing CRI's HPC system from off-campus locations

In order to securely access CRI's HPC system from off-campus locations, you need to set up a [virtual private network (VPN)](https://en.wikipedia.org/wiki/Virtual_private_network), to do this, go to: https://vpn.uchicago.edu and download the University's CISCO AnyConnect Secure Mobility Client (available for all operative systems), after downloading and installing the client, use **cvpn.uchicago.edu** as the server's name.

![cvpn](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/cvpn.uchicago.edu.png)

## How to get help

- Use the manual ($ man) command; to exit the manual type 'q'
- Ask for help ($ your_command --help)
- Use the command apropos ($ apropos text)

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
ls -a         ### List all the files in your current directory (do not ignore entries starting with .)
```
b. Change current directory to the root directory of the file system and explore the directory structure
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
cd ~
mkdir mydir  ### Make a new directory called 'mydir'
mkdir newfolder1
ls -l 
```
b. Create a new file in a directory
```bash
cd mydir
nano file1.txt        ### Use nano editor to create a new file, write some text and use Control-O to save and Control-X to exit.
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
pwd
```
d. Remove a file or directory

```bash
cd ~
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
cd /group/mscbmi/lecture1/
ls -l SRR*                     ### List all files that start with SRR
cp SRR001655.fastq ~/
cd ~/
gzip SRR001655.fastq           ### Compress a file, this command should create a compressed file named SRR001655.fastq.gz
ls -l                        
gunzip SRR001655.fastq.gz      ### Decompress a file
ls -l
```
h. Compress files protected by a password

```bash
zip -er SRR.zip SRR001655.fastq 
```
This will prompt you for a password. 
-e enables encryption for your zip file. This is what makes it ask for the password.
-r makes the command recursive, meaning that all the files inside a folder will be added to the zip file.

i. Unzip a file

```bash
unzip SRR.zip 
```

## :mortar_board:Your turn: File structure challange

1. Create the following file structure

![filestructure](https://raw.githubusercontent.com/MScBiomedicalInformatics/MSIB32500/master/cheatsheets/files.jpg)

2. Copy the file *SRR001655.fastq* from */group/mscbmi/lecture1/* to the *folder4*
3. Move a password protected encrypted copy of SRR001655.fastq from *folder4* to *folder1* 
4. Remove recursively the contents of the folder *ex1*

## :clap: Well done!


## File transfer between computers

a. Windows user download and install [WinSCP](http://winscp.net/eng/index.php). MacOS users open the Terminal

The *scp* (secure copy command) allows you to copy/move files between compters on the command line
```bash
scp example1.txt username@gardner.cri.uchicago.edu:.
```
We will now use the scp command to transfer a file from your local computer to your **home** directory on CRI's [GARDNER](https://en.wikipedia.org/wiki/Martin_Gardner) infrastructure. 

First, download the file available at the following link, to your local computer:

https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/data/GSE31736_RAW.tar

:red_circle: Using your computer's command line (**open a new terminal**), navigate to the directory where file 'GSE31736_RAW.tar' is located (use the **cd** command).

Now use the **scp** command with your username and password, to securely copy the file. Your command should look like:

```bash
scp ./GSE31736_RAW.tar t.cri.biowksp01@gardner.cri.uchicago.edu:~/
```

:bulb: The [hypexr.org](http://www.hypexr.org/linux_scp_help.php) website has a list of examples on how to use secure copy. 

b. You can use *wget* command to get a file from the internet directly to your working directory in Linux

```bash
cd ~
wget http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
```
## Input/Output redirect and pipe

a. Use the symbol '>' to redirect the output of a command to a file

```bash
cd ~
nano text1.txt                       ### Create a new file and write some text on it
 
cat text1.txt                        ### Print file1.txt to the screen
 
cat text1.txt > text2.txt            ### Print file1.txt to file2.txt
 
nano text2.txt                       ### Use nano to view file2.txt

```
b. Use the symbol '>>' to redirect the output of a command to a file using 'append'

```bash
cat text1.txt >> text2.txt           ### Append text1.txt to text2.txt
nano text2.txt                       ### View text2.txt
```
c. Use the symbol '<' to redirect the input
```bash
cat < text1.txt                      ### Print file1.txt to screen
```
d. Use the symbol '|' to pipe the output of one command as the input of another 

```bash
ls | wc -l              ### List the files and folders in the current directory and then count them
```
## Text extraction and manipulation

A standard Unix/Linux installation will have available several text editors (like: vi, vim, nano, emacs, and others) and text viewers (like: less, more, head and tail)

a. Use the symbol '>' to redirect the output of a command to a file

```bash
cd /group/mscbmi/lecture1/              
less SRR001655.fastq                 ### View a fastq file: use q to exit, space or f to the next page,
                                     ### b to the previous page, and / to search a word
head -20 SRR001655.fastq             ### Show the first 20 line of the file
tail -20 SRR001655.fastq             ### Show the last 20 lines of the file

```
b. Pattern Search with 'grep'. The grep command searches specified files or other input (stdin) for patterns matching a given expression(s).

```bash
cp /group/mscbmi/lecture1/list1.txt ~/
cp /group/mscbmi/lecture1/list2.txt ~/
cd ~
ls -l

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
 
 
$ cat list2.txt | sed 's/bean/button/g'       ### Replace bean with button *g*lobally, try also sed 's/bean/button/2'
Apple Sauce
wild rice
black buttons
kidney buttons
dry apples

$ cat list2.txt  ### why there is no changes in list2.txt?

$ cat list2.txt | tr /a-z/ /A-Z/              ### Change lower case to upper case (tr: translate)
thegeekstuff
THEGEEKSTUFF
APPLE SAUCE
WILD RICE
BLACK BEANS
KIDNEY BEANS
DRY APPLES
 
$ cat list2.txt 

$ echo 'attgctcgat' | tr /atgc/ /tacg/ | rev  ### Find a complement DNA sequence and then reverse it
atcgagcaat
```

e. Table manipulation: sort, uniq, cut, and paste

```bash
cp /group/mscbmi/lecture1/table.txt ~/
cd ~
ls -l
$ cat table.txt                                      ### show the contents of file table.txt
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
 ls -l
$ cat table_sorted1.txt | column -t
CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
 
$ sort -k8n table.txt  > table_sorted2.txt        ### sort by CHISQ value 
 
$ sort -k2 table.txt                               ### sort the table by the second column
 
  15   rs10401369   19268718    C      232      890    T      0.03232       0.2524       0.1157
  19   rs10401969   19268718    C      222      890    T      0.03462       0.8524       0.9857
  11   rs10873487  767334548    G      964     3811    A       0.5525       0.2356       0.2391
   1   rs10873883   76734548    G      934     3811    A       0.5325       0.4656       0.9691
   1   rs11589256  214196749    C      271     1084    T      0.01928       0.8896       0.9902
  13   rs11589552 2014196749    C      221     1184    T      0.01878       0.8796       0.1202
 CHR          SNP         BP   A1      C_A      C_U   A2        CHISQ            P           OR
 
 
$ cut -f1,2,3,5 table.txt | column -t               ####Extract columns (*f*ields) 1, 2, 3, and 5 from table.txt
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
 
cp /group/mscbmi/lecture1/p1.txt ~/
$ cat p1.txt           ### See the contents of file p1.txt
IBM
MSFT
Apple
SAP
Yahoo

cp /group/mscbmi/lecture1/p2.txt ~/
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
 
$ paste -d ":" p1.txt p2.txt  ### Paste the files delimited by ':'  
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

### Linux power tools: awk and sed

**awk** The linux command 'awk' is very useful and practical for text manipulation in bioinformatics, 'awk' works with data in **tabular format.** The name stands for Aho, Weinberger and Kernighan [Brian Kernighan](https://www.cs.princeton.edu/~bwk/)), the authors of the language, which started in 1977.

What is it that awk does?

**awk** is a utility/language designed for data extraction awk is often used with 'sed' to perform useful and practical text manipulation tasks in bioinformatics. One of the most simple and popular uses of 'awk' is selecting a column from a text file, or other command's output. 

The general syntax of awk is:

```
awk '/search pattern/{Actions}' filename
```

'awk' goes through each line in *filename* and if the line matches the *search pattern*, then action(s) is performed. For example, to solve the problem of the headers on a table (table.txt) being sorted whe we used the 'sort' command in the previous section, one can use the follwoing command:

```
awk 'NR==1; NR>1 {print $0 | "sort -k2"}' table.txt
```

Let's see more awk examples. Copy the file *mappingtools.txt* from /group/mscbmi/lecture2.txt. This is a tab delimited file containing a list of common mapping tools for NGS data.

```bash
cd ~
cp /group/mscbmi/lecture2/mappingtools.txt . 
```
To explore the contents of mappingtools.txt we can use:

```bash
cat mappingtools.txt | column -t 
```

We can use awk to select and print one column on the tab delimited file:

```bash
awk '{print $2}' mappingtools.txt
```

To select and print one row (line):

```bash
awk 'FNR == 2 {print}' mappingtools.txt
```

Find "BWA" in the first column:

```bash
$ awk '$1 ~ /BWA/' mappingtools.txt
BWA	Free_software	Illumina	Command_line	Output_adhere_to_SAM_format
BWA-SW	Free_software	454_IonTorrent_Illumina	Command_line	For_long_sequences_ranged_from_70bp_to_1Mbp
BWA-MEM	Free_software	454_IonTorrent_Illumina	Command_line	For_long_sequences_ranged_from_70bp_to_1Mbp
```

Find the exact string matches (use " "):

```bash
$ awk '$1== "BWA"' mappingtools.txt
BWA	Free_software	Illumina	Command_line	Output_adhere_to_SAM_format
```

Find "BWA" in the complete file:

```bash
$ awk '/BWA/' mappingtools.txt
BWA	Free_software	Illumina	Command_line	Output_adhere_to_SAM_format
SHRiMP2	Free_software	Illumina	Command_line	Higher_sensitivity_than_BWA
BWA-SW	Free_software	454_IonTorrent_Illumina	Command_line	For_long_sequences_ranged_from_70bp_to_1Mbp
BWA-MEM	Free_software	454_IonTorrent_Illumina	Command_line	For_long_sequences_ranged_from_70bp_to_1Mbp
```

If an exact match to the string "BWA" is found on the first collumn; then print the 3rd collumn

```bash
$ awk '$1== "BWA" {print $3}' mappingtools.txt
Illumina
```

How to omit the header for the column and get only the names printed?

```bash
awk 'NR!=1{print $1}' mappingtools.txt
```

Let's see which mapping tool can be used with Illumina technology:

```bash
awk '/Illumina/{print$1}' mappingtools.txt
```

Now, give me the name of the Software and the Notes for all mappers that can use 454 tecnology 

```bash
awk '/454/{print $1,"\t",$5}' mappingtools.txt
```
Note the effect of $1, and $5 (test to change the order in the command). 

The '\t' is a 'scape' sequence, used to represent a tab delimited. Other scape sequences are: '\n' = new line; '\r' = carriage return; '\\' = a literal backslash.

Now explore the result of the following command:

```bash
awk $5 '/Fast/{print $1,"\t",$2, "\t",$3}' mappingtools.txt
```

### More awk examples for bioinformatics:

Copy the file expression.txt from our repository to your home directory; explore the contents fo this file.

```bash
cp /group/mscbmi/lecture1/expression.txt .
head -20 expression.txt 
wc -l expression.txt

```
1. Choose rows where the gene expression of tumor replicate 1 is larger than the expression of normal replicate 1. Save the results on a file called expoutput.txt

```bash

awk '$2>$5' expression.txt > expoutput.txt
head -20 expoutput.txt 
wc -l expoutput.txt

```

2. Extract the expression value for the 3 normal replicates, columns 5,6,7. Save it on a text file

```bash

awk '{print $5,$6,$7}' expression.txt > tumor.txt
head -20 tumor.txt 

```

3. Extract the expression values for gene Rock2 on all samples

```bash
awk '/Rock2/{print}' expression.txt
```

4. Calculate the average expression value for all genes of normal replicate 3 (column 7)

```bash
awk '{x+=$7}END{print x/NR}' expression.txt
```

--------------
 
**sed** stands for **s**tream **ed**itor and is a stream oriented editor which was created exclusively for executing scripts. Thus all the input you feed into 'sed' passes through and goes to the screen (STDOUT). In other words,'sed' does not change the input file.

The general syntax for sed is:

```bash
sed /pattern/action
```
Where 'pattern' is a regular expression, and action is one of the following: 

- 'p'= Prints the line
- 'd'= Deletes the line
- 's/pattern1/pattern2/' = Substitutes the first occurrence of pattern1 with pattern2. 

If 'pattern' is omitted, action is performed for every line.

**Sed Command Examples**

Replacing or substituting string

```bash
cat mappingtools.txt | column -t
sed 's/Fast/Ultra_Fast/' mappingtools.txt | column -t
cat mappingtools.txt | column -t
```
The **sed** command also works on plain text files:

```bash
$ cd ~
$ cp /group/mscbmi/lecture2/bwa.txt . 
$ cat bwa.txt
Burrows-Wheeler Aligner


BWA is a software package for mapping low-divergent sequences against a large
reference genome, such as the human genome. It consists of three algorithms: 
BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for 
Illumina sequence reads up to 100bp, while the rest two for longer sequences 
ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as 
long-read support and split alignment, but BWA-MEM, which is the latest, is 
generally recommended for high-quality queries as it is faster and more accurate. 
BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina 
reads.
```
Replacing the *n*th occurrence of a pattern in a line

```bash
$ sed 's/BWA-/Burrows-Wheeler-Aligner/2' bwa.txt
Burrows-Wheeler Aligner


BWA is a software package for mapping low-divergent sequences against a large
reference genome, such as the human genome. It consists of three algorithms: 
BWA-backtrack, Burrows-Wheeler-AlignerSW and BWA-MEM. The first algorithm is designed for 
Illumina sequence reads up to 100bp, while the rest two for longer sequences 
ranged from 70bp to 1Mbp. BWA-MEM and Burrows-Wheeler-AlignerSW share similar features such as 
long-read support and split alignment, but BWA-MEM, which is the latest, is 
generally recommended for high-quality queries as it is faster and more accurate. 
BWA-MEM also has better performance than Burrows-Wheeler-Alignerbacktrack for 70-100bp Illumina 
reads.
```

Replacing all occurrences of the pattern in a line

```bash
sed 's/BWA-/Burrows-Wheeler-Aligner/g' bwa.txt
```

Using **sed** to delete a especific line

```bash
sed '5d' bwa.txt
```

Using **sed** to delete several lines

```bash
sed '2d;3d;4d;5d' bwa.txt
```

Using **sed** to print a especific line

```bash
sed '1p;7p' bwa.txt
```
________________________________________
## Linux File Permissions

**Permission Groups**

Each file and directory has three user based permission groups:

- **owner** - The Owner permissions apply only the owner of the file or directory, they will not impact the actions of other users.
- **group** - The Group permissions apply only to the group that has been assigned to the file or directory, they will not effect the actions of other users.
- **all users** - The All Users permissions apply to all other users on the system, this is the permission group that you want to watch the most.

**Permission Types**

Each file or directory has three basic permission types:

- **read** - The Read permission refers to a user's capability to read the contents of the file.
- **write** - The Write permissions refer to a user's capability to write or modify a file or directory.
- **execute** - The Execute permission affects a user's capability to execute a file or view the contents of a directory.

You can view the permissions by checking the file or directory permissions with the **ls -l** command.

The permission in the command line is displayed as: **_rwxrwxrwx 1 owner:group**

- The first character, marked with an *underscore* is the special permission flag that can vary. For a directory, that space will be marked with a *d*.
- The following set of three characters (rwx) is for the **owner** permissions
- The second set of three characters (rwx) is for the **Group** permissions
- The third set of three characters (rwx) is for the **All Users** permissions.

**Explicitly Defining Permissions**

For explicit definition of permissions, linux uses the following nomenclature:

- u - User (Owner)
- g - Group
- o - Others
- a - All users

**+** (plus) and **-** (minus) are used to tell the system whether to **add** or **remove** the specific permissions.

The Permission Types that are used are:

- r - Read
- w - Write
- x - Execute

To make this modification you would use the command: **chmod**. Let's see some examples:

Enable Read, Write and Execute for all users on file list1.txt

```bash
$ ls -l list1.txt
-rwxr-x--- 1 t.cri.biowksp01 t.cri.biowksp01 29 Mar 31 11:16 list1.txt

$ chmod a+rwx list1.txt
ls -l list1.txt
-rwxrwxrwx 1 t.cri.biowksp01 t.cri.biowksp01 29 Mar 31 11:16 list1.txt
```

Remove Read and Write permission for Others on file list1.txt

```bash
chmod o-rw list1.txt
ls -l list1.txt
-rwxrwx--x 1 t.cri.biowksp01 t.cri.biowksp01 29 Mar 31 11:16 list1.txt
```

Enable Execute permission for Group on list1.txt  

```bash
chmod g+x list1.txt
```

**Using Binary References to Set permissions**

The setting is done by entering three integers numbers. The first number represents the **Owner** permission; the second represents the **Group** permissions; and the last number represents the permissions for **all other** users. The numbers are a binary representation of the **rwx** string

- r = 4
- w = 2
- x = 1

You add the numbers to get the integer/number representing the permissions you wish to set. 

Let's see some examples:

```bash
$ ls -l list2.txt
-rwxr-x--- 1 t.cri.biowksp01 t.cri.biowksp01 58 Mar 31 11:16 list2.txt
$ chmod 740 list2.txt
$ ls -l list2.txt
-rwxr----- 1 t.cri.biowksp01 t.cri.biowksp01 58 Mar 31 11:16 list2.txt
```
This command enable **rwx** (7) to the user/owner, only **r** (4) to the group and no right to **rwx** for all other user.

Another example:

```bash
chmod 777 list2.txt
ls -l
```


## Week 1 Challange: :white_check_mark: 

:bulb: Download the [LinuxReference.pdf](https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/LinuxReference.pdf) file, review and practice at home.

- Follow the tutorial: [CloudEnviromentSetup.pdf] (https://github.com/MScBiomedicalInformatics/MSIB32500/blob/master/cheatsheets/CloudEnviromentSetup.pdf) to learn how to confirure your owm Elastic Cloud Computing instance (EC2) using  Amazon Web Services (AWS). 
- Send me an e-mal with and screenshot of your EC2 instance running. 

## Week 1 Suguested reading: :white_check_mark:
Read the Nature technology feature article **'Biology: The big challenges of big data'** available at: http://www.nature.com/nature/journal/v498/n7453/full/498255a.html 

Some [computer comics and cartoons](http://www.hypexr.org/comics.php)

    

