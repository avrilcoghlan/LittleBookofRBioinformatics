Answers to the exercises in chapter 1
=====================================                 

Q1.
---

*What are the last twenty nucleotides of the Bacteriophage lambda genome sequence?*

To answer this, you first need to install the "SeqinR" R library, and download
the Bacteriophage lambda genome sequence from the NCBI database and save it as
a file "lambda.fasta" in the "My Documents" folder. 

Then to find the length of the Bacteriophage lambda genome sequence, type in the R console:

::

    > library("seqinr")
    > lambda <- read.fasta(file="lambda.fasta")
    > lambdaseq <- lambda[[1]]
    > length(lambdaseq) 
    [1] 48502

This tells us that the length of the sequence is 48502 nucleotides.
Therefore, the last 20 nucleotides are from 48483 to 48502. You can
extract the sequence of these nucleotides by typing:

::

    > lambdaseq[48483:48502]
    [1] "c" "g" "g" "t" "g" "a" "t" "c" "c" "g" "a" "c" "a" "g" "g" "t" "t"
    [18] "a" "c" "g"




