Answers to the exercises on DNA Sequence Statistics (2)
=======================================================

Q1. *Draw a sliding window plot of GC content in the DEN-1 Dengue virus genome, using a window size of 200 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?*

To do this, you first need to download the DEN-1 Dengue virus sequence from the NCBI database. 
To do this follow the steps in the chapter `DNA Sequence Statistics (1) <./chapter1.html>`_.

Then read the sequence into R using the SeqinR library:

.. highlight:: r

::

    > library("seqinr")
    > dengue <- read.fasta(file = "den1.fasta")
    > dengueseq <- dengue[[1]]

Then write a function to make a sliding window plot:

::
   
    > slidingwindowplot <- function(windowsize, inputseq) 
      {
         starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
         n <- length(starts)    
         chunkGCs <- numeric(n)
         for (i in 1:n) { 
            chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
            chunkGC <- GC(chunk)
            chunkGCs[i] <- chunkGC 
         }  
         plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
      }

Then make a sliding window plot with a window size of 200 nucleotides: 

::

    > slidingwindowplot(200, dengueseq)

|image0|

The GC content varies from about 45% to about 50% throughout the DEN-1 Dengue virus genome, with
some noticeable troughs at about 2500 bases and at about 4000 bases along the sequence, where the
GC content drops to about 40%. There is no strong difference between the start and end of the
genome, although from around bases 4000-7000 the GC content is quite high (about 50%), and from
about 2500-3500 and 7000-9000 bases the GC content is relatively low (about 43-47%).   

We can also make a sliding window plot of GC content using a window size of 2000 nucleotides:

::

    > slidingwindowplot(2000, dengueseq)

|image1|

In this picture it is much more noticeable that the GC content is relatively high from around
4000-7000 bases, and lower on either side (from 2500-3500 and 7000-9000 bases).

Q2. *Draw a sliding window plot of GC content in the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677) using a window size of 20000 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?*

To do this, you first need to download the *Mycobacterium leprae* sequence from the NCBI 
database.
To do this follow the steps in the chapter `DNA Sequence Statistics (1) <./chapter1.html>`_.

Then read the sequence into R using the SeqinR library:

.. highlight:: r

::

    > leprae <- read.fasta(file = "leprae.fasta")
    > lepraeseq <- leprae[[1]]

Then make a sliding window plot with a window size of 20000 nucleotides:

::

    > slidingwindowplot(20000, lepraeseq) 

|image2|

We see the highest peak in GC content at about 1 Mb into the *M. leprae* genome. We also 
see troughs in GC content at about 1.1 Mb, and at about 2.6 Mb. 

With a window size of 200 nucleotides, the plot is very messy, and we cannot see the peaks and troughs
in GC content so easily:

::

    > slidingwindowplot(200, lepraeseq)

|image3|

With a window size of 200,000 nucleotides, the plot is very smooth, and we cannot see the peaks and troughs
in GC content very easily:

::

    > slidingwindowplot(200000, lepraeseq)
 
|image4|

Q3. *Write a function to calculate the AT content of a DNA sequence (ie. the fraction of the nucleotides in the sequence that are As or Ts). What is the AT content of the Mycobacterium leprae TN genome?*

Here is a function to calculate the AT content of a genome sequence:

::

    > AT <- function(inputseq)
      {
         mytable <- count(inputseq, 1) # make a table with the count of As, Cs, Ts, and Gs
         mylength <- length(inputseq) # find the length of the whole sequence
         myAs <- mytable[[1]] # number of As in the sequence
         myTs <- mytable[[4]] # number of Ts in the sequence
         myAT <- (myAs + myTs)/mylength
         return(myAT)
      }

We can then use the function to calculate the AT content of the *M. leprae* genome:

::

    > AT(lepraeseq)
    [1] 0.4220325
   
You should notice that the AT content is (1 minus GC content), ie. (AT content + GC content = 1):

::
   
    > GC(lepraeseq)
    [1] 0.5779675
    > 0.4220325 + 0.5779675
    [1] 1

Q4. *Write a function to draw a sliding window plot of AT content. Use it to make a sliding window plot of AT content along the Mycobacterium leprae TN genome, using a windowsize of 20000 nucleotides. Do you notice any relationship between the sliding window plot of GC content along the Mycobacterium leprae genome, and the sliding window plot of AT content?*

We can write a function to write a sliding window plot of AT content:

::

    > slidingwindowplotAT <- function(windowsize, inputseq) 
      {
         starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
         n <- length(starts)    
         chunkATs <- numeric(n)
         for (i in 1:n) { 
            chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
            chunkAT <- AT(chunk)
            chunkATs[i] <- chunkAT 
         }  
         plot(starts,chunkATs,type="b",xlab="Nucleotide start position",ylab="AT content")
     }

We can then use this function to make a sliding window plot with a window size of 20000 nucleotides:

::

    > slidingwindowplotAT(20000, lepraeseq)

|image5|

This is the mirror image of the plot of GC content (because AT equals 1 minus GC):
 
::

    > slidingwindowplot(20000, lepraeseq)

|image6|
  
Q5. *Is the 3-nucleotide word GAC GC over-represented or under-represented in the Mycobacterium leprae TN genome sequence?*

We can get the number of counts of each of the 3-nucleotide words by typing:

::

    > count(lepraeseq, 3)
       aaa   aac   aag   aat   aca   acc   acg   act   aga   agc   agg   agt   ata   atc   atg 
     32093 48714 36319 32592 44777 67449 57326 37409 31957 62473 38946 37470 25030 57245 44268 
       att   caa   cac   cag   cat   cca   ccc   ccg   cct   cga   cgc   cgg   cgt   cta   ctc 
     32973 52381 64102 64345 43838 64869 46037 87560 38504 78120 82057 89358 57451 29004 39954 
       ctg   ctt   gaa   gac   gag   gat   gca   gcc   gcg   gct   gga   ggc   ggg   ggt   gta 
     64730 36401 43486 61174 40728 58009 66775 80319 83415 62752 44002 81461 47651 69957 33139 
       gtc   gtg   gtt   taa   tac   tag   tat   tca   tcc   tcg   tct   tga   tgc   tgg   tgt 
     60958 65955 50421 21758 32971 29454 25076 48245 43166 78685 31424 49318 67270 67116 45595 
       tta   ttc   ttg   ttt 
     22086 43363 54346 32374

There are 61,174 GACs in the sequence. 

The total number of 3-nucleotide words is calculated by typing:

::

    > sum(count(lepraeseq,3))
    [1] 3268201

Therefore, the observed frequency of GAC is 61174/3268201 = 0.01871794.

To calculate the expected frequency of GAC, first we need to get the number of counts of 1-nucleotide words by typing:

::

    > count(lepraeseq, 1)
        a      c      g      t 
     687041 938713 950202 692247 

The sequence length is 3268203 bp.
The frequency of G is 950202/3268203 = 0.2907414.
The frequency of A is 687041/3268203 = 0.2102198.
The frequency of C is 938713/3268203 = 0.2872260.
The expected frequency of GAC is therefore 0.2907414*0.2102198*0.2872260 = 0.01755514.

The value of Rho is therefore the observed frequency/expected frequency = 0.01871794/0.01755514 = 1.066237.
That, is there are about 1.1 times as many GACs as expected. This means that GAC is slightly over-represented in this sequence.
The difference from 1 is so little however that it might not be statistically significant.

We can search for a function to calculate rho by typing:

::

    > help.search("rho")
      base::getHook                          Functions to Get and Set Hooks for Load, Attach, Detach and Unload
      seqinr::rho                            Statistical over- and under- representation of dinucleotides in a sequence
      stats::cor.test                        Test for Association/Correlation Between Paired Samples
      survival::pbc                          Mayo Clinic Primary Biliary Cirrhosis Dat

There is a function rho in the SeqinR library. For example, we can use it to calculate Rho for 
words of length 3 in the *M. leprae* genome by typing:

::

    > rho(lepraeseq, wordsize=3)
           aaa       aac       aag       aat       aca       acc       acg       act       aga 
      1.0570138 1.1742862 0.8649101 1.0653761 1.0793820 1.1899960 0.9991680 0.8949893 0.7610323 
           agc       agg       agt       ata       atc       atg       att       caa       cac 
      1.0888781 0.6706048 0.8856096 0.8181874 1.3695545 1.0462815 1.0697245 1.2626819 1.1309452 
           cag       cat       cca       ccc       ccg       cct       cga       cgc       cgg 
      1.1215062 1.0487995 1.1444773 0.5944657 1.1169725 0.6742135 1.3615987 1.0467726 1.1261261 
           cgt       cta       ctc       ctg       ctt       gaa       gac       gag       gat 
      0.9938162 0.6939044 0.6996033 1.1197319 0.8643241 1.0355868 1.0662370 0.7012887 1.3710523 
           gca       gcc       gcg       gct       gga       ggc       ggg       ggt       gta 
      1.1638601 1.0246015 1.0512300 1.0855155 0.7576632 1.0266049 0.5932565 1.1955191 0.7832457 
           gtc       gtg       gtt       taa       tac       tag       tat       tca       tcc 
      1.0544820 1.1271276 1.1827465 0.7112314 0.7888126 0.6961501 0.8135266 1.1542345 0.7558461 
           tcg       tct       tga       tgc       tgg       tgt       tta       ttc       ttg 
      1.3611325 0.7461477 1.1656391 1.1636701 1.1469683 1.0695410 0.7165237 1.0296334 1.2748168 
           ttt 
      1.0423929 
    
The Rho value for GAC is given as 1.0662370, in agreement with our calculation above.

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.

.. |image0| image:: ../_static/A2_image0.png
.. |image1| image:: ../_static/A2_image1.png
.. |image2| image:: ../_static/A2_image2.png
.. |image3| image:: ../_static/A2_image3.png
.. |image4| image:: ../_static/A2_image4.png
.. |image5| image:: ../_static/A2_image5.png
.. |image6| image:: ../_static/A2_image6.png

