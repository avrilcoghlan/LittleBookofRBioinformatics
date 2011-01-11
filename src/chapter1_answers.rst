Answers to the exercises on DNA Sequence Statistics
===================================================   

Q1. *What are the last twenty nucleotides of the DEN-1 Dengue virus genome sequence?*

To answer this, you first need to install the "SeqinR" R library, and download
the DEN-1 Dengue genome sequence from the NCBI database and save it as
a file "den1.fasta" in the "My Documents" folder. 

Then to find the length of the DEN-1 Dengue virus genome sequence, type in the R console:

.. highlight:: r

::

    > library("seqinr")
    > dengue <- read.fasta(file="den1.fasta")
    > dengueseq <- dengue[[1]]
    > length(dengueseq) 
    [1] 10735

This tells us that the length of the sequence is 10735 nucleotides.
Therefore, the last 20 nucleotides are from 10716 to 10735. You can
extract the sequence of these nucleotides by typing:

::

    > dengueseq[10716:10735]
    [1] "c" "t" "g" "t" "t" "g" "a" "a" "t" "c" "a" "a" "c" "a" "g" "g" "t" "t" "c"
    [20] "t"
    
Q2. *What is the length in nucleotides of the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC\_002677)?*

To get the *Mycobacterium leprae* TN genome, and find out its length, it’s necessary to first go to the NCBI website (www.ncbi.nlm.nih.gov) and search for NC\_002677 and download it as a fasta format file (eg. "leprae.fasta") and save it in the "My Documents" folder. Then in R type:

::

    > leprae <- read.fasta(file="leprae.fasta")
    > lepraeseq <- leprae[[1]]
    > length(lepraeseq)
    [1] 3268203

Q3. *How many of each of the four nucleotides A, C, T and G, and any other symbols, are there in the Mycobacterium leprae TN genome sequence?*

Type:

::

    > table(lepraeseq)
    lepraeseq
         a      c      g      t 
    687041 938713 950202 692247 

Q4. *What is the GC content of the Mycobacterium leprae TN genome sequence, when (i) all non-A/C/T/G nucleotides are included, (ii) non-A/C/T/G nucleotides are discarded?*  

Find out how the GC function deals with non-A/C/T/G nucleotides, type:

::

    > help("GC")

Type:

::

    > GC(lepraeseq)
    [1] 0.5779675
    > GC(lepraeseq, exact=FALSE)
    [1] 0.5779675

This gives 0.5779675 or 57.79675%. This is the GC content when non-A/C/T/G nucleotides are not taken into account.  

The length of the *M. leprae* sequence is 3268203 bp, and it has 938713 Cs and 950202 Gs, and 687041 As and 692247 Ts. So to calculating the GC content when only considering As, Cs, Ts and Gs, we can also 
type:

::

    > (938713+950202)/(938713+950202+687041+692247)
    [1] 0.5779675

To take non-A/C/T/G nucleotides into account when calculating GC, type:

::

    > GC(lepraeseq, exact=TRUE)
    [1] 0.5779675

We get the same answer as when we ignored non-A/C/G/T nucleotides. This is actually because the *M. leprae* TN sequence does not have any non-A/C/G/T nucleotides. 

However, many other genome sequences do contain non-A/C/G/T nucleotides. Note that under ‘Details’ in the box that appears when you type ‘help(‘GC’)’, it says : "When exact is set to TRUE the G+C content is estimated with ambiguous bases taken into account. Note that this is time expensive. A first pass is made on non-ambiguous bases to estimate the probabilities of the four bases in the sequence. They are then used to weight the contributions of ambiguous bases to the G+C content."

Q5. *How many of each of the four nucleotides A, C, T and G are there in the complement of the Mycobacterium leprae TN genome sequence?*

First you need to search for a function to calculate reverse complement, eg. by typing:

::

    > help.search("complement")

You will find that there is a function seqinr::comp that complements a nucleic acid sequence. This means it is a function in the SeqinR library.

Find out how to use this function by typing:

::

    > help("comp")

The help says "Undefined values are returned as NA". This means that the complement of non-A/C/T/G symbols will be returned as NA.

To find the number of A, C, T, and G in the reverse complement type:

::

    > complepraeseq <- comp(lepraeseq)
    > table(complepraeseq)
     complepraeseq
          a      c      g      t 
     692247 950202 938713 687041 

Note that in the *M. leprae* sequence we had 687041 As, in the complement have 687041 Ts.
In the *M. leprae* sequence we had 938713 Cs, in the complement have 938713 Gs.
In the *M. leprae* sequence we had 950202 Gs, in the complement have 950202 Cs.
In the *M. leprae* sequence we had 692247 Ts, in the complement have 692247 As.

Q6. *How many occurrences of the DNA words CC, CG and GC occur in the Mycobacterium leprae TN genome sequence?*

::

    > count(lepraeseq, 2)
        aa     ac     ag     at     ca     cc     cg     ct     ga     gc     gg 
     149718 206961 170846 159516 224666 236971 306986 170089 203397 293261 243071 
        gt     ta     tc     tg     tt 
     210473 109259 201520 229299 152169 

Get count for CC is 236,971; count for CG is 306,986; count for GC is 293,261.

Q7. *How many occurrences of the DNA words CC, CG and GC occur in the (i) first 1000 and (ii) last 1000 nucleotides of the Mycobacterium leprae TN genome sequence?*

Type:

::

    > length(lepraeseq)
    [1] 3268203

to find the length of the *M. leprae* genome sequence.  It is 3,268,203 bp. Therefore the first 1000 nucleotides will have indices 1-1000, and the last thousand nucleotides will have indices 3267204-3268203. We find the count of DNA words of length 2 by typing:

::

    > count(lepraeseq[1:1000],2)
     aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt 
     78 95 51 49 85 82 92 54 68 63 39 43 42 73 31 54 
    > count(lepraeseq[3267204:3268203],2)
     aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt 
     70 85 44 55 94 81 87 50 53 75 49 51 36 72 48 49 

To check that the subsequences that you looked at are 1000 nucleotides long, you can type:

::

    > length(lepraeseq[1:1000])
    [1] 1000
    > length(lepraeseq[3267204:3268203])
    [1] 1000

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.


