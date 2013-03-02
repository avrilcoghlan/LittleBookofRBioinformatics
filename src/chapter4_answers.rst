Answers to the exercises on Sequence Alignment 
==============================================

Q1. 
---
*Download FASTA-format files of the Brugia malayi Vab-3 protein (UniProt accession A8PZ80) and the Loa loa Vab-3 protein (UniProt accession E1FTG0) sequences from UniProt.*

We can use SeqinR to retrieve these sequences by typing:

::

    > library("seqinr")                           # load the SeqinR package
    > choosebank("swissprot")                     # select the ACNUC sub-database to be searched
    > query("brugia", "AC=A8PZ80")                # search for the Brugia sequence
    > brugiaseq <- getSequence(brugia$req[[1]])   # get the Brugia sequence
    > query("loa", "AC=E1FTG0")                   # search for the Loa sequence
    > loaseq <- getSequence(loa$req[[1]])         # get the Loa sequence
    > closebank()                                 # close the connection to the ACNUC sub-database

Q2. 
---
*What is the alignment score for the optimal global alignment between the Brugia malayi Vab-3 protein and the Loa loa Vab-3 protein, when you use the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?*

We can use the Biostrings R package to answer this, by typing:

::

    > library("Biostrings")                       # load the Biostrings package
    > data(BLOSUM50)                              # load the BLOSUM50 scoring matrix
    > brugiaseqstring <- c2s(brugiaseq)           # convert the Brugia sequence to a string 
    > loaseqstring <- c2s(loaseq)                 # convert the Loa sequence to a string
    > brugiaseqstring <- toupper(brugiaseqstring) # convert the Brugia sequence to uppercase
    > loaseqstring <- toupper(loaseqstring)       # convert the Loa sequence to a string
    > myglobalAlign <- pairwiseAlignment(brugiaseqstring, loaseqstring, substitutionMatrix = "BLOSUM50", 
      gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE) # align the two sequences
    > myglobalAlign
      Global PairwiseAlignedFixedSubject (1 of 1)
      pattern: [1] MK--LIVDSGHTGVNQLGGVFVNGRPLPDSTRQKI...IESYKREQPSIFAWEIRDKLLHEKVCSPDTIPSA 
      subject: [1] SSSNLFADSGHTGVNQLGGVFVNGRPLPDSTRQKI...IESYKREQPSIFAWEIRDKLLHEKVCSPDTIPSV 
      score: 777.5

The alignment score is 777.5.

Q3. 
---
*Use the printPairwiseAlignment() function to view the optimal global alignment between Brugia malayi Vab-3 protein and the Loa loa Vab-3 protein, using the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5.*

To do this, first you must copy and paste the printPairwiseAlignment() function into R. 

Then you can use it to view the alignment that you obtained in Q2:

::

    > printPairwiseAlignment(myglobalAlign)
      [1] "MK--LIVDSGHTGVNQLGGVFVNGRPLPDSTRQKIVDLAHQGARPCDISRILQVSNGCVS 58"
      [1] "SSSNLFADSGHTGVNQLGGVFVNGRPLPDSTRQKIVDLAHQGARPCDISRILQVSNGCVS 60"
      [1] " "
      [1] "KILCRYYESGTIRPRAIGGSKPRVATVSVCDKIESYKREQPSIFAWEIRDKLLHEKVCSP 118"
      [1] "KILCRYYESGTIRPRAIGGSKPRVATVSVCDKIESYKREQPSIFAWEIRDKLLHEKVCSP 120"
      [1] " "
      [1] "DTIPSA 178"
      [1] "DTIPSV 180"
      [1] " "
     
The two proteins are very similar over their whole lengths, with few gaps and mostly identities (few mismatches).

Q4. 
---
*What global alignment score do you get for the two Vab-3 proteins, when you use the BLOSUM62 alignment matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?*

Again, we can use the Biostrings R package to answer this, by typing:

::

    > data(BLOSUM62)                              # load the BLOSUM62 scoring matrix
    > myglobalAlign2 <- pairwiseAlignment(brugiaseqstring, loaseqstring, substitutionMatrix = "BLOSUM62", 
      gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE) # align the two sequences
    > myglobalAlign2
      Global PairwiseAlignedFixedSubject (1 of 1)
      pattern: [1] MK--LIVDSGHTGVNQLGGVFVNGRPLPDSTRQKI...IESYKREQPSIFAWEIRDKLLHEKVCSPDTIPSA 
      subject: [1] SSSNLFADSGHTGVNQLGGVFVNGRPLPDSTRQKI...IESYKREQPSIFAWEIRDKLLHEKVCSPDTIPSV 
      score: 593.5 

The alignment score when BLOSUM62 is used is 593.5, while the score when BLOSUM50 is used is 777.5 (from Q2). 

We can print out the alignment and see if the alignment made using BLOSUM62 is different from that
when BLOSUM50 is used:

::

    > printPairwiseAlignment(myglobalAlign2)
      [1] "MK--LIVDSGHTGVNQLGGVFVNGRPLPDSTRQKIVDLAHQGARPCDISRILQVSNGCVS 58"
      [1] "SSSNLFADSGHTGVNQLGGVFVNGRPLPDSTRQKIVDLAHQGARPCDISRILQVSNGCVS 60"
      [1] " "
      [1] "KILCRYYESGTIRPRAIGGSKPRVATVSVCDKIESYKREQPSIFAWEIRDKLLHEKVCSP 118"
      [1] "KILCRYYESGTIRPRAIGGSKPRVATVSVCDKIESYKREQPSIFAWEIRDKLLHEKVCSP 120"
      [1] " "
      [1] "DTIPSA 178"
      [1] "DTIPSV 180"
      [1] " "

The alignment made using BLOSUM62 is actually the same as that made using BLOSUM50, so it doesn't
matter which scoring matrix we use in this case.

Q5.
---
*What is the statistical significance of the optimal global alignment for the Brugia malayi and Loa loa Vab-3 proteins made using the BLOSUM50 scoring matrix, with a gap opening penalty of -10 and a gap extension penalty of -0.5?*

To answer this, we can first make 1000 random sequences using a multinomial model in which the probabilities
of the 20 amino acids are set equal to their frequencies in the *Brugia malayi* Vab-3 protein.

First you need to first copy and paste the generateSeqsWithMultinomialModel() function into R,
and then you can use it as follows:

::

    > randomseqs <- generateSeqsWithMultinomialModel(brugiaseqstring,1000) 

This makes a vector *randomseqs*, containing 1000 random sequences, each of 
the same length as the *Brugia malayi* Vab-3 protein.

We can then align each of the 1000 random sequences to the *Loa loa* Vab-3 protein, and store
the scores for each of the 1000 alignments in a vector *randomscores*:

::

    > randomscores <- double(1000) # Create a numeric vector with 1000 elements
    > for (i in 1:1000) 
      {
         score <- pairwiseAlignment(loaseqstring, randomseqs[i], substitutionMatrix = "BLOSUM50", 
           gapOpening = -9.5, gapExtension = -0.5, scoreOnly = TRUE)
         randomscores[i] <- score
      }

The score for aligning the *Brugia malayi* and *Loa loa* Vab-3 proteins using BLOSUM50 with a 
gap opening penalty of -10 and gap extension penalty of -0.5 was 777.5 (from Q2).

We can see what fraction of the 1000 alignments between the random sequences (of the same
composition as *Brugia malayi* Vab-3) and *Loa loa* Vab-3 had scores equal to or higher than 777.5:

::

    > sum(randomscores >= 777.5)
    [1] 0 

We see that none of the 1000 alignments had scores equal to or higher than 777.5.

Thus, the *p*-value for the alignment of *Brugia malayi* and *Loa loa* Vab-3 proteins is 0, and 
we can therefore conclude that the alignment score is statistically significant (as it is less than 0.05).
Therefore, it is very likely that the *Brugia malayi* Vab-3 and *Loa loa* Vab-3 proteins are
homologous (related).

Q6.
---
*What is the optimal global alignment score between the Brugia malayi Vab-6 protein and the Mycobacterium leprae chorismate lyase protein?*

To calculate the optimal global alignment score, we must first retrieve the *M. leprae* 
chorismate lyase sequence:

::

    > choosebank("swissprot")
    > query("leprae", "AC=Q9CD83")
    > lepraeseq <- getSequence(leprae$req[[1]])
    > closebank()
    > lepraeseqstring <- c2s(lepraeseq)     
    > lepraeseqstring <- toupper(lepraeseqstring)

We can then align the *Brugia malayi* Vab-3 protein sequence to the *M. leprae* chorismate
lyase sequence:

::

    > myglobalAlign3 <- pairwiseAlignment(brugiaseqstring, lepraeseqstring, substitutionMatrix = "BLOSUM50", 
      gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE) # align the two sequences
    > myglobalAlign3
      Global PairwiseAlignedFixedSubject (1 of 1)
      pattern: [1] M-----------------KLIVDSGHTGVNQLGGV...------INYAKQNNNLL----DRFILP---FSKL 
      subject: [1] MTNRTLSREEIRKLDRDLRILVATNGT-LTRVLNV...DTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL 
      score: 67.5 

The alignment score is 67.5. 

We can print out the alignment as follows:

::

    > printPairwiseAlignment(myglobalAlign3)
      [1] "M-----------------KLIVDSGHTGVNQLGGVFVNGRPLPDSTRQKIVDLAHQGARP 43"
      [1] "MTNRTLSREEIRKLDRDLRILVATNGT-LTRVLNVVANEEIVVDIINQQLLDVA-----P 54"
      [1] " "
      [1] "-------CDISRILQ---VSNGCVSKILCRYYESGTI---RPRAIGG-----SKPRVATV 85"
      [1] "KIPELENLKIGRILQRDILLKGQKSGILFVAAESLIVIDLLPTAITTYLTKTHHP-IGEI 113"
      [1] " "
      [1] "SVCDKIESYKREQ-------PSIFA----WEIRDKLLHEKVCSPDTIPSAVV-------- 126"
      [1] "MAASRIETYKEDAQVWIGDLPCWLADYGYWDL---------------PKRAVGRRYRIIA 158"
      [1] " "
      [1] "--EAIIV-----------------INYAKQNNNLL----DRFILP---FSKL 160"
      [1] "GGQPVIITTEYFLRSVFQDTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL 218"
      [1] " "

The alignment does not look very good, it contains many gaps and mismatches and few matches.

In Q5, we made a vector *randomseqs* that contains 1000 random sequences generated using a multinomial
model in which the probabilities of the 20 amino acids are set equal to their frequencies in 
the *Brugia malayi* Vab-3 protein.

To calculate a statistical significance for the alignment between *Brugia malayi* Vab-3 and
*M. leprae* chorismate lyase, we can calculate the alignment scores for the 1000 random sequences
to *M. leprae* chorismate lyase:

::

    > randomscores <- double(1000) # Create a numeric vector with 1000 elements
    > for (i in 1:1000) 
      {
         score <- pairwiseAlignment(lepraeseqstring, randomseqs[i], substitutionMatrix = "BLOSUM50", 
           gapOpening = -9.5, gapExtension = -0.5, scoreOnly = TRUE)
         randomscores[i] <- score
      }

We can then see how many of the 1000 alignment score exceed the actual alignment score for
*B. malayi* Vab-3 and *M. leprae* chorismate lyase (67.5):

::

    > sum(randomscores >= 67.5)
    [1] 22

We see that 22 of the 1000 scores for the 1000 random sequences to *M. leprae* chorismate lyase
are higher than the actual alignment score of 67.5. Therefore the *P-value* for the alignment score
is 22/1000 = 0.022. This is just under 0.05, and so is quite near to the general cutoff for statistical
significance (0.05). However, in fact it is close enough to 0.05 that we should have some doubt
about whether the alignment is statistically significant.

In fact, the *B. malayi* Vab-3 and *M. leprae* chorismate lyase proteins are not known to be 
homologous (related), and so it is likely that the relatively high alignment score (67.5) is
just due to chance alone. 

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.sanger.ac.uk/research/projects/parasitegenomics/>`_) corrections or suggestions for improvements to
my email address alc@sanger.ac.uk

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.


