Practical 4 for 2009/2010 - Sequence Alignment
==============================================

A little more about R
---------------------

In previous practicals (for example, Practical 2,
`www.ucc.ie/microbio/MB6301/practical2\_words\_dbs\_v2.html <http://www.ucc.ie/microbio/MB6301/practical2_words_dbs_v2.html>`_)
you learnt that it is possible to create your own functions in R.
For example, if you want to create a function to find the square of
a number and add twenty to that, you can type:

::

    > myfunction <- function(x) { return((x*x) + 20) }

You can then use this function to find the value of 20 added to the
square of 100, by typing:

::

    > myfunction(100)
    [1] 10020

If you are likely to be using the same function again and again,
you can type it into a text file, for example, using Notepad.
(Note: Notepad is a program for writing text that is available on
all Windows computers. Go to the "Start" menu at the bottom left of
your Windows screen, and select "Notepad" from the list of programs
to start up Notepad.)

For example, say you use Notepad to type this line into a text
file, which you save as file "myfunctions.R" (to save a file in
Notepad, go to the "File" menu in Notepad and select "Save as", and
choose "All Files" from the "Save as type" drop-down list in the
"Save as" window that appears):

myfunction <- function(x) { return((x\*x) + 20) }

Whenever you start R next (for example, next week sometime), you
can then load the functions in file "myfunctions.R" into R by
typing:

::

    > source("myfunctions.R")

and then you can use them straight away. This is handy as it saves
you having to type your function into R each time that you start up
R again.

Note that when you save the file in Notepad, you must select "All
Files" from the "Save as type" drop-down list in the "Save as"
window. (If you don't do this, Notepad will automatically add
".txt" to the end of the file name that you specify, so for
example, your file will be renamed "myfunctions.R.txt").

Another very useful R function is the scan() function, which is for
reading data files into R. As its arguments (input), the scan()
function expects you to give it the name of the input file, and the
type of the data in that file, which is specified with the "what"
argument. The type of data can be numbers, in which case you need
to say "what=numeric", or can be characters/letters, in which case
you need to say "what=character". For example, if you have a data
file containing a column of numbers (for example, see the file
"Data1.txt" at
`http://www.ucc.ie/microbio/MB6301/Data1.txt <http://www.ucc.ie/microbio/MB6301/Data1.txt>`_),
you can read it into R by typing:

::

    > mydata <- scan("Data1.txt", what="numeric")
    > mydata # Print out the vector "mydata"
    [1] "33"   "22.9" "2233" "222"  "123"  "662" 

The command above reads in the column of numbers in the input file
"Data1.txt" into a vector *mydata*, so that each of the numbers is
stored in one element of the *mydata* vector. Similarly, you can
use scan() to read in an input data file in which the input data
consists of characters/letter (for example, see the file
"Data2.txt" at
`http://www.ucc.ie/microbio/MB6301/Data2.txt <http://www.ucc.ie/microbio/MB6301/Data2.txt>`_):

::

    > mydata2 <- scan("Data2.txt", what="character")
    > mydata2 # Print out the vector "mydata2"
    [1] "AFGAGADD" "DSDSDSDD" "CDSDGSDC"

Note that the scan() function expects that the input data file will
be in the "My Documents" folder on your computer, so you need to
copy your input file there so that scan() will be able to find it.

UniProt
-------

In previous practicals you learnt how to retrieve DNA and protein
sequences from the NCBI database. The NCBI database is a key
database in bioinformatics because it contains essentially all DNA
sequences ever sequenced.

As mentioned in previous practicals, a subsection of the NCBI
database called "RefSeq" consists of high quality DNA and protein
sequence data. Furthermore, the NCBI entries for the RefSeq
sequences have been *manually curated*, which means that biologists
employed by NCBI have added additional information to the NCBI
entries for those sequences, such as details of scientific papers
that describe the sequences.

Another extremely important manually curated database is UniProt
(`www.uniprot.org <http://www.uniprot.org>`_), which focuses on
protein sequences. UniProt aims to contains manually curated
information on all known protein sequences. While many of the
protein sequences in UniProt are also present in RefSeq, the amount
and quality of manually curated information in UniProt is much
higher than that in RefSeq.

For each protein in UniProt, the UniProt curators read all the
scientific papers that they can find about that protein, and add
information from those papers to the protein's UniProt entry. For
example, for a human protein, the UniProt entry for the protein
usually includes information about the biological function of the
protein, in what human tissues it is expressed, whether it
interacts with other human proteins, and much more. All this
information has been manually gathered by the UniProt curators from
scientific papers, and the papers in which the found the
information are always listed in the UniProt entry for the
protein.

Just like NCBI, UniProt also assigns an *accession* to each
sequence in the UniProt database. Although the same protein
sequence may appear in both the NCBI database and the UniProt
database, it will have different NCBI and UniProt accessions.
However, there is usually a link on the NCBI entry for the protein
sequence to the UniProt entry, and vice versa.

Viewing the UniProt page for a protein sequence
-----------------------------------------------

If you are given the UniProt accession for a protein, to find the
UniProt entry for the protein, you first need to go the UniProt
website, `www.uniprot.org <http://www.uniprot.org>`_. At the top of
the UniProt website, you will see a search box, and you can type
the accession of the protein that you are looking for in this
search box, and then click on the "Search" button to search for
it.

For example, if you want to find the sequence for the chorismate
lyase protein from *Escherichia coli* strain K12, which has UniProt
accession P26602, you would type just "P26602" in the search box
and press "Search":

|image0|

The UniProt entry for UniProt accession P26602 will then appear in
your web browser. The picture below shows the top part of the
UniProt entry for accession P26602. You can see there is a lot of
information about the protein in its UniProt entry.

Beside the heading "Organism" you can see the organism is given as
*Escherichia coli*(strain K12). Beside the heading "Taxonomic
lineage", you can see "Bacteria > Proteobacteria >
Gammaproteobacteria > Enterobacteriales > Enterobacteriaceae >
Escherichia". This tells us that *Escherichia* is a species of
bacteria, which belongs to a group of related bacteria called the
Enterobacteriaceae, which itself belongs to a larger group of
related bacteria called the Enterobacteriales, which itself belongs
to an even larger group of related bacteria called the
Enterobacteriales, which itself belongs to the Gammaproteobacteria,
which itself belongs to a huge group of bacteria called the
Proteobacteria.

Beside the heading "Sequence length" we see that the sequence is
165 amino acids long (165 letters long). Further down, beside the
heading "Function", it says that the function of this protein is
that it "Removes the pyruvyl group from chorismate, with
concomitant aromatization of the ring, to provide 4-hydroxybenzoate
(4HB) for the ubiquinone pathway". This tells us this protein is an
enzyme (a protein that increases the rate of a specific biochemical
reaction), and tells us what is the particular biochemical reaction
that this enzyme is involved in.

Further down the UniProt page for this protein, you will see a lot
more information, as well as many links to webpages in other
biological databases, such as NCBI. The huge amount of information
about proteins in UniProt means that if you want to find out about
a particular protein, the UniProt page for that protein is a great
place to start.

|image1|

Retrieving a protein sequence from UniProt as a FASTA-format file
-----------------------------------------------------------------

To retrieve a FASTA-format file containing the sequence for a
particular protein, you need to look at the top right of the
UniProt entry for the protein on the UniProt website. You will see
a small orange button labelled "FASTA", which you should click on:

|image2|

The FASTA-format sequence for the accession will now appear in your
web browser. To save it as a file, go to the "File" menu of your
web browser, choose "Save page as", and save the file. Remember to
give the file a sensible name (eg. "P26602.fasta" for accession
P26602), and in a place that you will remember (eg. in the "My
Documents" folder).

For example, you can retrieve the protein sequences for the
chorismate lyase protein from *Escherichia coli* strain K12 (which
has UniProt accession P26602) and for the chorismate lyase protein
from *Salmonella typhi* (UniProt accession Q8Z1T7), and save them
as FASTA-format files (eg. "P26602.fasta" and "Q8Z1T7.fasta", as
described above.

Note that the *Escherichia coli* and *Salmonella typhi* chorismate
lyase proteins are an example of a pair of homologous (related)
proteins in two related species of bacteria. *Escherichia coli* is
often studied in microbiology laboratories as a model organism, and
is also one of the foremost causes of food poisoning.
*Salmonella typhi* is relatively closely related bacterium, and is
the bacterium that causes typhoid fever, a common illness
worldwide.

Once you have downloaded the protein sequences for UniProt
accessions P26602 and Q8Z1T7 and saved them as FASTA-format files
(eg. "P26602.fasta" and "Q8Z1T7.fasta"), you can read them into R
using the read.fasta() function in the SeqinR R library (as
described in Practical 1,
`www.ucc.ie/microbio/MB6301/practical1\_words\_v2.html#NCBI <http://www.ucc.ie/microbio/MB6301/practical1_words_v2.html#NCBI>`_).
Remember that the read.fasta() function expects that you have put
your FASTA-format files in the "My Documents" folder on your
computer.

For example, the following commands will read the FASTA-format
files P26602.fasta and Q8Z1T7.fasta into R, and store the two
protein sequences in two vectors *coliseq* and *typhiseq*:

::

    > library("seqinr")
    > coli <- read.fasta(file = "P26602.fasta")
    > typhi <- read.fasta(file = "Q8Z1T7.fasta")
    > coliseq <- coli[[1]]
    > typhiseq <- typhi[[1]]
    > coliseq # Display the contents of the vector "coliseq"
      [1] "m" "s" "h" "p" "a" "l" "t" "q" "l" "r" "a" "l" "r" "y" "c" "k" "e" "i"
     [19] "p" "a" "l" "d" "p" "q" "l" "l" "d" "w" "l" "l" "l" "e" "d" "s" "m" "t"
     [37] "k" "r" "f" "e" "q" "q" "g" "k" "t" "v" "s" "v" "t" "m" "i" "r" "e" "g"
     [55] "f" "v" "e" "q" "n" "e" "i" "p" "e" "e" "l" "p" "l" "l" "p" "k" "e" "s"
     [73] "r" "y" "w" "l" "r" "e" "i" "l" "l" "c" "a" "d" "g" "e" "p" "w" "l" "a"
     [91] "g" "r" "t" "v" "v" "p" "v" "s" "t" "l" "s" "g" "p" "e" "l" "a" "l" "q"
    [109] "k" "l" "g" "k" "t" "p" "l" "g" "r" "y" "l" "f" "t" "s" "s" "t" "l" "t"
    [127] "r" "d" "f" "i" "e" "i" "g" "r" "d" "a" "g" "l" "w" "g" "r" "r" "s" "r"
    [145] "l" "r" "l" "s" "g" "k" "p" "l" "l" "l" "t" "e" "l" "f" "l" "p" "a" "s"
    [163] "p" "l" "y"

Pairwise global alignment of DNA sequences using the Needleman-Wunsch algorithm
-------------------------------------------------------------------------------

If you are studying a particular pair of genes or proteins, an
important question is to what extent the two sequences are similar.
To quantify similarity, it is necessary to *align* the two
sequences, and then you can calculate a similarity score based on
the alignment.

There are two types of alignment in general. A *global* alignment
is an alignment of the full length of two sequences, for example,
of two protein sequences or of two DNA sequences. A *local*
alignment is an alignment of part of one sequence to part of
another sequence.

The first step in computing a alignment (global or local) is to
decide on a scoring system. For example, we may decide to give a
score of +2 to a match and a penalty of -1 to a mismatch, and a
penalty of -2 to a gap. Thus, for the alignment:

::

    G A A T T C
    G A T T - A

we would compute a score of 2 + 2 -1 + 2 -2 - 1 = 2. Similarly, the
score for the following alignment is 2 + 2 -2 + 2 + 2 -1 = 5:

::

    G A A T T C
    G A - T T A

The scoring system above can be represented by a *scoring matrix*,
*σ* (also known as a *substitution matrix*). The matrix *σ*
(pronounced "sigma") has one row and one column for each possible
letter in our alphabet of letters (eg. 4 rows and 4 columns for DNA
sequences). The *(i,j)* element of *σ*, *σ(i,j)* has a value of +2
in case of a match and -1 in case of a mismatch.

We can make a scoring matrix in R by using the
nucleotideSubstitutionMatrix() function in the Biostrings()
library. The Biostrings library is part of a set of R libraries for
bioinformatics analysis known as Bioconductor
(`www.bioconductor.org/ <http://www.bioconductor.org/>`_). The
arguments (inputs) for the nucleotideSubstitutionMatrix() function
are the score that we want to assign to a match and the score that
we want to assign to a mismatch. We can also specify that we want
to use only the four letters representing the four nucleotides (ie.
A, C, G, T) by setting 'baseOnly=TRUE', or whether we also want to
use the letters that represent ambiguous cases where we are not
sure what the nucleotide is (eg. 'N' = A/C/G/T).

To make a scoring matrix which assigns a score of +2 to a match and
-1 to a mismatch, we type:

::

    > library(Biostrings)
    > sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
    > sigma # Print out the matrix
       A  C  G  T
    A  2 -1 -1 -1
    C -1  2 -1 -1
    G -1 -1  2 -1
    T -1 -1 -1  2

Instead of assigning the same penalty (eg. -8) to every gap
position, it is common instead to assign a *gap opening penalty* to
the first position in a gap (eg. -8), and a smaller
*gap extension penalty* to every subsequent position in the same
gap. The reason for doing this is that it is likely that adjacent
gap positions were created by the same insertion or deletion event,
rather than by several independent insertion or deletion events.
Therefore, we don't want to penalise a 3-letter gap as much as we
would penalise three separate 1-letter gaps, as the 3-letter gap
may have arisen due to just one insertion or deletion event, while
the 3 separate 1-letter gaps probably arose due to three
independent insertion or deletion events.

For example, if we want to compute the score for a global alignment
of two short DNA sequences 'GAATTC' and 'GATTA', we can use the
Needleman-Wunsch algorithm to calculate the highest-scoring
alignment using a particular scoring function. The
pairwiseAlignment() function in the Biostrings R library finds the
score for the optimal global alignment between two sequences using
the Needleman-Wunsch algorithm, given a particular scoring system.
As arguments (inputs), the pairwiseAlignment() function takes the
two sequences that you want to align, the scoring matrix, the gap
opening penalty, and the gap extension penalty. You can also tell
the function that you want to just have the optimal global
alignment's score by setting "scoreOnly = TRUE", or that you want
to have both the optimal global alignment and its score by setting
"scoreOnly = FALSE". For example, to find the score for the optimal
global alignment between the sequences 'GAATTC' and 'GATTA', we
type:

::

    > s1 <- "GAATTC"
    > s2 <- "GATTA"
    > globalAligns1s2 <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -2, 
    gapExtension = -8, scoreOnly = FALSE)
    > globalAligns1s2 # Print out the optimal alignment and its score
    Global Pairwise Alignment (1 of 1)
    pattern: [1] GAATTC 
    subject: [1] GA-TTA 
    score: -3

The above commands print out the optimal global alignment for the
two sequences and its score.

Pairwise global alignment of protein sequences using the Needleman-Wunsch algorithm
-----------------------------------------------------------------------------------

As well as DNA alignments, it is also possible to make alignments
of protein sequences. In this case it is necessary to use a scoring
matrix *σ* for amino acids rather than for nucleotides. There are
several well known scoring matrices that come with R, such a the
BLOSUM series of matrices. Different BLOSUM matrices exist, named
with different numbers. BLOSUM with high numbers are designed for
comparing closely related sequences, while BLOSUM with low numbers
are designed for comparing distantly related sequences. For
example, BLOSUM62 is used for less divergent alignments (alignments
of sequences that differ little), and BLOSUM30 is used for more
divergent alignments (alignments of sequences that differ a lot).

Many R libraries come with example data sets or data files. The
data() function is used to load these data files. You can use the
data() function in R to load a data set of BLOSUM matrices that
comes with R Biostrings() library. To load the BLOSUM50 matrix, we
type:

::

    > data(BLOSUM50)
    > BLOSUM50 # Print out the data
       A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -5
    R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1 -5
    N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1 -5
    D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1 -5
    C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2 -5
    Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1 -5
    E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1 -5
    G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2 -5
    H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1 -5
    I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1 -5
    L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1 -5
    K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1 -5
    M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1 -5
    F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2 -5
    P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2 -5
    S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1 -5
    T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0 -5
    W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3 -5
    Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1 -5
    V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1 -5
    B -2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  2 -1 -5
    Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  2  5 -1 -5
    X -1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1 -5
    * -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1

You can get a list of the available scoring matrices that come with
the Biostrings library by using the data() function, which takes as
an argument the name of the library for which you want to know the
data sets that come with it:

::

    > data(package="Biostrings")
    Data sets in package 'Biostring':
    BLOSUM100                                 Scoring matrices
    BLOSUM45                                  Scoring matrices
    BLOSUM50                                  Scoring matrices
    BLOSUM62                                  Scoring matrices
    BLOSUM80                                  Scoring matrices

To find the optimal global alignment between the protein sequences
"PAWHEAE" and "HEAGAWGHEE" using the Needleman-Wunsch algorithm
using the BLOSUM50 matrix, we type:

::

    > data(BLOSUM50)
    > s3 <- "PAWHEAE"
    > s4 <- "HEAGAWGHEE"
    > globalAligns3s4 <- pairwiseAlignment(s3, s4, substitutionMatrix = "BLOSUM50", gapOpening = -2,
    gapExtension = -8, scoreOnly = FALSE)
    > globalAligns3s4 # Print out the optimal global alignment and its score
    Global Pairwise Alignment (1 of 1)
    pattern: [1] P---AWHEAE 
    subject: [1] HEAGAWGHEE 
    score: -5 

Retreiving sequences from a database, and aligning them
-------------------------------------------------------

In previous practicals and this one, you learnt how to retrieve
sequences from sequence databases such as NCBI and UniProt, to save
them as FASTA-format files, and then to read them into R using the
read.fasta() function.

For example, earlier in this practical, you learnt how to retrieve
the sequences for the chorismate lyase proteins from
*Escherichia coli* strain K12 (UniProt P26602) and
*Salmonella typhi* (UniProt Q8Z1T7), and read them into R, and
store them as vectors *coliseq* and *vectorseq*.

Once you have read in sequences using read.fasta(), you can align
them using pairwiseAlignment() from the Biostrings library.

As its input, the pairwiseAlignment() function requires that the
sequences be in the form of a single string (eg. "ACGTA"), rather
than as a vector of characters (eg. a vector with the first element
as "A", the second element as "C", etc.). Therefore, to align the
*E. coli* and *S. typhi* chorismate lyase proteins, we first need
to convert the vectors *coliseq* and *vectorseq* into strings. We
can do this using the c2s() function in the SeqinR library:

::

    > coliseq # Print the content of vector "coliseq"
      [1] "m" "s" "h" "p" "a" "l" "t" "q" "l" "r" "a" "l" "r" "y" "c" "k" "e" "i"
     [19] "p" "a" "l" "d" "p" "q" "l" "l" "d" "w" "l" "l" "l" "e" "d" "s" "m" "t"
     [37] "k" "r" "f" "e" "q" "q" "g" "k" "t" "v" "s" "v" "t" "m" "i" "r" "e" "g"
     [55] "f" "v" "e" "q" "n" "e" "i" "p" "e" "e" "l" "p" "l" "l" "p" "k" "e" "s"
     [73] "r" "y" "w" "l" "r" "e" "i" "l" "l" "c" "a" "d" "g" "e" "p" "w" "l" "a"
     [91] "g" "r" "t" "v" "v" "p" "v" "s" "t" "l" "s" "g" "p" "e" "l" "a" "l" "q"
    [109] "k" "l" "g" "k" "t" "p" "l" "g" "r" "y" "l" "f" "t" "s" "s" "t" "l" "t"
    [127] "r" "d" "f" "i" "e" "i" "g" "r" "d" "a" "g" "l" "w" "g" "r" "r" "s" "r"
    [145] "l" "r" "l" "s" "g" "k" "p" "l" "l" "l" "t" "e" "l" "f" "l" "p" "a" "s"
    [163] "p" "l" "y"
    > coliseqstring <- c2s(coliseq) # Make a string that contains the sequence in "coliseq"
    > coliseqstring # Print the content of string coliseqstring
    [1] "mshpaltqlralryckeipaldpqlldwllledsmtkrfeqqgktvsvtmiregfveqneipeelpllpkesrywlreillcadgepwlagrtvvpvstlsgpelalqklgktplgrylftsstltrdfieigrdaglwgrrsrlrlsgkpllltelflpasply"
    > typhiseqstring <- c2s(typhiseq) # Make a string that contains the sequence in "typhiseq"

Furthermore, pairwiseAlignment() requires that the sequences be
stored as uppercase characters. Therefore, we need to use the
toupper() function to convert *coliseqstring* and *typhiseqstring*
to uppercase:

::

    > coliseqstring <- toupper(coliseqstring)
    > coliseqstring # Print out the content of vector "coliseqstring"
    [1] "MSHPALTQLRALRYCKEIPALDPQLLDWLLLEDSMTKRFEQQGKTVSVTMIREGFVEQNEIPEELPLLPKESRYWLREILLCADGEPWLAGRTVVPVSTLSGPELALQKLGKTPLGRYLFTSSTLTRDFIEIGRDAGLWGRRSRLRLSGKPLLLTELFLPASPLY"
    > typhiseqstring <- toupper(typhiseqstring) 

We can now align the the *E. coli* and *S. tytphi* chorismate lyase
protein sequences using the pairwiseAlignment() function:

::

    > globalAlignColiTyphi <- pairwiseAlignment(coliseqstring, typhiseqstring,
    substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
    > globalAlignColiTyphi # Print out the optimal global alignment and its score
    Global Pairwise Alignment (1 of 1)
    pattern: [1] MSHPALTQLRALRYCKEIPALDPQLLDWLLLEDSMTKRFEQQGKTVSVTMI...GRYLFTSSTLTRDFIEIGRDAGLWGRRSRLRLSGKPLLLTELFLPASPLY 
    subject: [1] MSHPALTQLRALRYFDAIPALEPHLLDWLLLEDSVTKRFEQQGKRVSVTLI...GRYLFTSSTLTRDFIEIGRDATLWGRRSRLRLSGKPLLLTELFLPASPLY 
    score: 931 

As the alignment is very long, when you type
*globalAlignColiTyphi*, you only see the start and the end of the
alignment (see above). Therefore, we need to have a function to
print out the whole alignment (see below).

Viewing a long pairwise alignment
---------------------------------

If you want to view a long pairwise alignment such as that between
the *E. coli* and *S. typhi* chorismate lyase proteins, it is
convenient to print out the alignment in blocks. As mentioned
above, it is possible to save R functions in a file, and to load
them for us at a later date. The file Rfunctions.R (which you can
download from the web at
`www.ucc.ie/microbio/MB6301/Rfunctions.R <http://www.ucc.ie/microbio/MB6301/Rfunctions.R>`_)
contains a function printPairwiseAlignment() that contains a
function for printing out blocks of an alignment. Download the file
to your "My Documents" folder, and load it into R by using the
source() function:

::

    > source("Rfunctions.R")

We can then use our function printPairwiseAlignment() to print out
the alignment between the *E. coli* and *S. typhi* chorismate lyase
proteins (which we stored this alignment in the
*globalAlignColiTyphi* variable, see above), in blocks of 60
alignment columns:

::

    > printPairwiseAlignment(globalAlignColiTyphi, 60)
    [1] "MSHPALTQLRALRYCKEIPALDPQLLDWLLLEDSMTKRFEQQGKTVSVTMIREGFVEQNE 60"
    [1] "MSHPALTQLRALRYFDAIPALEPHLLDWLLLEDSVTKRFEQQGKRVSVTLIREAFVGQSE 60"
    [1] " "
    [1] "IPEELPLLPKESRYWLREILLCADGEPWLAGRTVVPVSTLSGPELALQKLGKTPLGRYLF 120"
    [1] "VEEASGLLPSESRYWLREILLCADGEPWLAGRTVVPESTLCGPEQVLQHLGKTPLGRYLF 120"
    [1] " "
    [1] "TSSTLTRDFIEIGRDAGLWGRRSRLRLSGKPLLLTELFLPASPLY 180"
    [1] "TSSTLTRDFIEIGRDATLWGRRSRLRLSGKPLLLTELFLPASPLY 180"
    [1] " "

The position in the protein of the amino acid that is at the end of
each line of the printed alignment is shown after the end of the
line. For example, the first line of the alignment above finishes
at amino acid position 60 in the *E. coli* protein and also at
amino acid position 60 in the *S. typhi* protein.

If we were printing out an alignment that contained gaps in the
first 60 alignment columns, the first 60 alignment columns may end
before the 60th amino acid in the two sequences that were aligned.

Calculating the statistical significance of a pairwise global alignment
-----------------------------------------------------------------------

We have seen that when we align the 'PAWHEAE' and 'HEAGAWGHEE'
protein sequences, they have some similarity, and the score for
their optimal global alignment is -5. But is this alignment
*statistically significant*? In other words, is this alignment
better than we would expect between any two random proteins? The
Needleman-Wunsch alignment algorithm will produce a global
alignment even if we give it two unrelated random protein
sequences, although the alignment score would be low. Therefore, we
should ask: is the score for our alignment better than expected
between two random sequences of the same lengths and amino acid
compositions?

It is reasonable to expect that if the alignment score is
statistically significant, then it will be higher than the scores
obtained from aligning pairs of random protein sequences that have
the same lengths and amino acid compositions as our original two
sequences. Therefore, to assess if the score for our alignment
between the 'PAWHEAE' and 'HEAGAWGHEE' protein sequence is
statistically significant, a first step is to make some random
sequences that have the same amino acid composition and length as
one of our initial two sequences, for example, as the same amino
acid composition and length as the sequence 'PAWHEAE'.

How can we obtain random sequences of the same amino acid
composition and length as the sequence 'PAWHEAE'? One way is to
generate sequences using a
*multinomial model for protein sequences* in which the
probabilities of the different amino acids set to be equal to their
frequencies in the sequence 'PAWHEAE'. That is, we can generate
sequences using a multinomial model for proteins, in which the
probability of 'P', *p\ :sub:`P`\ *, is set to 0.1428571 (1/7); the
probability of 'A', *p\ :sub:`A`\ *, is set to 0.2857143 (2/7); the
probability of 'W', *p\ :sub:`W`\ *, is set to 0.1428571 (1/7); the
probability of 'H', *p\ :sub:`H`\ *, is set to 0.1428571 (1/7); and
the probabilty of 'E', *p\ :sub:`E`\ *, is set to 0.2857143 (2/7),
and the probabilities of the other 15 amino acids are set to 0.

To generate a sequence with this multinomial model, we choose the
letter for each position in the sequence according to those
probabilities. This is as if we have made a roulette wheel in which
1/7*th* of the circle is taken up by a pie labelled "P", 2/7*ths*
by a pie labelled "A", 1/7*th* by a pie labelled "W", 1/7*th* by a
pie labelled "H", and 2/7*ths* by a pie labelled "E":

|image3|

To generate a sequence using the multinomial model, we keep
spinning the arrow in the centre of the roulette wheel, and write
down the letter that the arrow stops on after each spin. To
generate a sequence that is 7 letters long, we can spin the arrow 7
times. To generate 1000 sequences that are each 7 letters long, we
can spin the arrow 7000 times, where the letters chosen form 1000
7-letter amino acid sequences.

The procedure above was used to generate 1000 7-letter amino acid
sequences, using a multinomial model in which
*p\ :sub:`P`\ *=0.1428571, *p\ :sub:`A`\ *=0.2857143,
*p\ :sub:`W`\ *=0.1428571, *p\ :sub:`H`\ *=0.1428571, and
*p\ :sub:`E`\ *=0.2857143. The file containing the sequences is
called "SeqsFromMultinomial1", and you can download it from
`http://www.ucc.ie/microbio/MB6301/SeqsFromMultinomial1 <http://www.ucc.ie/microbio/MB6301/SeqsFromMultinomial1>`_.
(Note: I have not shown you here how to use R to generate a
sequence according to a particular multinomial model, because this
is one of the problems you need to solve for MB6301 Assignment 1).

You can then read the sequences in file "SeqsFromMultinomial1" into
R using the scan() function:

::

    > randomseqs1 <- scan("SeqsFromMultinomial1", what="character")

This reads the sequences into a vector *randomseqs1*, where each
element in the vector contains one of the sequences generated using
the multinomial model. Therefore, to print out the first ten
sequences, we can print out the first ten elements of the vector
*randomseqs1*:

::

    > randomseqs1[1:10]
     [1] "EEHAAAE" "AWPHPHA" "AEAPHWE" "WAHAAHA" "HAEEPHP" "APPWAWA" "WEPPPPH"
     [8] "HEAEHWA" "EEEHWPP" "WWAAEAW"

We can then use the Needleman-Wunsch algorithm to align the
sequence 'HEAGAWGHEE' to one of the 1000 random sequences generated
using the multinomial model with *p\ :sub:`P`\ *=0.1428571,
*p\ :sub:`A`\ *=0.2857143, *p\ :sub:`W`\ *=0.1428571,
*p\ :sub:`H`\ *=0.1428571, and *p\ :sub:`E`\ *=0.2857143. For
example, to align 'HEAGAWGHEE' to the first of the 1000 random
sequences ('EEHAAAE'), we type:

::

    > s4 <- "HEAGAWGHEE"
    > pairwiseAlignment(s4, randomseqs1[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
    gapExtension = -8, scoreOnly = FALSE)
    Global Pairwise Alignment (1 of 1)
    pattern: [1] HEAGAWGHEE 
    subject: [1] EEHAA---AE 
    score: -12 

If we use the pairwiseAlignment() function with the argument
'scoreOnly=TRUE', it will just give us the score for the
alignment:

::

    > pairwiseAlignment(s4, randomseqs1[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
    gapExtension = -8, scoreOnly = TRUE)
    [1] -12

If we repeat this 1000 times, that is, for each of the 1000 random
sequences in vector *randomseqs1*, we can get a distribution of
alignment scores expected for aligning 'HEAGAWGHEE' to random
sequences of the same length and (approximately the same) amino
acid composition as 'PAWHEAE'. We can then compare the actual score
for aligning 'PAWHEAE' to 'HEAGAWGHEE' (ie. -5) to the distribution
of scores for aligning 'HEAGAWGHEE' to the random sequences.

::

    > randomscores <- double(1000) # Create a numeric vector with 1000 elements
    > for (i in 1:1000) {
            score <- pairwiseAlignment(s4, randomseqs1[i], substitutionMatrix = "BLOSUM50", gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
            randomscores[i] <- score
     }

The code above first uses the double() function to create a numeric
vector *randomscores* for storing real numbers (ie. not integers),
with 1000 elements. This will be used to store the alignment scores
for 1000 alignments between 'HEAGAWGHEE' and the 1000 different
random sequences generated using the multinomial model. The 'for
loop' takes each of the 1000 different random sequences, aligns
each one to 'HEAGAWGHEE', and stores the 1000 alignment scores in
the *randomscores* vector. We can make a histogram plot of the 1000
scores in vector *randomscores* by typing:

::

    > hist(randomscores, col="red") # Draw a red histogram

|image4| We can see from the histogram that a lot of the random
sequences seem to have higher alignment scores than -3 when aligned
to 'HEAGAWGHEE' (where -3 is the alignment score for 'PAWHEAE' and
'HEAGAWGHEE').

We can use the vector *randomscores* of scores for 1000 alignments
of random sequences to 'HEAGAWGHEE' to calculate the probability of
getting a score as large as the real alignment score for 'PAWHEAE'
and 'HEAGAWGHEE' (ie. -5) by chance.

::

    > sum(randomscores >= -5)
    [1] 289

We see that 289 of the 1000 alignments of random sequences to
'HEAGAWGHEE' had alignment scores that were equal to or greater
than -5. Thus, we can estimate that the probability of getting a
score as large as the real alignment score by chance is (289/1000
=) 0.289. In other words, we can calculate a *P-value* of 0.289.
This probability or *P*-value is quite high (almost 30%, or 1 in
3), so we can conclude that it is quite probable that we could get
an alignment score as high as -5 by chance alone. This indicates
that the sequences 'HEAGAWGHEE' and 'PAWHEAE' are not more similar
than any two random sequences, and so they are probably not related
sequences.

.. include:: <isoamsr.txt>

Another way of saying this is that the *P*-value that we calculated
is high (0.289), and as a result we conclude that the alignment
score for the sequences 'HEAGAWGHEE' and 'PAWHEAE' is not
*statistically significant*. Generally, if the *P*-value that we
calculate for an alignment of two sequences is >0.05, we conclude
that the alignment score is not statistically significant, and that
the sequences are probably not related. On the other hand, if the
*P*-value is |les| 0.05, we conclude that the alignment score is
statistically significant, and the sequences are very probably
related (homologous).

Summary
-------

In this practical, you will have learnt to use the following R
functions:


#. data() for reading in data that comes with an R library
#. double() for creating a numeric vector for storing real
   (non-integer) numbers
#. toupper() for converting a string of characters from lowercase
   to uppercase
#. scan() for reading in data from a file

All of these functions belong to the standard installation of R.

You have also learnt the following R functions that belong to the
bioinformatics libraries:


#. nucleotideSubstitutionMatrix() in the Biostrings library for
   making a nucleotide scoring matrix
#. pairwiseAlignment() in the Biostrings library for making a
   global alignment between two sequences
#. c2s() in the SeqinR library for converting a sequence stored in
   a vector to a string of characters

Links and Further Reading
-------------------------

Some links are included here for further reading, which will be
especially useful if you need to use the R package and SeqinR
library for your project or assignments.

For background reading on sequence alignment, it is recommended to
read Chapter 3 of
*Introduction to Computational Genomics: a case studies approach*
by Cristianini and Hahn (Cambridge University Press;
`www.computational-genomics.net/book/ <http://www.computational-genomics.net/book/>`_).

For more in-depth information and more examples on using the SeqinR
library for sequence analysis, look at the SeqinR documentation,
`seqinr.r-forge.r-project.org/seqinr\_2\_0-1.pdf <http://seqinr.r-forge.r-project.org/seqinr_2_0-1.pdf>`_.

For more information on and examples using the Biostrings library,
see the Biostrings documentation at
`bioconductor.org/packages/2.5/bioc/html/Biostrings.html <http://bioconductor.org/packages/2.5/bioc/html/Biostrings.html>`_.

There is also a very nice chapter on "Analyzing Sequences", which
includes examples of using the SeqinR and Biostrings libraries for
sequence analysis, in the book
*Applied statistics for bioinformatics using R* by Krijnen
(available online at
`cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf <http://cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf>`_).

Acknowledgements
----------------

Many of the ideas for the examples and exercises for this practical
were inspired by the Matlab case study on the Eyeless protein
(`www.computational-genomics.net/case\_studies/eyeless\_demo.html <http://www.computational-genomics.net/case_studies/eyeless_demo.html>`_)
from the website that accompanies the book
*Introduction to Computational Genomics: a case studies approach*
by Cristianini and Hahn (Cambridge University Press;
`www.computational-genomics.net/book/ <http://www.computational-genomics.net/book/>`_).

The examples of DNA sequences and protein sequences to align
('GAATTC' and 'GATTA', and sequences 'PAWHEAE' and 'HEAGAWGHEE'),
as well as some ideas related to finding the statistical
significance of a pairwise alignment, were inspired by the chapter
on "Analyzing Sequences" in the book
*Applied statistics for bioinformatics using R* by Krijnen
(`cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf <http://cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf>`_).

Thank you to Jean Lobry and Simon Penel for helpful advice on using
the SeqinR library.

Exercises
---------

Answer the following questions, using the R package. For each
question, please record your answer, and what you typed into R to
get this answer.

Q1. Download FASTA-format files of the *Drosophila melanogaster* Eyeless protein (UniProt accession O18381) and the human Aniridia protein (UniProt accession Q66SS1) sequences from UniProt (`www.uniprot.org <http://www.uniprot.org>`_). 
    Note: the *eyeless* gene of the fruitfly *Drosophila melanogaster*
    and the human gene *aniridia* are distantly related genes that
    control eye development in these two species. Some regions of the
    fruitfly *eyeless* gene and human *aniridia* gene are almost
    identical.
Q2. What is the alignment score for the optimal global alignment between the *Drosophila melanogaster* Eyeless protein and the human Aniridia protein, when you use the *BLOSUM50* scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5? Q3. Use the printPairwiseAlignment() function to view the optimal global alignment between *Drosophila melanogaster* Eyeless protein and the human Aniridia protein, using the *BLOSUM50* scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5. 
    Do you see any regions where the alignment is very good (lots of
    identities and few gaps)?
    If so, what are the coordinates of these long regions of good
    alignment with respect to the Aniridia and Eyeless proteins? It is
    sufficient to give approximate coordinates.
Q4. What global alignment score do you get for the Aniridia and Eyeless proteins, when you use the *BLOSUM62* alignment matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5? 
    Which scoring matrix do you think is more appropriate for using for
    this pair of proteins: BLOSUM50 or BLOSUM62?
Q5. What is the statistical significance of the optimal global alignment for the Aniridia and Eyeless proteins made using the *BLOSUM50* scoring matrix, with a gap opening penalty of -10 and a gap extension penalty of -0.5? 
    In other words, what is the probability of getting a score as large
    as the real alignment score for Eyeless and Aniridia by chance?
    Hint: to answer your question, you may find it useful to analyse
    this file containing 1000 random sequences generated using a
    multinomial model with the probabilities of the 20 amino acids set
    equal to their frequencies in the *D. melanogaster* Eyeless protein
    sequence:
    `http://www.ucc.ie/microbio/MB6301/SeqsFromMultinomial2 <http://www.ucc.ie/microbio/MB6301/SeqsFromMultinomial2>`_.
    Each of these 1000 random sequences is the same length as Eyeless.
Q6. What is the optimal global alignment score between the *Drosophila melanogaster* Eyeless protein and the *E. coli* chorismate lyase protein? 
    Is the alignment score statistically significant (what is the
    *P-*value?)?
    Does this surprise you?

Other ways to do the same thing
-------------------------------

It is possible to carry out some of the analyses that you have
carried out in the practicals via websites. For example, it is
possible to calculate the optimal global alignment between two
sequences using the Needleman-Wunsch algorithm using the Needle
program, via the website
`mobyle.rpbs.univ-paris-diderot.fr/cgi-bin/portal.py?form=needle <http://mobyle.rpbs.univ-paris-diderot.fr/cgi-bin/portal.py?form=needle>`_.
The Needle is also available to download as part of the EMBOSS
package
(`emboss.sourceforge.net <http://emboss.sourceforge.net/>`_), and
so can also be run on your own computer.




.. |image0| image:: ../_static/P4_image7.png
.. |image1| image:: ../_static/P4_image8.png
.. |image2| image:: ../_static/P4_image9.png
.. |image3| image:: ../_static/P4_image10.png
.. |image4| image:: ../_static/P4_image11.png
