Practical 5 for 2009/2010 - Phylogenetic trees
==============================================



A little more about R
---------------------

In previous practicals (see practical 2, at
`http://www.ucc.ie/microbio/MB6301/practical2\_words\_dbs\_v2.html <http://www.ucc.ie/microbio/MB6301/practical2_words_dbs_v2.html>`_),
you learnt that you can define your own functions in R using the
function() command to define a function. We can define a function
that carries out several functions one after another. For example,
if we wanted to define a function that calculates the sum of the
elements in a vector, and then calculates the log to the base 10 of
that sum, we would type:

::

    > myfunction2 <- function(x) { return(log10(sum(x))) } 
    > myvector <- c(12, 133, 423, 333)
    > myfunction2(myvector)
    [1] 2.954725

Note that the function can be defined over several lines, and will
do exactly the same thing:

::

    > myfunction <- function(x)
         {
            mysum <- sum(x)        # Calculate the sum of numbers in vector x
            mylog <- log10(mysum)  # Calculate the log of mysum
            return(mylog)
         }
    > myvector <- c(12, 133, 423, 333)
    > myfunction(myvector)
    [1] 2.954725

If you are using a computer that is directly connected to the
internet, the procedure for installing most R libraries is simple.
You can use the function install.packages() for this, which takes
the name of the R library that you want to install as its argument
(input). For example, to install the Ape R library, which is a
library for phylogenetic analysis that you will use in this
practical, you can type:

::

    > install.packages("ape")

When you type the command above, you may be asked which website to
download the library from, as the R libraries are available for
download from various websites worldwide. It is a good idea to
select the Irish website (managed by the Irish Higher Education
Authority, or HEA), as this will be quickest to download from.

Note that the install.packages() command can be used to install a
large number of R libraries. However, there are some R libraries
that cannot be installed using this function, including the
libraries in the Bioconductor set of R libraries, which have a
special installation procedure. Furthermore, the install.packages()
command can only be used if you are using a computer that has a
direct internet connection (and so cannot be used in many of the
computer labs on the UCC campus, which connect to the internet
indirectly, via a proxy).

Retrieving sequences from a database using R
--------------------------------------------

In previous practicals (see Practical 3, at
`http://www.ucc.ie/microbio/MB6301/practical3\_revised\_v2.html <http://www.ucc.ie/microbio/MB6301/practical3_revised_v2.html>`_),
you learnt how to search for DNA or protein sequences in sequence
databases such as FASTA-format files. You can also search and
download the sequences using R, with the SeqinR R library. The file
"Rfunctions.R" (which can be downloaded from
`www.ucc.ie/microbio/MB6301/Rfunctions.R <http://www.ucc.ie/microbio/MB6301/Rfunctions.R>`_)
contains a function retrieveuniprotseqs() which uses the SeqinR R
library to search for and retrieve sequences from the UniProt
database. As its arguments (inputs), the retrieveuniprotseqs()
function takes a vector containing the names of the sequences. The
retrieveuniprotseqs() function returns a list variable, in which
each element is a vector containing one of the sequences.

For example, to retrieve the protein sequences for UniProt
accessions Q9NVS1 and P62286, you can type:

::

    > source("Rfunctions.R")
    > seqnames <- c("Q9NVS1", "P62286") # Make a vector containing the names of the sequences
    > seqs <- retrieveuniprotseqs(seqnames)    # Retrieve the sequences and store them in list variable "seqs"
    > length(seqs)                      # Print out the number of sequences retrieved
    [1] 2
    > seq1 <- seqs[[1]]                 # Get the first sequence
    > seq1[1:20]                        # Print out the first 20 letters of the first sequence            
     [1] "M" "A" "N" "R" "R" "V" "G" "R" "G" "C" "W" "E" "V" "S" "P" "T" "E" "R" "R" "P"
    > seq2 <- seqs[[2]]                 # Get the second sequence
    > seq2[1:20]                        # Print out the first 20 letters of the second sequence
     [1] "M" "A" "T" "R" "R" "A" "G" "R" "S" "W" "E" "V" "S" "P" "S" "G" "P" "R" "P" "A"

The commands above use the function retrieveuniprotseqs() to
retrieve two UniProt sequences. The sequences are returned in a
list variable *seqs*. As you learnt in previous practicals (see
Practical 1 at
`http://www.ucc.ie/microbio/MB6301/practical1\_words\_v2.html <http://www.ucc.ie/microbio/MB6301/practical1_words_v2.html>`_),
to access the elements in an R list variable, you need to use
double square brackets. Therefore, the second element of the list
variable is accessed by typing *seqs[[2]]*. Each element of the
list variable *seqs* contains a vector which stores one of the
sequences.

Once you have retrieved the sequences using retrieveuniprotseqs(),
you can then use the function write.fasta() from the SeqinR library
to write the sequences to a FASTA-format file. As its arguments
(inputs), the write.fasta() function takes the list variable
containing the sequences, and a vector containing the names of the
sequences, and the name that you want to give to the FASTA-format
file. For example:

::

    > write.fasta(seqs, seqnames, file="myseqs.fasta")

The command above will write the sequences in list variable *seqs*
to a FASTA-format file called "myseqs.fasta" in the "My Documents"
folder on your computer.

Note that the retrieveuniprotseqs() function will not work if you
are using a computer that does not have a direct internet
connection (for example, it will not work if you are using some of
the computers in UCC campus teaching labs, as do not connect
directly to the internet, but rather connect through a proxy). In
this case, you will have to search and retrieve the DNA or protein
sequences via the website for the sequence database (for example,
via the NCBI or UniProt website).

The retrieveuniprotseqs() function can be used to retrieve protein
sequences from the UniProt protein sequence database. The file
"Rfunctions.R" also contains a very similar function
retrievegenbankseqs() to retrieve mRNA or DNA sequences from the
NCBI sequence database (the subsection of the NCBI sequence
database that contains DNA sequences is commonly known as
"GenBank").

Creating a multiple alignment of protein, DNA or mRNA sequences
---------------------------------------------------------------

For most phylogenetic analyses, it is necessary to first make a
multiple alignment of the sequences that you are including in your
analysis. You can make a multiple alignment of DNA, mRNA, or
protein sequences using the T-Coffee software, which can be run via
the website
`www.ebi.ac.uk/Tools/t-coffee/ <http://www.ebi.ac.uk/Tools/t-coffee/index.html>`_.

For example, say you want to make a multiple alignment of the human
ASPM protein (a protein involved in brain development; UniProt
accession Q9NVS1) and its homologs in other mammals (eg. UniProt
P62293 is the chimp homolog; Q8CJ27 is the mouse homolog, UniProt
P62286 is the dog homolog; UniProt P62285 is the cow homolog; and
P62297 is the sheep homolog).

We can first retrieve these sequences from UniProt using the
retrieveuniprotseqs() function, and then write them to a
FASTA-format file using the write.fasta() function:

::

    > aspmnames <- c("Q9NVS1", "P62293", "Q8CJ27", "P62286", "P62285", "P62297")
    > aspmseqs <- retrieveuniprotseqs(aspmnames)  # Retrieve the sequences and store them in list variable "aspmseqs"
    > length(aspmseqs)                            # Find the number of sequences that were retrieved
    [1] 6

Note that the function retrieveuniprotseqs() may take a minute or
so to run if you are retrieving a list of sequences, as it takes
time to retrieve each sequence from UniProt. This will write the
sequences in list variable *aspmseqs* to a FASTA-format file called
"aspm.fasta" in the "My Documents" folder on your computer.

You can then make a multiple alignment of these protein sequences
using the T-Coffee multiple alignment software. To do this, go to
the website
`www.ebi.ac.uk/Tools/t-coffee/ <http://www.ebi.ac.uk/Tools/t-coffee/index.html>`_.
On this website you will see a box for pasting in a sequence, and
below the box there will be some text saying "Upload a file". Click
on the "Browse" button to the right of "Upload a file", and select
the FASTA file containing your sequences (eg. "aspm.fasta"). Then
press the red "Run" button to the right of the "Browse" button, to
run T-Coffee.

When T-Coffee has finished running, you will see a page entitled
"T-Coffee Results". You can view the multiple alignment by clicking
on the "Start Jalview button" in the middle of the page, which will
display the alignment in the Jalview alignment viewer. The
alignment displayed in Jalview has a row for each of your
sequences. Jalview colours sets of chemically similar amino acids
in similar colours. For example, tyrosine (Y) is coloured
blue-green, while the chemically similar amino acid phenylalanine
(F) is coloured blue. You can scroll to the right and left along
the alignment using the scrollbar at the bottom of the Jalview
window.

|image0|

To download the T-Coffee alignment, on the "T-Coffee Results"
webpage, right-click (ie. click with your right mouse button) on
the link beside "Phylip tree file" on the fifth line down of the
table at the top of the Results page. Save the alignment file with
a sensible name, eg. "aspm.phy", in the "My Documents" folder on
your computer.

Reading a multiple alignment file into R
----------------------------------------

To read a protein sequence alignment into R from a file, you can
use the read.alignment() function in the SeqinR library. For
example, to read in the multiple sequence alignment of ASPM
proteins, we type:

::

    > aspmaln  <- read.alignment(file = "aspm.phy", format = "phylip")

The *aspmaln* variable is a list variable that stores the
alignment.

As you learnt in previous practicals (see Practical 1, at
`http://www.ucc.ie/microbio/MB6301/practical1\_words\_v2.html <http://www.ucc.ie/microbio/MB6301/practical1_words_v2.html>`_),
an R list variable can have named elements, and you can access the
named elements of a list variable by typing the variable name,
followed by "$", followed by the name of the named element.

The list variable *aspmaln* has named elements "nb", "nam", "seq",
and "com". In fact, the named element "seq" contains the alignment,
which you can view by typing:

::

    > aspmaln$seq
    [[1]]
    [1] "---------------------------------------------------------------------pnee...
    
    [[2]]
    [1] "matrragr-swevspsgprpaa------geaaaasppvlslshfcrspflcfgdvrlggsrtlplllhnpnde...
    
    [[3]]
    [1] "manrrvgrgcwevspterrppaglrgpaaeeeassppvlslshfcrspflcfgdvllgasrtlslaldnpnee...
    
    [[4]]
    [1] "matmqaas-cpeergrrarp-------dpeagdpsppvlllshfcgvpflcfgdvrvgtsrtrslvlhnphee...
    
    [[5]]
    [1] "manrrvgrgcwevspterrppaglrgpaaeeeassppvlslshfcrspflcfgdvllgasrtlslaldnpnee...
    
    [[6]]
    [1] "---------------------------------------------------------------------pnee...

Only the first part of the alignment stored in *aspm$seq* is shown
here, as it is very long.

Calculating genetic distances between protein sequences
-------------------------------------------------------

A common first step in performing a phylogenetic analysis is to
calculate the pairwise genetic distances between sequences. The
genetic distance is an estimate of the divergence between two
sequences, and is usually measured in quantity of evolutionary
change (eg. number of mutations).

We can calculate the genetic distances between protein sequences
using the dist.alignment() function in the SeqinR library. The
dist.alignment() function takes a multiple alignment as input.
Based on the multiple alignment that you give it, dist.alignment()
calculates the genetic distance between each pair of proteins in
the multiple alignment. For example, to calculate genetic distances
between the ASPM proteins based on the multiple sequence alignment,
we type:

::

    > aspmaln  <- read.alignment(file = "aspm.phy", format = "phylip") # Read in the alignment
    > aspmdist <- dist.alignment(aspmaln)                              # Calculate the genetic distances
    > aspmdist                                                         # Print out the genetic distance matrix
               P62285     P62286     P62293     P62297     Q8CJ27    
    P62286     0.32615432                                            
    P62293     0.32062588 0.27259970                                 
    P62297     0.13809534 0.32468430 0.32651352                      
    Q8CJ27     0.37921767 0.36982995 0.35654846 0.38301128           
    Q9NVS1     0.32294090 0.27419231 0.06345439 0.32878574 0.35789648

The genetic distance matrix above shows the genetic distance
between each pair of proteins. Based on the genetic distance matrix
above, we can see that the genetic distance between the cow and
sheep ASPM proteins (P62285 and P62297) by looking at the cell at
the intersection of the column for P62285 (the first column) and
the row for P62297 (the third row), and see that it is 0.13809534.
Similarly, the genetic distance between the human and cow ASPM
proteins (Q9NVS1 and P62285) is in the cell at the intersection of
the column for P62285 (the first column) and the row for Q9NVS1
(the last row), and is 0.32422538.

The larger the genetic distance between two sequences, the more
amino acid changes that have occurred since they shared a common
ancestor, and the longer ago their common ancestor probably lived.
The genetic distance between the human and cow ASPM proteins is
larger than the genetic distance between the sheep and cow ASPM
proteins, indicating that more amino acid changes have occurred in
the human and cow ASPM proteins since they shared a common
ancestor, compared in the sheep and cow ASPM proteins since they
shared a common ancestor.

Building an unrooted phylogenetic tree for protein sequences based on a distance matrix
---------------------------------------------------------------------------------------

Once we have a distance matrix that gives the pairwise distances
between all our protein sequences, we can build a phylogenetic tree
based on that distance matrix. One method for using this is the
neighbour-joining algorithm.

You can build a phylogenetic tree using the neighbour-joining
algorithm with the nj() function the Ape R library. The nj()
function takes a distance matrix as its argument (input), and
builds a phylogenetic tree.

::

    > library("ape")
    > aspmaln  <- read.alignment(file = "aspm.phy", format = "phylip") # Read in the alignment
    > aspmdist <- dist.alignment(aspmaln)                              # Calculate the genetic distance matrix
    > aspmtree <- nj(aspmdist)                                         # Calculate the neighbour-joining tree   

After building a neighbour-joining tree, we can then plot a picture
of the tree using the plot.phylo() function from the Ape library.
The plot.phylo() function has an argument "type", which tells it
what sort of tree you want. For example, if a tree does not contain
an outgroup, then it is an unrooted tree, and you can tell
plot.phylo() to draw an unrooted tree by using the type="u"
argument.

For example, to plot a picture of the unrooted phylogenetic tree of
ASPM proteins, we type:

::

                            
    > plot.phylo(aspmtree, type="u") # Plot the tree

|image1|

In the plot of the phylogenetic tree, pairs of sequences that
dist.alignment() calculated as having small pairwise genetic
distances should be close together in the tree, while pairs of
sequences that dist.alignment() calculated as having large pairwise
genetic distances should be further apart in the tree. For example,
the human and cow ASPM proteins (Q9NVS1 and P62285), which are
separated by a relatively large genetic distance, are further apart
in the tree than the sheep and cow ASPM proteins (P62297 and
P62285), which are separated by a relatively small genetic distance
(see above).

Furthermore, the lengths of the branches in the plot of the tree
are proportional to the evolutionary change along the the branches.
Thus, we can see from the tree above that the human and chimp ASPM
proteins (Q9NVS1 and P62293) are more closely related to each other
than to any other ASPM proteins, and that the genetic distances
between these two proteins and their last common ancestor node are
relatively small compared to the other genetic distances in the
tree (ie. lengths of the branches, shown in red in the plot above,
are short compared to other branch lengths in the tree).

Finding an outgroup to make a rooted phylogenetic tree
------------------------------------------------------

The tree above of the ASPM proteins is an *unrooted* phylogenetic
tree as it does not contain an outgroup sequence. As a result, we
cannot tell which direction evolutionary time ran in along the
internal branches of the tree.

In order to convert the unrooted tree into a rooted tree, we need
to add an outgroup sequence. Normally, the outgroup sequence is a
sequence that we know from some prior knowledge to be more
distantly related to the other sequences under study than they are
to each other.

For example, as an outgroup to the ASPM proteins, we could add the
ASPM homolog from the zebrafish (UniProt accession Q1L925). It is
well known from fossil and morphological evidence that the
zebrafish is more distantly related to the other species under
study (human, chimp, dog, mouse, cow and sheep) than they are to
each other, so we can assume that the zebrafish ASPM protein is
more distantly related to the ASPM human, chimp, dog, mouse, cow
and sheep ASPM proteins than they are to each other. Therefore, the
zebrafish ASPM protein is a suitable outgroup for the tree of ASPM
proteins.

To add an outgroup sequence to the tree, we need to first retrieve
the outgroup sequence from the database (UniProt here). To do this
for the ASPM proteins, we type:

::

    > aspmnames2 <- c("Q9NVS1", "P62293", "Q8CJ27", "P62286", "P62285", "P62297", "Q1L925")
    > aspmseqs2 <- retrieveuniprotseqs(aspmnames2)  # Retrieve the sequences and store them in list variable "aspmseqs2"
    > length(aspmseqs2)                             # Find the number of sequences that were retrieved
    [1] 7 
    > write.fasta(aspmseqs2, aspmnames2, file="aspm2.fasta")

We then need to build a new alignment of the sequences including
the outgroup sequences, for example, using T-Coffee as described
above (eg. to make file "aspm2.phy").

Building a rooted phylogenetic tree for protein sequences based on a distance matrix
------------------------------------------------------------------------------------

To build a rooted phylogenetic tree that contains the outgroup
sequence, we need to build a new distance matrix based on the new
alignment containing the outgroup, and then a new tree based on the
new distance matrix:

::

    > aspmaln2  <- read.alignment(file = "aspm2.phy", format = "phylip")
    > aspmdist2 <- dist.alignment(aspmaln2)
    > aspmdist2 # Print out the genetic distance matrix
               P62285     P62286     P62293     P62297     Q1L925     Q8CJ27    
    P62286     0.32574573                                                       
    P62293     0.32104239 0.27210625                                            
    P62297     0.13809534 0.32427367 0.32737556                                 
    Q1L925     0.51374081 0.51107841 0.50982780 0.51253811                      
    Q8CJ27     0.37928066 0.36994931 0.35604085 0.38307486 0.52414242           
    Q9NVS1     0.32335374 0.27370220 0.06345439 0.32964116 0.51044540 0.35739034
    > aspmtree2 <- nj(aspmdist2)
    > aspmtree2$tip.label # Print out the names of the sequences in the tree
    [1] "P62285    " "P62286    " "P62293    " "P62297    " "Q1L925    " "Q8CJ27    " "Q9NVS1    "

The last line of the commands above prints out the names of the
sequences in the tree *aspmtree2* (this is because *aspmtree2* is a
list variable that has a named element "tip.label" containing the
names of the sequences in the tree). The sequence names may include
some extra spaces when they are stored in a phylogenetic tree such
as *aspmtree2*, for example, the zebrafish protein's name is stored
as "Q1L925 " (with 4 spaces after the accession).

Once we have built a new tree based on the new distance matrix, we
need to tell R that it is a tree with an outgroup, that is, a
rooted tree. This can be done using the root() function from the R
Ape library. The root() function takes as its argument (input) the
name of the sequence that you want to be the outgroup in the tree
(the zebrafish protein Q1L925 here). We need to give the root()
function the name for the outgroup that is used in the tree, for
example, "Q1L925 " (with the 4 extra spaces after the accession).
This is necessary so that R realises which sequence in the tree you
want to be the outgroup

For example, to make a bifurcating rooted tree of the ASPM
proteins, we type:

::

    > rootedaspmtree2 <- root(aspmtree2,"Q1L925    ",r=TRUE) # Specify that Q1L925 is the outgroup.
    > plot.phylo(rootedaspmtree2) 

|image2|

The above tree shows the zebrafish protein Q1L925 as the outgroup
to the ASPM protein tree. As this is a rooted tree, we know the
direction that evolutionary time ran: from left to right in this
case. Thus, we can infer from the tree that the human, chimp and
dog proteins (human Q9NVS1, chimp P62293, and dog P62286) shared a
common ancestor with each other more recently than they did with
the other ASPM proteins in the tree. In addition, the sheep and cow
ASPM proteins (sheep P62297 and cow P62285) shared a common
ancestor with each other more recently than they did with the other
ASPM proteins in the tree.

The lengths of branches in this tree are proportional to the amount
of evolutionary change that occurred along the branches. The
branches leading back from the sheep and cow ASPM proteins to their
last common ancestor (coloured blue) are slightly longer than the
branches leading back from the chimp and human ASPM proteins to
their last common ancestor (coloured red). This indicates that
there has been more evolutionary change in the sheep and cow ASPM
proteins since they diverged, than there has been in the chimp and
human ASPM proteins since they diverged.

Building a phylogenetic tree with bootstrap values
--------------------------------------------------

The above tree gives us an idea of the evolutionary relationships
between the ASPM proteins. However, if we want to know how
confident we are in each part of the tree, it is necessary to build
a phylogenetic tree with bootstrap values.

The bootstrap values are calculated by making many (for example,
100) random "resamples" of the alignment that the phylogenetic tree
was based upon. Each "resample" of the alignment consists of a
certain number *x* (eg. 200) of randomly sampled columns from the
alignment. Each "resample" of the alignment (eg. 200 randomly
sampled columns) forms a sort of fake alignment of its own, and a
phylogenetic tree can be based upon the "resample". We can make 100
random resamples of the alignment, and build 100 phylogenetic trees
based on the 100 resamples. These 100 trees are known as the
"bootstrap trees". For each clade that we see in our original
phylogenetic tree, we can count in how many of the 100 bootstrap
trees it appears. This is known as the "bootstrap value" for the
clade in our original phylogenetic tree.

For example, if we calculate 100 random resamples of the ASPM
protein alignment, and build 100 phylogenetic trees based on these
resamples, we can calculate the bootstrap values for each clade in
the ASPM phylogenetic tree. For example, in the tree above, we saw
a clade consisting of chimp ASPM, P62293, and human ASPM, Q9NVS1
(shown in red). The bootstrap value for this clade is the number of
the bootstrap trees that this clade appears in.

The bootstrap values for a phylogenetic tree can be calculated in R
using the boot.phylo() function in the Ape R library. By default,
the boot.phylo() function calculates the bootstrap values based on
100 bootstrap trees. The boot.phylo() function takes as arguments:


-  the original tree that we want to add bootstrap values to (eg.
   the tree *rootedaspmtree2*)
-  the alignment in the form of a matrix of characters
-  the function to use to build both the original tree (eg.
   *rootedaspmtree2*) and the bootstrap trees

If you look at the help page for the boot.phylo() function, you
will see that it requires its second argument (input) to be the
alignment in the form of a matrix of characters with one row per
sequence and one column per alignment column. Normally, when you
read in an alignment using the read.alignment() function, it is
stored as a list variable that has named elements "nb", "nam",
"seq", and "com". As discussed above, the named element "seq"
stores the alignment. To convert this list variable into an
alignment in the form of a matrix of characters, we can use the
as.matrix.alignment() function from the SeqinR library:

::

    > aspmaln2        <- read.alignment(file = "aspm2.phy", format = "phylip")      # Read in the alignment
    > aspmaln2mat     <- as.matrix.alignment(aspmaln2)                              # Convert alignment to a matrix of characters
    > aspmaln2mat                                                                   # Print out aspmaln2mat
              1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20 
    P62285     "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    P62286     "m" "a" "t" "r" "r" "a" "g" "r" "-" "s" "w" "e" "v" "s" "p" "s" "g" "p" "r" "p"
    P62293     "m" "a" "n" "r" "r" "v" "g" "r" "g" "c" "w" "e" "v" "s" "p" "t" "e" "r" "r" "p"
    P62297     "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    Q1L925     "m" "s" "f" "k" "v" "a" "k" "s" "e" "c" "l" "d" "f" "s" "p" "p" "l" "d" "s" "h"
    Q8CJ27     "m" "a" "t" "m" "q" "a" "a" "s" "-" "c" "p" "e" "e" "r" "g" "r" "r" "a" "r" "p"
    Q9NVS1     "m" "a" "n" "r" "r" "v" "g" "r" "g" "c" "w" "e" "v" "s" "p" "t" "e" "r" "r" "p"
    ...

In the matrix of characters representing the alignment, each column
of the matrix represents one column of the alignment, and each row
represents one row in the alignment. Only the start of the matrix
of characters *aspmaln2mat* is printed out above, as it is very
large. If we have an alignment in the format of a matrix of
characters, we can convert it back into a list variable by using
the as.alignment function from the Ape library, for example:

::

    > aspmaln2b       <- as.alignment(aspmaln2mat) # Convert the matrix of characters into a list variable aspmaln2b

To build a rooted phylogenetic tree with bootstrap values using
boot.phylo(), we can first define the function that we want to use
to build the tree, for example:

::

    > myrootedprotnjtree <- function(alignmentmat)
     {
        alignment  <- as.alignment(alignmentmat)   # Convert alignmentmat into the format required by dist.alignment() 
        distmat    <- dist.alignment(alignment)    # Calculate the genetic distance matrix
        tree       <- nj(distmat)                  # Calculate the neighbour-joining tree
        rootedtree <- root(tree, "Q1L925    " , r=TRUE) # Convert the tree into a rooted tree
        return(rootedtree)
     }

This function builds a rooted phylogenetic tree using the zebrafish
sequence Q1L925 as the outgroup. It takes as its argument a matrix
of characters representing the alignment, that is, the alignment in
the format produced by function as.matrix.alignment(). The
dist.alignment() function requires as its argument the alignment in
the form of a list variable, so we use the as.alignment() function
to convert the matrix of characters representing the alignment
*alignmentmat* into a list variable format.

Once we have defined the function that we want to use to build a
phylogenetic tree, we can then build a rooted phylogenetic tree of
ASPM proteins by typing:

::

    > aspmaln2        <- read.alignment(file = "aspm2.phy", format = "phylip")        # Read in the alignment
    > aspmaln2mat     <- as.matrix.alignment(aspmaln2)                                # Convert the alignment to the format required by boot.phylo()
    > rootedaspmtree2 <- myrootedprotnjtree(aspmaln2mat)                              # Build a rooted phylogenetic tree

We can then calculate bootstrap values for the rooted phylogenetic
tree of ASPM proteins using the boot.phylo() function, by typing:

::

    > aspmboot        <- boot.phylo(rootedaspmtree2, aspmaln2mat, myrootedprotnjtree) # Calculate the bootstrap values as percentages
    > aspmboot                                                                        # Print out the bootstrap values as percentages
    [1] 100 100 100 100 100 100

We can then plot the tree using the plot.phylo() function, and
display the bootstrap values as percentages on the nodes of the
tree using the nodelabels() function from the Ape library, by
typing:

::

    > plot.phylo(rootedaspmtree2)
    > nodelabels(aspmboot)

|image3|

In the plot above, the bootstrap value for the clade containing the
human ASPM (Q9NVS1) and chimp ASPM (P62293) proteins is 100 (100%),
which means that this clade occurred in 100% of the bootstrap
trees.

Calculating genetic distances between DNA or mRNA sequences
-----------------------------------------------------------

In the example above, a phylogenetic tree was built of ASPM protein
sequences from vertebrates. The genomes of distantly related
organisms such as vertebrates will have accumulated many mutations
since they diverged. Sometimes, so many mutations have occurred
since the organisms diverged that their DNA sequences are hard to
align correctly and it is also hard to accurately estimate
evolutionary distances from alignments of those DNA sequences. In
contrast, as many mutations at the DNA level are synonymous at the
protein level, protein sequences diverge at a slower rate than DNA
sequences. This is why for reasonably distantly related organisms
such as vertebrates, it is usually preferable to use protein
sequences for phylogenetic analyses.

If you are studying closely related organisms such as primates, few
mutations will have occurred since they diverged. As a result, if
you use protein sequences for a phylogenetic analysis, there may be
too few amino acid substitutions to provide enough 'signal' to use
for the phylogenetic analysis. Therefore, it is often preferable to
use DNA sequences for a phylogenetic analysis of closely related
organisms such as primates.

One example where this is the case is phylogenetic analysis of
*aspm* genes from primates. The NCBI Sequence Database contains
*aspm* gene mRNA sequences from various primates, including human
(NCBI accession AF509326), chimp (NCBI accession AY367066),
orangutan (NCBI accession AY367067), and gorilla (NCBI accession
AY508451). We can retrieve these sequences and save them to a FASTA
format file "aspm3.fasta" using SeqinR, by using the
retrievegenbankseqs() function:

::

    > aspmnames3 <- c("AF509326","AY367066","AY367067","AY508451")
    > aspmseqs3 <- retrievegenbankseqs(aspmnames3)
    > length(aspmseqs3) # Print out the number of sequences
    [1] 4 
    > write.fasta(aspmseqs3, aspmnames3, file="aspm3.fasta")

We can then use T-Coffee to build a multiple alignment of these DNA
sequences, and save the alignment file (eg. as "aspm3.phy"), as
described above.

Building a phylogenetic tree for DNA or mRNA sequences based on a distance matrix
---------------------------------------------------------------------------------

To carry out a phylogenetic analysis based on DNA sequences, you
need to use slightly different methods for calculating a genetic
distance matrix than used for protein sequences. You can calculate
a genetic distance for DNA sequences using the dist.dna() function
in the Ape R library. dist.dna() takes a multiple alignment of DNA
sequences as its input, and calculates the genetic distance between
each pair of DNA sequences in the multiple alignment. The
dist.dna() function requires the input alignment to be in a special
format known as "DNAbin" format, so we must use the as.DNAbin()
function to convert our DNA alignment into this format before using
the dist.dna() function. For example, to calculate the genetic
distance between each pair of DNA sequences in an alignment file
"aspm3.phy", we type:

::

    > aspmaln3        <- read.alignment(file = "aspm3.phy", format = "phylip") # Read in the alignment
    > aspmaln3DNAbin  <- as.DNAbin(aspmaln3)      # Convert the alignment into "DNAbin" format
    > aspmdist3       <- dist.dna(aspmaln3DNAbin) # Calculate the genetic distance matrix
    > aspmdist3                                   # Print out the genetic distance matrix
                AF509326    AY367066    AY367067
    AY367066 0.005824745                        
    AY367067 0.013771317 0.014465548            
    AY508451 0.007675178 0.009243747 0.014559095

Once you have built a distance matrix that gives the pairwise
distances between all your DNA sequences, you can use the nj()
function to build a phylogenetic tree based on that distance matrix
using the neighbour-joining algorithm. For example, to build a
phylogenetic tree of the *aspm* mRNA sequences, using the orangutan
sequence AY367067 as the outgroup, we type:

::

    > aspmaln3tree       <- nj(aspmdist3)                                         # Calculate the neighbour-joining tree
    > rootedaspmtree3    <- root(aspmaln3tree, "AY367067", r=TRUE)                # Convert the tree to a rooted tree
    > plot.phylo(rootedaspmtree3)                                                 # Plot the tree

|image4|

The orangutan sequence (accession AY367067) is used as the
outgroup, as we know from prior knowledge that orangutans are more
distantly related to chimp, human and gorilla than they are to each
other.

If you want to add bootstrap values to a rooted phylogenetic tree
based on a DNA or mRNA sequence alignment, you can easily do that
using the boot.phylo() function. You would first need to define a
function that built a rooted tree with a certain outgroup, for
example, here is a function to build a rooted tree using the
orangutan "AY367067" sequence as the outgroup:

::

    > myrooteddnanjtree <- function(alnbin)
     {
        distmat    <- dist.dna(alnbin)                # Calculate the genetic distance matrix
        tree       <- nj(distmat)                     # Calculate the neighbour-joining tree
        rootedtree <- root(tree, "AY367067", r=TRUE)  # Convert the tree into a rooted tree
        return(rootedtree)
     }

For example, to add bootstrap values to the phylogenetic tree of
*aspm* mRNA sequences, we type:

::

    > aspmboot3          <- boot.phylo(rootedaspmtree3, aspmaln3DNAbin, myrooteddnanjtree) # Calculate the bootstrap values as percentages
    > aspmboot3                                                                            # Print the bootstrap values as percentages 
    [1] 100  39  99 
    > plot.phylo(rootedaspmtree3)                                                          # Make a plot of the rooted tree
    > nodelabels(aspmboot3)                                                                # Add the bootstrap values as labels to the nodes

|image5|

We can see from the tree that the human and chimp *aspm* genes
(accessions AF509326 and AY367066) shared a common ancestor with
each other more recently than they did with the gorilla *aspm* gene
(AY508451).

Summary
-------

In this practical, you have learnt the following R functions that
belong to the bioinformatics libraries:


#. install.packages() for installing an R library (except for
   Bioconductor R libraries), if you have a direct internet connection
#. retrieveuniprotseqs() from "Rfunctions.R", which uses SeqinR to
   retrieve protein sequences from UniProt
#. retrievegenbankseqs() from "Rfunctions.R", which uses SeqinR to
   retrieve DNA or mRNA sequences from NCBI
#. write.fasta() from the SeqinR library for writing sequences to a
   FASTA-format file
#. read.alignment() from the SeqinR library for reading in a
   multiple alignment
#. dist.alignment() from the SeqinR library for calculating genetic
   distances between protein sequences
#. nj() from the Ape library for building a neighbour-joining tree
#. plot.phylo() from the Ape library for plotting a phylogenetic
   tree
#. root() from the Ape library for converting an unrooted tree to a
   rooted tree
#. as.matrix.alignment() from the SeqinR library for converting an
   alignment in the form of a list variable to an alignment in the
   form of a matrix of characters
#. as.alignment() from the Ape library for converting an alignment
   in the form of a matrix of characters to an alignment int he form
   of a list variable
#. boot.phylo() from the Ape library for calculating bootstrap
   values for a tree
#. nodelabels() from the Ape library for adding labels to the nodes
   of a tree in a tree plot
#. as.DNAbin() from the Ape library for convering an alignment in
   the form a a list to the "DNAbin" format required by the dist.dna()
   function
#. dist.dna() from the Ape library for calculating genetic
   distances between DNA or mRNA sequences

Links and Further Reading
-------------------------

Some links are included here for further reading, which will be
especially useful if you need to use the R package and SeqinR and
Ape libraries for your project or assignments.

For background reading on phylogenetic trees, it is recommended to
read Chapter 7 of
*Introduction to Computational Genomics: a case studies approach*
by Cristianini and Hahn (Cambridge University Press;
`www.computational-genomics.net/book/ <http://www.computational-genomics.net/book/>`_).

For more in-depth information and more examples on using the SeqinR
library for sequence analysis, look at the SeqinR documentation,
`seqinr.r-forge.r-project.org/seqinr\_2\_0-1.pdf <http://seqinr.r-forge.r-project.org/seqinr_2_0-1.pdf>`_.

For more in-depth information and more examples on the Ape library
for phylogenetic analysis, look at the Ape documentation,
`ape.mpl.ird.fr/ <http://ape.mpl.ird.fr/>`_.

If you are using the Ape library for a phylogenetic analysis
project, it would be worthwhile to obtain a copy of the book
*Analysis of Phylogenetics and Evolution with R* by Emmanuel
Paradis, published by Springer, which has many nice examples of
using R for phylogenetic analyses.

Acknowledgements
----------------

Many of the ideas for the examples and exercises for this practical
were inspired by the Matlab case study on SARS
(`www.computational-genomics.net/case\_studies/sars\_demo.html <http://www.computational-genomics.net/case_studies/eyeless_demo.html>`_)
from the website that accompanies the book
*Introduction to Computational Genomics: a case studies approach*
by Cristianini and Hahn (Cambridge University Press;
`www.computational-genomics.net/book/ <http://www.computational-genomics.net/book/>`_).

Thank you to Jean Lobry and Simon Penel for helpful advice on using
the SeqinR library.

Thank you to Emmanuel Paradis and François Michonneau for help in
using the Ape library.

Exercises
---------

Answer the following questions, using the R package. For each
question, please record your answer, and what you typed into R to
get this answer.

Q1. Calculate the genetic distances between the following Spike
proteins from different coronaviruses:

-  bovine coronavirus CoV1 Spike protein (UniProt Q8V436)
-  bovine coronavirus CoV2 Spike protein (UniProt Q91A26)
-  human coronavirus OC43 Spike protein (UniProt P36334)
-  porcine coronavirus HEV3 Spike protein (UniProt Q8BB25)
-  murine coronavirus HV2 Spike protein (UniProt P11224)
-  avian coronavirus IBV3 Spike protein (UniProt P11223)
-  porcine coronavirus PEDV Spike protein (UnniProt Q91AV1)
-  canine coronavirus CoV1 Spike protein (UniProt Q65984)
-  feline coronavirus CoV4 Spike protein (UniProt Q66951)
-  human SARS coronavirus CoV Spike protein (UniProt P59594)
-  palm civet coronavirus Spike protein (UniProt Q5GDB3)

Which protein is has the smallest genetic distance from the human
SARS Spike protein?
SARS (Severe Acute Respiratory Syndrome) is a human illness that
first appeared in late 2002 in Guangdong Province, China. It is now
known that the disease is caused by the SARS coronavirus
(SARS-CoV), a novel coronavirus.
Q2. Build an unrooted phylogenetic tree with bootstrap values of
the proteins from Q1, using the neighbour-joining algorithm. Based
on the phylogenetic tree for the coronavirus Spike proteins, which
coronavirus do you think that human SARS is most closely related
to?
Based on the bootstrap values in the tree, how confident are you of
this?
Q3. Calculate an unrooted phylogenetic tree with bootstrap values
of the following Spike gene DNA sequences from human SARS viruses
that were isolated from infected patients:

-  isolated from a patient in Guangzhou (Guangdong Province, China)
   on 16th Dec 2002
-  isolated from a patient in Zhongshan (Guangdong Province, China)
   on 26th Dec 2002
-  isolated from a patient in Zhongshan (Guangdong Province, China)
   on 4th Jan 2003
-  isolated from a patient in Guangzhou (Guangdong Province, China)
   on 24th Jan 2003
-  isolated from a patient in Guangzhou Hospital (Guangdong
   Province, China) on 31st Jan 2003
-  isolated from a patient in Guangzhou (Guangdong Province, China)
   on 18th Feb 2003
-  isolated from a patient in Hong Kong on 21st Feb 2003
-  isolated from a patient in Hanoi, Vietnam on 26th Feb 2003
-  isolated from a patient in Toronto, Canada on 27th Feb 2003
-  isolated from a patient in Singapore on 1st Mar 2003
-  isolated from a patient in Taiwan, on 8th Mar 2003
-  isolated from a patient in Hong Kong, on 19th Mar 2003
-  isolated from a patient in Hong Kong, on 15th Mar 2003

To save you time, we have already made a FASTA-format file containing these DNA sequences, called `sars\_spike.fasta <http://www.ucc.ie/ucc/depts/microbio/MB6301/sars_spike.fasta>`_, which you can download to use for this analysis. Based on the Spike gene DNA phylogenetic tree, what is the relationship between the palm civet coronavirus and the human SARS isolates? 
    Would the Spike gene sequence from palm civet make a suitable
    outgroup, and why?
    Make a rooted tree using the palm civet Spike gene as the outgroup.
Q4. Based on your phylogenetic tree from Q3, is the palm civet coronavirus more closely related to human SARS isolates that were isolated early or late in the epidemic? 
    Note: the date and place that each sample was collected should be
    recorded in its name, eg. the sample labelled '03Feb26Han' was
    collected on 26th February 2003 in Hanoi.
    What does this tell us about the history of the epidemic (eg. place
    and time of origin of the human SARS virus)?
Q5. Based on your phylogenetic tree from Q3, what is the relationship between the human SARS isolates from the Metropole Hotel in Hong Kong and (i) those in Guangdong province? (ii) those in other world cities (Taiwan, Hanoi, Toronto, Singapore)? 
    Note: the samples collected in the Metropole hotel in Hong Kong are
    labelled 'YearMonthDateHon' ('Hon' stands for 'HongKong' here).
    What role did people who stayed in the Metropole hotel probably
    play in the spread of SARS?

Other ways to do the same thing
-------------------------------

It is possible to carry out some of the analyses that you have
carried out in the practicals via websites. You can download the
T-Coffee alignment program from
`www.tcoffee.org/Projects\_home\_page/t\_coffee\_home\_page.html <http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html>`_
and run it on your own computer.

It is possible to calculate genetic distances between protein
sequences using the Protdist program, via the website
`mobyle.pasteur.fr/cgi-bin/portal.py?form=protdist <http://mobyle.pasteur.fr/cgi-bin/portal.py?form=protdist>`_.
Similarly, it is possible to calculate genetic distances between
DNA or mRNA sequences using the DNAdist program, via the website
`mobyle.pasteur.fr/cgi-bin/portal.py?form=dnadist <http://mobyle.pasteur.fr/cgi-bin/portal.py?form=dnadist>`_.
You can build a phylogenetic tree based on a genetic distance
matrix with the neighbour-joining algorithm by using the Neighbor
program, via the website
`mobyle.pasteur.fr/cgi-bin/portal.py?form=neighbor <http://mobyle.pasteur.fr/cgi-bin/portal.py?form=neighbor>`_.
A nice program for plotting a phylogenetic tree produced by the
Neighbor program is Phylodendron, available at the website
`iubio.bio.indiana.edu/treeapp/treeprint-form.html <http://iubio.bio.indiana.edu/treeapp/treeprint-form.html>`_.

The Protdist, Dnadist and Neighbor programs are also available for
download as part of the PHYLIP package for phylogenetic analysis
(`evolution.genetics.washington.edu/phylip.html <http://evolution.genetics.washington.edu/phylip.html>`_)
and so can also be run on your own computer.

As well as the R Ape library and PHYLIP, there are a large number
of other software packages available for phylogenetic analyses. Joe
Felsenstein maintains a very useful list of phylogenetic software
packages on his website at
`evolution.gs.washington.edu/phylip/software.html <http://evolution.gs.washington.edu/phylip/software.html>`_.




.. |image0| image:: ../_static/P5_image0.png
.. |image1| image:: ../_static/P5_image2b.png
.. |image2| image:: ../_static/P5_image3b.png
.. |image3| image:: ../_static/P5_image4.png
.. |image4| image:: ../_static/P5_image7b.png
.. |image5| image:: ../_static/P5_image7.png
