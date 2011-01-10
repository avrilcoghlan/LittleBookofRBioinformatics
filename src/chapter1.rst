Chapter 1: DNA Sequence Statistics
==================================

Using R for Bioinformatics 
--------------------------

This booklet tells you how to use the R software to carry out some simple analyses
that are common in bioinformatics. In particular, the focus is on computational analysis
of biological sequence data such as genome sequences and protein sequences.

This booklet assumes that the reader has some basic knowledge of biology, but not
necessarily of bioinformatics. The focus of the booklet is to explain simple bioinformatics
analysis, and to explain how to carry out these analyses using R.

To use R, you first need to start the R program on your computer.
You should have already installed R on your computer (if not, for instructions on how to
install R, see `How to install R <./installr.html>`_).

Neglected tropical diseases
---------------------------

Neglected tropical diseases are serious diseases that affect many people in
tropical countries and which have been relatively little studied. The World
Health Organisation lists the following as neglected tropical diseases:
trachoma, leprosy, schistosomiasis, soil transmitted helminths, lymphatic
filariasis, and onchocerciasis (see `http://www.neglectedtropicaldiseases.org <http://www.neglectedtropicaldiseases.org/>`_). 

Other diseases that are sometimes included as neglected
diseases are Buruli ulcer, yaws, Chagas disease, African trypanosomiasis,
leishmaniasis, Dengue fever, rabies, Dracunculiasis (guinea-worm disease),
and Fascioliasis.

The genomes of many of the organisms that cause neglected tropical diseases have
been fully sequenced, or are currently being sequenced, including:
* the bacterium <i>Chlamydia trachomatis</i>, which causes trachoma
* the bacterium <i>Mycobacterium leprae</i>, which causes leprosy
* the bacterium <i>Treponema pallidum</i>, which causes yaws
* the bacterium <i>Mycobacterium ulcerans</i>, which causes Buruli ulcer
* the schistosome worm <i>Schistosoma mansoni</i>, which causes schistosomiasis
* the nematode worm <i>Brugia malayi</i>, which causes lymphatic filariasis
* the protist <i>Trypanosoma cruzi</i>, which causes Chagas disease
* the protist <i>Trypanosoma brucei</i>, which causes African trypanosomiasis
* the protist <i>Leishmania major</i>, which causes leishmaniasis
* the Dengue virus, which causes Dengue fever
* the Rabies virus, which causes Rabies

R packages for bioinformatics: Bioconductor and SeqinR
------------------------------------------------------

Many authors have written R packages for performing a wide variety
of analyses. These do not come with the standard R installation,
but must be installed and loaded as "add-ons".

Bioinformaticians have written several specialised *packages* for
R. In this practical, you will learn to use the SeqinR package to
retrieve sequences from a DNA sequence database, and to carry out
simple analyses of DNA sequences.

Some well known bioinformatics packages for R are the Bioconductor
set of R packages  
(`www.bioconductor.org <http://www.bioconductor.org/>`_), which
contains several packages with many R functions for analysing
biological data sets such as microarray data; and the SeqinR
package
(`pbil.univ-lyon1.fr/software/seqinr/home.php?lang=eng <http://pbil.univ-lyon1.fr/software/seqinr/home.php?lang=eng>`_),
which contains R functions for obtaining sequences from DNA and protein
sequence databases, and for analysing DNA and protein sequences.

To use function from the SeqinR package, 
we first need to install the SeqinR package (for instructions on how to
install an R package, see `How to install an R package <./installr.html#how-to-install-an-r-package>`_).
Once you have installed the "SeqinR" R package, you can load the "SeqinR" R package by typing:

.. highlight:: r

::

    > library("seqinr")

Remember that you can ask for more information about a particular R
command by using the help() function. For example, to ask for more
information about the library() function, you can type:

::

    > help("library")

FASTA format
------------

The FASTA format is a simple and widely used format for storing
biological (DNA or protein) sequences. It was first used by the
FASTA program for sequence alignment. It begins with a single-line
description starting with a ">" character, followed by lines of
sequences. Here is an example of a FASTA file:

::

    > A06852 183 residues
    MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEE
    PQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSN
    LNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIAR
    LC

The NCBI sequence database
--------------------------

The National Centre for Biotechnology Information (NCBI)
(`www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov/>`_) in the US
maintains a huge database of all the DNA and protein sequence data
that has been collected, the NCBI Sequence Database. This also a
similar database in Europe, the European Molecular Biology
Laboratory (EMBL) Sequence Database
(`www.ebi.ac.uk/embl <http://www.ebi.ac.uk/embl/>`_), and also a
similar database in Japan, the DNA Data Bank of Japan (DDBJ;
`www.ddbj.nig.ac.jp <http://www.ddbj.nig.ac.jp/>`_). These three
databases exchange data every night, so at any one point in time,
they contain almost identical data.

Each sequence in the NCBI Sequence Database is stored in a separate
*record*, and is assigned a unique identifier that can be used to
refer to that sequence record. The identifier is known as an
*accession*, and consists of a mixture of numbers and letters. For
example, 

xxx
Bacteriophage lambda infects the bacterium
*Escherichia coli*, and was one of the first viral genomes to be
completely sequenced (in 1982). The NCBI accession for the DNA
sequence of the Bacteriophage lambda is NC\_001416.

Note that because the NCBI Sequence Database, the EMBL Sequence
Database, and DDBJ exchange data every night, the Bacteriophage
lambda sequence will be present in all three databases, but it will
have different accessions in each database, as they each use their
own numbering systems for referring to their own sequence records.

Links and Further Reading
-------------------------

Some links are included here for further reading.

For a more in-depth introduction to R, a good online tutorial is
available on the "Kickstarting R" website,
`cran.r-project.org/doc/contrib/Lemon-kickstart <http://cran.r-project.org/doc/contrib/Lemon-kickstart/>`_.

There is another nice (slightly more in-depth) tutorial to R
available on the "Introduction to R" website,
`cran.r-project.org/doc/manuals/R-intro.html <http://cran.r-project.org/doc/manuals/R-intro.html>`_.

Acknowledgements
----------------

Thank you to Noel O'Boyle for helping in using Sphinx, `http://sphinx.pocoo.org <http://sphinx.pocoo.org>`_, to create
this document, and github, `https://github.com/ <https://github.com/>`_, to store different versions of the document
as I was writing it, and readthedocs, `http://readthedocs.org/ <http://readthedocs.org/>`_, to build and distribute
this document.

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.

.. |image4| image:: ../_static/image4.png

