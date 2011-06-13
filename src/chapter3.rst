Sequence Databases
==================

The NCBI Sequence Database
--------------------------

All published genome sequences are available over the internet, as
it is a requirement of every scientific journal that any published
DNA or RNA or protein sequence must be deposited in a public
database. The main resources for storing and distributing sequence
data are three large databases: the NCBI database
(`www.ncbi.nlm.nih.gov/ <http://www.ncbi.nlm.nih.gov/>`_), the
European Molecular Biology Laboratory (EMBL) database
(`www.ebi.ac.uk/embl/ <http://www.ebi.ac.uk/embl/>`_, and the DNA
Database of Japan (DDBJ) database
(`www.ddbj.nig.ac.jp/ <http://www.ddbj.nig.ac.jp/>`_). These
databases collect all publicly available DNA, RNA and protein
sequence data and make it available for free. They exchange data
nightly, so contain essentially the same data.

In this chapter we will discuss the NCBI database. Note however
that it contains essentially the same data as in the EMBL/DDBJ
databases. The NCBI database contains several sub-databases,
including the NCBI Nucleotide database and the NCBI Protein
database. These are both sequence databases. The NCBI database also
contains sub-databases that contain non-sequence data, such as
PubMed, which contains data on scientific publications.

Sequences in the NCBI Sequence Database (or EMBL/DDBJ) are
identified by an accession number. This is a unique number that is
only associated with one sequence. For example, the accession
number NC\_001477 is for the DEN-1 Dengue virus genome
sequence. The accession number is what identifies the sequence. It
is reported in scientific papers describing that sequence.

As well as the sequence itself, for each sequence the NCBI database
(or EMBL/DDBJ databases) also stores some additional *annotation*
data, such as the name of the species it comes from, references to
publications describing that sequence, etc. Some of this annotation
data was added by the person who sequenced a sequence and submitted
it to the NCBI database, while some may have been added later by a
human curator working for NCBI.

Searching for an accession number in the NCBI database
------------------------------------------------------

In the `DNA Sequence Statistics chapter (1) <chapter1.html>`_, 
you learnt how to obtain a FASTA file containing the DNA sequence
corresponding to a particular accession number, eg. accession
number NC\_001477 (the DEN-1 Dengue virus genome sequence), either
`via the NCBI website <./chapter1.html#retrieving-genome-sequence-data-via-the-ncbi-website>`_
or `using the getncbiseq() function in R <./chapter1.html#retrieving-genome-sequence-data-using-seqinr>`_.

As explained in the `DNA Sequence Statistics (1) chapter <chapter1.html#fasta-format>`_, 
the FASTA format is a file format commonly used to store sequence information. The first line starts
with the character '>' followed by a name and/or description for
the sequence. Subsequent lines contain the sequence itself.

::

    >mysequence1
    ACATGAGACAGACAGACCCCCAGAGACAGACCCCTAGACACAGAGAGAG
    TATGCAGGACAGGGTTTTTGCCCAGGGTGGCAGTATG

A FASTA file can contain more than one sequence. If a FASTA file
contains many sequences, then for each sequence it will have a
header line starting with '>' followed by the sequence itself.

::

    >mysequence1
    ACATGAGACAGACAGACCCCCAGAGACAGACCCCTAGACACAGAGAGAG
    TATGCAGGACAGGGTTTTTGCCCAGGGTGGCAGTATG
    >mysequence2
    AGGATTGAGGTATGGGTATGTTCCCGATTGAGTAGCCAGTATGAGCCAG
    AGTTTTTTACAAGTATTTTTCCCAGTAGCCAGAGAGAGAGTCACCCAGT
    ACAGAGAGC

NCBI Sequence Format (NCBI Format)
----------------------------------

As mentioned above, for each sequence the NCBI database stores some
extra information such as the species that it came from,
publications describing the sequence, etc. This information is
stored in the NCBI entry or NCBI record for the sequence. The NCBI
entry for a sequence can be viewed by searching the NCBI database
for the accession number for that sequence. The NCBI entries for
sequences are stored in a particular format, known as NCBI format.

To view the NCBI entry for the DEN-1 Dengue virus (which has
accession NC\_001477), follow these steps:


#. Go to the NCBI website
   (`www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_).
#. Search for the accession number.
#. On the results page, if your sequence corresponds to a
   nucleotide (DNA or RNA) sequence, you should see a hit in the
   Nucleotide database, and you should click on the word 'Nucleotide'
   to view the NCBI entry for the hit. Likewise, if your sequence
   corresponds to a protein sequence, you should see a hit in the
   Protein database, and you should click on the word 'Protein' to
   view the NCBI entry for the hit.
#. After you click on 'Nucleotide' or 'Protein' in the previous
   step, the NCBI entry for the accession will appear.

For example, the NCBI entry for the DEN-1 Dengue virus genome sequence
(NCBI accession NC\_001477) looks like this:

|image2|

The NCBI entry for an accession contains a lot of information about
the sequence, such as papers describing it, features in the
sequence, etc. The 'DEFINITION' field gives a short description for
the sequence. The 'ORGANISM' field in the NCBI entry identifies the
species that the sequence came from. The 'REFERENCE' field contains
scientific publications describing the sequence. The 'FEATURES'
field contains information about the location of features of
interest inside the sequence, such as regulatory sequences or genes
that lie inside the sequence. The 'ORIGIN' field gives the
sequence itself.

RefSeq
------

When carrying out searches of the NCBI database, it is important to
bear in mind that the database may contain redundant sequences for
the same gene that were sequenced by different laboratories. Some
of these sequences may be better quality than others.

There are also many different types of nucleotide sequences and
protein sequences in the NCBI database. With respect to nucleotide
sequences, some many be entire genomic DNA sequences, some may be
mRNAs, and some may be lower quality sequences such as expressed
sequence tags (ESTs, which are derived from parts of mRNAs), or DNA
sequences of contigs from genome projects. Furthermore, some
sequences may be manually curated so that the associated entries
contain extra information, but the majority of sequences are
uncurated.

As mentioned above, the NCBI database often contains redundant
information for a gene (because many different labs have sequenced
the gene, and submitted their sequences to the NCBI database), some
of which may be low quality. As a result, NCBI has made a special
database called RefSeq (reference sequence database), which is a
subset of the NCBI database. The data in RefSeq is manually
curated, is high quality sequence data, and is non-redundant; this means
that each gene (or splice-form of a gene, in the case of eukaryotes),
protein, or genome sequence is only represented once. 

The data in RefSeq is curated and is of much higher quality than the rest of the NCBI Sequence
Database. However, unfortunately, because of the high level of
manual curation required, RefSeq does not cover all species, and is
not comprehensive for the species that are covered so far.

You can easily tell that a sequence comes from RefSeq because its
accession number starts with particular sequence of letters. That
is, accessions of RefSeq sequences corresponding to protein records usually start with
'NP\_', and accessions of RefSeq curated complete genome sequences usually start with
'NC\_' or 'NS\_'.

Querying the NCBI Database
--------------------------

You may need to interrogate the NCBI Database
to find particular sequences or a set of sequences matching given
criteria, such as:
  
-  The sequences published in *Nature* **460**:352-358
-  All sequences from *Chlamydia trachomatis*
-  Sequences submitted by Matthew Berriman
-  Flagellin or fibrinogen sequences
-  The glutamine synthetase gene from *Mycobacteriuma leprae*
-  The upstream control region of the *Mycobacterium leprae dnaA* gene
-  The sequence of the *Mycobacterium leprae* DnaA protein
-  The genome sequence of *Trypanosoma cruzi*
-  All human nucleotide sequences associated with malaria

There are two main ways that you can query the NCBI database to find these
sets of sequences. The first possibility is to carry out searches on the 
`NCBI website <http://www.ncbi.nlm.nih.gov>`_.
The second possiblity is to carry out searches from R. 

In the examples below, we will show how to use both methods to carry out
queries on the NCBI database. In general, the two methods should give the
same result, but in some cases they do not, for various reasons, as shall be explained below.

Querying via the NCBI Website
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are carrying out searches on the `NCBI website <http://www.ncbi.nlm.nih.gov>`_, 
to narrow down your searches to specific types of sequences or to specific organisms, 
you will need to use "search tags".

For example, the search tags "[PROP]" and "[ORGN]" 
let you restrict your search to a specific subset of the
NCBI Sequence Database, or to sequences from a particular taxon,
respectively. We will explain how to use these search tags below.

Other useful NCBI search tags, which shall be illustrated in examples below, are:

-  "[JOUR]": to restrict your search to sequences described in a
   paper published in a particular journal
-  "[VOL]": to restrict your search to sequences described in a
   paper published in a particular volume of a journal
-  "[PAGE]": to restrict your search to sequences described in a
   paper with a particular start-page in a journal
-  "[AU]": to restrict your search to sequences submitted to the
   NCBI Database by a particular person, or described in a journal
   paper by a particular person. The person's name should be in the
   form: surname first-initial (eg. Bloggs J[AU])
-  "[ORGN]": to restrict your search to sequences from a particular
   species or taxon (eg. *Mycobacterium leprae* or *Mycobacterium* or Bacteria or
   Archaea)
-  "[PROP]": to restrict your search to a particular subset of the
   NCBI database (eg. "srcdb\_refseq[PROP]" restricts your search to
   RefSeq) or to a particular type of molecule (eg. "biomol
   mrna[PROP]" restrict your search to mRNA sequences).

Querying via R
^^^^^^^^^^^^^^

Instead of carrying out searches of the NCBI database on the NCBI website, you can
carry out searches directly from R by using the SeqinR R package.

It is possible to use the SeqinR R library to retrieve sequences from these databases.
The SeqinR library was written by the group that created the ACNUC database in Lyon, France
(http://pbil.univ-lyon1.fr/databases/acnuc/acnuc.html).
The ACNUC database is a database that contains most of the data from the NCBI Sequence Database,
as well as data from other sequence databases such as UniProt and Ensembl. 

An advantage of the ACNUC database is that it brings together data from various different sources, and makes
it easy to search, for example, by using the SeqinR R library.

As will be explained below, the ACNUC database is organised into various different ACNUC (sub)-databases,
which contain different parts of the NCBI database, and when you want to search the NCBI database
via R, you will need to specify which ACNUC sub-database the NCBI data that you want to query is stored in.

To obtain a full list of the ACNUC sub-databases that you can access using SeqinR, you
can use the "choosebank()" function from SeqinR:

::

    > library("seqinr") # Load the SeqinR R package
    > choosebank()      # List all the sub-databases in ACNUC
      [1] "genbank"       "embl"          "emblwgs"       "swissprot"    
      [5] "ensembl"       "hogenom"       "hogenomdna"    "hovergendna"  
      [9] "hovergen"      "hogenom4"      "hogenom4dna"   "homolens"     
      [13] "homolensdna"   "hobacnucl"     "hobacprot"     "phever2"      
      [17] "phever2dna"    "refseq"        "nrsub"         "greviews"     
      [21] "bacterial"     "protozoan"     "ensbacteria"   "ensprotists"  
      [25] "ensfungi"      "ensmetazoa"    "ensplants"     "mito"         
      [29] "polymorphix"   "emglib"        "taxobacgen"    "refseqViruses"

Two of the most important sub-databases in ACNUC which can be searched from R are:

- genbank: this contains DNA and RNA sequences from the NCBI Sequence Database, except for certain
classes of sequences (eg. draft genome sequence data from genome sequencing projects)
- refseq: this contains DNA and RNA sequences from `Refseq <./chapter3.html#refseq>`_, the curated part of the NCBI Sequence Database
- refseqViruses: this contains DNA, RNA and proteins sequences from viruses from RefSeq 

You can find more information about what each of these ACNUC databases contains by
looking at the `ACNUC website <http://pbil.univ-lyon1.fr/search/releases.php>`_. 

You can carry out complex queries using the "query()" function from
the SeqinR library. If you look at the help page for the query() function (by
typing "help(query)", you will see that it allows you to specify criteria that you
require the sequences to fulfill. 

For example, to search for a sequence with a particular NCBI accession, you can use the "AC=" argument in "query()".
The "query()" function will then search for 
sequences in the NCBI Sequence Database that match your criteria. 

In the examples below, we will explain how to carry out searches of the NCBI database
both by searching the ACNUC database via R, and by going directly to the NCBI website to carry out the search.

Example: finding the sequences published in *Nature* **460**:352-358
--------------------------------------------------------------------------------

.. rubric:: Method 1: searching via the NCBI website
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To find the sequences published in *Nature* **460**:352-358, one method is to 
go to the `NCBI website <http://www.ncbi.nlm.nih.gov>`_ and type in the search 
box on the top: "Nature"[JOUR] AND 460[VOL] AND 352[PAGE]

|image3|

Here [JOUR] specifies the journal name, [VOL] the volume of the journal the paper is in, and [PAGE] the page number.

This should bring up a results page with "50890" beside the word "Nucleotide", and "1" beside the word
"Genome", and "25701" beside the word "Protein", indicating that there were 50890 hits to sequence records in the Nucleotide database, 
which contains DNA and RNA sequences, and 1 hit to the Genome database, which contains genome sequences, and 25701
hits to the Protein database, which contains protein sequences:

|image4|

If you click on the word "Nucleotide", it will bring up a webpage with a list of links to the NCBI sequence 
records for those 50890 hits. The 50890 hits are all contigs from the schistosome worm *Schistosoma mansoni*.

Likewise, if you click on the word "Protein", it will bring up a webpage with a list of links to the NCBI
sequence records for the 25701 hits, and you will see that the hits are all predicted proteins for *Schistosoma
mansoni*.

If you click on the word "Genome", it will bring you to the NCBI record for the *Schistosoma mansoni* genome
sequence, which has NCBI accession NS\_00200. Note that the accession starts with "NS\_", which indicates that
it is a RefSeq accession. 

Therefore, in *Nature* volume 460, page 352, the *Schistosoma mansoni* genome sequence was published, along
with all the DNA sequence contigs that were sequenced for the genome project, and all the predicted proteins
for the gene predictions made in the genome sequence. You can view the original paper on the *Nature* website
at `http://www.nature.com/nature/journal/v460/n7253/abs/nature08160.html <http://www.nature.com/nature/journal/v460/n7253/abs/nature08160.html>`_.

Note: *Schistmosoma mansoni* is a parasitic worm that is responsible for causing 
`schistosomiasis <http://apps.who.int/tdr/svc/diseases/schistosomiasis>`_, 
which is classified by the WHO as a neglected tropical disease.

.. rubric:: Method 2: searching via R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To find the sequences published in *Nature* **460**:352-358, a second method is to use the SeqinR
R package to search the ACNUC databases (which contain the NCBI sequence data) from R. 

If you look at the help page the "query()" function, you see that you can query for sequences 
published in a particular paper using R=refcode, specifying the reference as refcode 
such as in jcode/volume/page (e.g., JMB/13/5432 or R=Nature/396/133). For the
paper *Nature* **460**:352-358, we would need to use the refcode 'R=Nature/460/352'.

As explained above, the ACNUC database contains the NCBI sequence data organised into several
sub-databases, and you can view the list of those sub-databases by using the "choosebank()"
function from the SeqinR package. When you want to use "query()" to carry out a particular 
sub-database (eg. "genbank", which contains DNA and RNA sequences from the NCBI Sequence Database), you
need to first specify the database that you want to search by using the "choosebank()" function,
for example:

::

    > choosebank("genbank") # Specify that we want to search the 'genbank' ACNUC sub-database

We can then search the 'genbank' database for sequences that match a specific set of criteria
by using the "query()" function. For example, to search for sequences that were published in
*Nature* **460**:352-358, we type:

::

    > query('naturepaper', 'R=Nature/460/352')

The line above tells R that we want to store the results of the query in an R list variable called
*naturepaper*. Remember that a list is an R object that is like a vector, but can contain elements
that are numeric and/or contain characters. In this case, the list *naturepaper* contains information
on the NCBI records that match the query (ie. information on the NCBI records that contain sequences
published in Nature* **460**:352-358). 

If you look at the help page for "query()", the details of the arguments are given under the heading "Arguments",
and the details of the results (outputs) are given under the heading "Value". If you read this now, you
will see that it tells us that the result of the "query()" function is a list with six different named
elements, named "call", "name", "nelem", "typelist", "req", and "socket". The content of each of these
six named elements is explained, for example, the "nelem" element contains the number of sequences
that match the query, and the "req" element contains their accession numbers.

In our example, the list object *naturepaper* is an output of the "query()" function, and
so has each of these six named elements, as we can find out by using the "attributes()" function, 
and looking at the named elements listed under the heading "$names":

::

    > attributes(naturepaper)
      $names
      [1] "call"     "name"     "nelem"    "typelist" "req"      "socket"  
      $class
      [1] "qaw"

As explained in the `brief introduction to R <./installr.html#a-brief-introduction-to-r>`_, we can retrieve 
the value for each of the named elements in the list *naturepaper* by using "$", followed by the element's name, 
for example, to get the value of the element named "nelem" in the list *naturepaper*, we type:

::

    > naturepaper$nelem
      [1] 19022

This tells us that there were 19022 sequences in the 'genbank' ACNUC database that matched the query.
The 'genbank' ACNUC database contains DNA or RNA sequences from the NCBI Nucleotide database. 
Why don't we get the same number of sequences as found by carrying out the search on the NCBI website
(where we found 50890 hits to the NCBI Nucleotide database)? The reason is that the ACNUC 'genbank'
database does not contain all the sequences in the NCBI Nucleotide database, for example, it does
not contain sequences that are in RefSeq or many short DNA sequences from sequencing projects.

To obtain the accession numbers of the first five of the 19022 sequences, we can type:

::

    > accessions <- naturepaper$req
    > accessions[1:5]
      [[1]]
           name     length      frame     ncbicg 
      "FN357292"  "4179495"        "0"        "1" 
      [[2]]
           name     length      frame     ncbicg 
      "FN357293"  "2211188"        "0"        "1" 
      [[3]]
           name     length      frame     ncbicg 
      "FN357294"  "1818661"        "0"        "1" 
      [[4]]
           name     length      frame     ncbicg 
      "FN357295"  "2218116"        "0"        "1" 
      [[5]]
           name     length      frame     ncbicg 
      "FN357296"  "3831198"        "0"        "1" 

This tells us that the NCBI accessions of the first five sequences (of the 19022
DNA or RNA sequences found that were published in *Nature* **460**:352-358) are FN357292,
FN357293, FN357294, FN357295, and FN357296. 

Example: finding all human nucleotide sequences associated with malaria
-----------------------------------------------------------------------

Say for example that you want to find all high-quality 
nucleotide sequences associated with malaria. Here were are not looking
for sequences from the malaria organism itself, but sequences of human
genes and human proteins that somehow interact with or respond to the malaria organism.

To find all nucleotide sequences associated with malaria, follow these
steps:

#. Go to the NCBI website
   (`www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_).
#. As you want to search for nucleotide sequences, select
   'Nucleotide' from the drop-down list above the search box at the
   top of the NCBI homepage.
#. Type **malaria** in the search box. (Note that if you are searching for
   a phrase such as 'colon cancer', you would need to 
   include the inverted commas, ie. type **"colon cancer"** and not
   **colon cancer**. This is because if you type just
   **colon cancer**, the search will be for records that contain the
   words 'colon' or 'cancer' (not necessarily both words), while you
   want records that contain the phrase 'colon cancer'.) Press 'Search'.
#. The search results will include all nucleotide sequences for
   which the phrase 'malaria' appears somewhere in their NCBI
   records. The phrase may appear in the 'DEFINITION' field of the
   NCBI record (which gives a short description), in the title of a
   journal article on the nucleotide sequence, or elsewhere in the
   NCBI record.

The search above should have identified thousands of sequences from
many different species. Some of these may be of low quality. To
limit your search to high quality sequences, you may decide to
restrict your search to RefSeq sequences. You can do this using
NCBI search tags. NCBI search tags allow you to limit your restrict
your search to a specific data set, such as the RefSeq data set. It
also allows us to limit searches to retrieve records with certain
attributes, such as molecule type (eg. mRNAs) or species.

The NCBI search tag "[PROP]" allows you to restrict your search to
sequences form a particular subset of the NCBI Sequence Database,
such as RefSeq. To use NCBI search tags to restrict your search to
nucleotide sequences from RefSeq that are associated with malaria, follow these steps:

#. Go to the NCBI website, and select 'Nucleotide' from the
   drop-down list above the search box.
#. In the search box, type
   **malaria AND srcdb\_refseq[PROP]**, and press 'Search'.

This should give you all RefSeq nucleotide sequences for which the phrase
malaria appears somehwere in the NCBI record.

Note that you should find fewer sequences than when you just
searched for **malaria**, but these should be higher quality
sequences (since they are RefSeq sequences), 
and their NCBI entries will contain manually curated
information about the sequences (eg. details of publications about
the sequences and features in them).

The search above should have identified RefSeq sequences from
several species (eg. malaria itself, human, mouse, etc.) that are associated with
malaria (or more precisely, where the word 'malaria'
appears somewhere in the NCBI records). 
What if you are only interested in human sequences
associated with malaria?

One way to solve this problem is to use NCBI search tags to
restrict your search to human sequences. The "[ORGN]" search tag
allows you to restrict your search to sequences from a particular
species (eg. *Mycobacteriuma leprae*, the bacterium that causes
leprosy, or set of species (eg. Bacteria). To use NCBI search tags to retrieve human RefSeq
sequences associated with malaria, follow these steps:

#. Go to the NCBI website, and select 'Nucleotide' from the
   drop-down list above the search box.
#. In the search box, type
   **malaria AND srcdb\_refseq[PROP] AND "Homo sapiens"[ORGN]**,
   and press 'Search'.

This will give you a list of all human nucleotide sequences from
RefSeq that are associated with malaria (or more precisely, all
the human nucleotide sequences from Refseq for which the word 'malaria'
appears somewhere in the NCBI record).

Finding the genome sequence for a particular species
----------------------------------------------------

Microbial genomes are generally smaller than eukaryotic genomes
(*Escherichia coli* has about 5 million base pair in its genome,
while the human genome is about 3 billion base pairs). Because they
are considerably less expensive to sequence, many microbial genome
sequencing projects have been completed.

If you don't know the accession number for a genome sequence (eg.
for *Mycobacterium leprae*, the bacterium that causes leprosy), how can you find it out? One way to
do this is to look at the NCBI Genome website, which lists all
fully sequenced genomes and gives the accession numbers for the
corresponding DNA sequences.

If you didn't know the accession number for the
*Mycobacterium leprae* genome, you could find it on the NCBI
Genome website by following these steps:

#. Go to the NCBI Genome website
   (`http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome <http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome>`_)
#. On the homepage of the NCBI Genome website, it gives links to the
   major subdivisions of the Genome database, which include
   Eukaryota, Prokaryota (Bacteria and Archaea), and Viruses.
   Click on 'Prokaryota', since
   *Mycobacterium leprae* is a bacterium. This will bring up a list
   of all fully sequenced bacterial genomes, with the corresponding
   accession numbers. Note that more than one genome (from various
   strains) may have been sequenced for a particular species.
#. Use 'Find' in the 'Edit' menu of your web browser to search for
   'Mycobacterium leprae' on the webpage. You should find that the
   genomes of several different *M. leprae* strains have been
   sequenced. One of these is *M. leprae* TN, which has
   accession number NC\_002677.

The list of sequenced genomes on the NCBI Genomes website is not a
definitive list; that is, some sequenced genomes may be missing
from this list. If you want to find out whether a particular genome
has been sequenced, but you don't find it NCBI Genomes website's
list, you should search for it by following these steps:

#. Go to the NCBI website
   (`www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_).
#. Select 'Genome' from the drop-down list above the search box.
#. Type the name of the species you are interested in in the search
   box (eg. **"Mycobacterium leprae"[ORGN]**). Press 'Search'.

Note that you could also have found the *Mycobacterium leprae*
genome sequence by searching the NCBI Nucleotide database, as the
NCBI Genome database is just a subset of the NCBI Nucleotide
database.

How many genomes have been sequenced, or are being sequenced now?
-----------------------------------------------------------------

On the NCBI Genome website
(`http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome <http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome>`_),
the front page gives a link to a list of all sequenced genomes in the
groups Eukaryota, Prokaryota (Bacteria and Archaea) and Viruses.
If you click on one of these links (eg. Prokaryota), at the top of the
page it will give the number of sequenced genomes in that group (eg. number of sequenced
prokaryotic genomes). For example, in this screenshot (from January 2011), we see that there
were 1409 complete prokaryotic genomes (94 archaeal, 1315 bacterial):

|image1| 

Another useful website that lists genome sequencing projects is the
Genomes OnLine Database (GOLD), which lists genomes that have been
completely sequenced, or are currently being sequenced. To find the
number of complete or ongoing bacterial sequencing projects, follow
these steps:

#. Go to the GOLD website
   (`http://genomesonline.org/ <http://genomesonline.org/>`_).
#. Click on the yellow 'Enter GOLD' button in the centre of the
   webpage. On the subsequent page, it will give the number of ongoing
   bacterial, archaeal and eukaryotic genome sequencing projects.
#. Click on the 'Bacterial Ongoing' link to see the list of
   ongoing bacterial genome sequencing projects. By default, just the
   first 100 projects are listed, and the rest are listed on subsequent pages.
   In one of the columns
   of the page, this gives the university or institute that the genome
   was sequenced in. Other columns give the taxonomic information for
   the organism, and links to the sequence data.
#. Find the number of published genome sequencing projects. Go back
   one page, to the page with the 'Bacterial Ongoing' link. 
   You will see that this page also lists the number of complete published
   genomes. To see a list of these genomes, click on 'Complete Published'.
   This will bring up a page that gives the number of published
   genomes at the top of the page. In one column of the page, this
   gives the university or institute that the genome was sequenced
   in.

As explained above, it is possible to identify genome sequence data
in the NCBI Genome database. The GOLD database also gives some
information about ongoing genome projects. Often, the GOLD database
lists some ongoing projects that are not yet present in the NCBI
Genome Database, because the sequence data has not yet been
submitted to the NCBI Database. If you are interested in finding
out how many genomes have been sequenced or are currently being
sequenced for a particular species (eg. *Mycobacterium leprae*), it
is a good idea to look at both the NCBI Genome database and at
GOLD.

Summary
-------

In this chapter, you have learnt how to retrieve sequences from
the NCBI Sequence database, as well as to find out how many genomes
have been sequenced or are currently being sequenced for a
particular species.

Links and Further Reading
-------------------------

There is detailed information on how to search the NCBI database on
the NCBI Help website at
`http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=helpentrez?part=EntrezHelp <http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=helpentrez%26part=EntrezHelp>`_.

There is more information about the GOLD database in the paper
describing GOLD by Liolios *et al*, which is available at
`http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808860/?tool=pubmed <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808860/?tool=pubmed>`_.

Acknowledgements
----------------

Thank you to Noel O'Boyle for helping in using Sphinx, `http://sphinx.pocoo.org <http://sphinx.pocoo.org>`_, to create
this document, and github, `https://github.com/ <https://github.com/>`_, to store different versions of the document
as I was writing it, and readthedocs, `http://readthedocs.org/ <http://readthedocs.org/>`_, to build and distribute
this document.

Thank you to Andrew Lloyd and David Lynn, who generously shared their practical on sequence databases 
with me, which inspired many of the examples in this practical. 

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.

Exercises
---------

Answer the following questions. For each question, please record
your answer, and what you did/typed to get this answer.

Model answers to the exercises are given in the chapter entitled
`Answers to the exercises on Sequence Databases <./chapter3_answers.html>`_.

Q1. What information about the DEN-1 Dengue virus sequence (NCBI accession NC\_001477) can you obtain from its annotations in the NCBI Sequence Database? 
    What does it say in the DEFINITION and ORGANISM fields of its NCBI
    record?
Q2. What were the nucleotide sequences published in *Nature* volume 460, page 352?
    What are their accession numbers in the NCBI Sequence Database?

Q3. How many nucleotide sequences are there from the bacterium *Chlamydia trachomatis* in the NCBI Sequence Database? 
    Remember to type **"Chlamydia trachomatis"** including the inverted commas.

Q4. How many nucleotide sequences are there from the bacterium *Chlamydia trachomatis* in the *RefSeq* part of the 
    NCBI Sequence Database? 

Q5. How many nucleotide sequences were submitted to NCBI by Matthew Berriman?
    Note that the name of the person who submitted a sequence is stored
    in the author field of the NCBI record, as is the name of people
    who published papers on the sequence. There may be more than one
    author fields in the NCBI record for a sequence, corresponding to
    the person who submitted the sequence and/or people who published
    papers on the sequence.

Q6. How many nucleotide sequences from nematode worms are there in the NCBI Database? 

Q7. How many nucleotide sequences for collagen genes from nematode worms are there in the NCBI Database? 
    Hint: look at the examples above for malaria-related genes.

Q8. How many *mRNA sequences* for collagen genes from nematode worms are there in the NCBI Database? 
    Hint: look at the notes about the "[PROP]" search tag above.

Q9. How many *protein sequences* for collagen proteins from nematode worms are there in the NCBI database? 

Q10. What is the accession number for the *Trypanosoma cruzi* genome in NCBI? 
    Do you see genome sequences for more than one strain of *Trypanosoma cruzi*?

Q11. How many fully sequenced nematode worm species are represented in the NCBI Genome database? 

xxx


Q12. How many ongoing genome sequencing projects are there for Bacteria, Archaea, and Eukarotes, respectively, in the GOLD database? Q13. Are there any genome sequencing projects ongoing at University College Cork, acccording to the GOLD database? 
    Hint: Use the 'Find' option in the 'Edit' menu of your web browser
    to search for 'Cork' in the GOLD database's webpage listing ongoing
    genome sequencing projects.
Q14. How many genome sequences are there for *Lactobacillus salivarius* in the NCBI Genomes database? 
    Why are there more than one?
Q15. How many complete or ongoing genome sequencing projects for *Lactobacillus salivarius* are listed in GOLD? 
    Does GOLD or NCBI Genomes have more sequencing projects for this
    species? If not, can you suggest an explanation why?

.. |image1| image:: ../_static/P3_image1.png
            :width: 900
.. |image2| image:: ../_static/P3_image2.png
            :width: 500
.. |image3| image:: ../_static/P3_image3.png
            :width: 500
.. |image4| image:: ../_static/P3_image4.png
            :width: 500
