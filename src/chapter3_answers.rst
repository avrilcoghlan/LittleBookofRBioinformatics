Answers to the exercises on Sequence Databases
==============================================   

Q1. *What information about the rabies virus sequence (NCBI accession NC\_001542) can you obtain from its annotations in the NCBI Sequence Database?*

To do this, you need to go to the `www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_ website 
and type the rabies virus genome sequence accession (NC\_001542) in the search box, and press 'Search'. 

On the search results page, you should see '1' beside the word 'Nucleotide', meaning that there was one hit to a sequence record in the NCBI Nucleotide database, which contains DNA and RNA sequences. If you click on the word 'Nucleotide', it will bring you to the sequence record, which should be the NCBI sequence record for the rabies virus' genome (ie. for accession NC\_001542):

|image7|

On the webpage (above), you can see the DEFINITION, ORGANISM and REFERENCE fields of the NCBI record:

DEFINITION: Rabies virus, complete genome.

ORGANISM: Rabies virus 

REFERENCE: There are several papers (the first is):
AUTHORS: Tordo,N., Poch,O., Ermine,A., Keith,G. and Rougeon,F.

TITLE: Completion of the rabies virus genome sequence determination: highly conserved domains among the L (polymerase) proteins of unsegmented negative-strand RNA viruses

JOURNAL: Virology 165 (2), 565-576 (1988)

There are also some other references, for papers published about the rabies virus genome sequence. 

An alternative way of retrieving the annotations for the rabies virus sequence is to use the SeqinR R package.
As the rabies virus is a virus, its genome sequence should be in the "refseqViruses" ACNUC sub-database.
Therefore, we can perform the following query to retrieve the annotations for the rabies virus
genome sequence (accession NC\_001542):

::

    > library("seqinr")                                 # load the SeqinR R library
    > choosebank("refseqViruses")                       # select the ACNUC sub-database to be searched
    > query("rabies", "AC=NC_001542")                   # specify the query
    > annots <- getAnnot(rabies$req[[1]])               # retrieve the annotations
    > annots[1:20]                                      # print out the first 20 lines of the annotations
      [1] "LOCUS       NC_001542              11932 bp ss-RNA     linear   VRL 08-DEC-2008"
      [2] "DEFINITION  Rabies virus, complete genome."                                     
      [3] "ACCESSION   NC_001542"                                                          
      [4] "VERSION     NC_001542.1  GI:9627197"                                            
      [5] "DBLINK      Project: 15144"                                                     
      [6] "KEYWORDS    ."                                                                  
      [7] "SOURCE      Rabies virus"                                                       
      [8] "  ORGANISM  Rabies virus"                                                       
      [9] "            Viruses; ssRNA negative-strand viruses; Mononegavirales;"           
      [10] "            Rhabdoviridae; Lyssavirus."                                         
      [11] "REFERENCE   1  (bases 5388 to 11932)"                                           
      [12] "  AUTHORS   Tordo,N., Poch,O., Ermine,A., Keith,G. and Rougeon,F."              
      [13] "  TITLE     Completion of the rabies virus genome sequence determination:"      
      [14] "            highly conserved domains among the L (polymerase) proteins of"      
      [15] "            unsegmented negative-strand RNA viruses"                            
      [16] "  JOURNAL   Virology 165 (2), 565-576 (1988)"                                   
      [17] "   PUBMED   3407152"                                                            
      [18] "REFERENCE   2  (bases 1 to 5500)"                                               
      [19] "  AUTHORS   Tordo,N., Poch,O., Ermine,A., Keith,G. and Rougeon,F."              
      [20] "  TITLE     Walking along the rabies genome: is the large G-L intergenic region"
    > closebank()

Q2. *How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the NCBI Sequence Database?*

To answer this, you need to go to `www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_, 
and select "Nucleotide" from the drop-down list at the top 
of the webpage, as you want to search for nucleotide (DNA or RNA) sequences.

Then in the search box, type "Chlamydia trachomatis"[ORGN] and press 'Search':

|image8|

Here [ORGN] specifies the organism you are interested in, that is, the species name in Latin.

The results page should give you a list of the hits to sequence records in the NCBI Nucleotide database: 

|image9|

It will say "Found 35577 nucleotide sequences.   Nucleotide (35429)   GSS (148)". 
This means that 35,577 sequences were found, of which 35429 are DNA or RNA sequences, and 
148 are DNA sequences from the Genome Sequence Surveys (GSS), that is, from 
genome sequencing projects [as of 15-Jun-2011]. Note that there are new sequences 
being added to the database continuously, so if you check this again in a couple of months, you will 
probably find a higher number of sequences (eg. 36,000 sequences).

Note: if you just go to the www.ncbi.nlm.nih.gov database, and search for "Chlamydia trachomatis"[ORGN] 
(without choosing "Nucleotide" from the drop-down list), you will see 35429 hits to the Nucleotide 
database and 148 to the GSS (Genome Sequence Survey) database:

|image10|

Note also that if you search for "Chlamydia trachomatis", without using [ORGN] to specify the organism, 
you will get 56032 hits to the Nucleotide database and 149 to the GSS database, but some of these might 
not be *Chlamydia trachomatis* sequences - some could be sequences from other species for which the NCBI sequence 
record contains the phrase "Chlamydia trachomatis" somewhere.

An alternative way to search for nucleotide sequences from the bacterium *Chlamydia trachomatis* is to
use the SeqinR package. We want to find nucleotide sequences, so the correct ACNUC sub-database to search
is the "genbank" sub-database. Thus, we can carry out our search by typing:

::

    > library("seqinr")                                 # load the SeqinR R library
    > choosebank("genbank")                             # select the ACNUC sub-database to be searched
    > query("Ctrachomatis", "SP=Chlamydia trachomatis") # specify the query
    > Ctrachomatis$nelem                                # print out the number of matching sequences
      [1] 35471
    > closebank()

We find 35,471 nucleotide sequences from *Chlamydia trachomatis*. We do not get exactly the same number
of sequences as we got when we searched via the NCBI website (35,577 sequences), but the numbers are very close.
The likely reasons for the differences could be that the ACNUC "genbank" sub-database excludes some sequences from
whole genome sequencing projects from the NCBI Nucleotide database, and in addition, the ACNUC databases
are updated very regularly, but may be missing a few sequences that were added to the NCBI database
in the last day or two.

Q3. *How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the RefSeq part of the NCBI Sequence Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov and select "Nucleotide" from the drop-down list 
at the top of the webpage, as you want to search for nucleotide sequences.

Then in the search box, type "Chlamydia trachomatis"[ORGN] AND srcdb_refseq[PROP] and press 'Search':

|image11|

Here [ORGN] specifies the organism, and [PROP] specifies a property of the sequences (in this case that 
they belong to the RefSeq subsection of the NCBI database).

At the top of the results page, it should say "Results: 1 to 20 of 29 sequences", so there were
29 matching sequences [as of 15-Jun-2011]. 
As for Q2, if you try this again in a couple of months, the number will probably be higher, due to extra 
sequences added to the database. 

Note that the sequences in Q2 are all *Chlamydia trachomatis* DNA and RNA sequences in the NCBI database. 
The sequences in Q3 gives the *Chlamydia trachomatis* DNA and RNA sequences in the RefSeq part of the NCBI 
database, which is a subsection of the database for high-quality manually-curated data. 

The number of sequences in RefSeq is much fewer than the total number of *C. trachomatis* sequences, 
partly because low quality sequences are never added to RefSeq, but also because RefSeq curators have 
probably not had time to add all high-quality sequences to RefSeq (this is a time-consuming process, 
as the curators add additional information to the NCBI Sequence records in RefSeq, such as references to 
papers that discuss a particular sequence). 

An alternative way to search for nucleotide sequences from the bacterium *Chlamydia trachomatis* in RefSeq
use the SeqinR package. We want to find RefSeq sequences, so the correct ACNUC sub-database to search
is the "refseq" sub-database. Thus, we can carry out our search by typing:

::

    > library("seqinr")                                  # load the SeqinR R library
    > choosebank("refseq")                               # select the ACNUC sub-database to be searched
    > query("Ctrachomatis2", "SP=Chlamydia trachomatis") # specify the query
    > Ctrachomatis2$nelem                                # print out the number of matching sequences
      [1] 1
    > closebank()

We find 1 RefSeq sequence from *Chlamydia trachomatis*. We do not get exactly the same number
of sequences as we got when we searched via the NCBI website (29 sequences). This is because the
29 sequences found via the NCBI website include whole genome sequences, but the whole genome sequences
from bacteria are stored in the ACNUC "bacterial" sub-database, and so are not in the ACNUC "refseq" 
sub-database.

Q4. *How many nucleotide sequences were submitted to NCBI by Matthew Berriman?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Nucleotide" from the drop-down list, 
as you want to search for nucleotide sequences.

Then in the search box, type "Berriman M"[AU] and press 'Search'.

Here [AU] specifies the name of the person who either submitted the sequence to the NCBI database, 
or wrote a paper describing the sequence. 

On the results page, it should say at the top: "Found 460052 nucleotide sequences.   Nucleotide (250328)   EST (121075)   GSS (88649)". This means that 460052 DNA/RNA sequences were either submitted to the NCBI database by someone called M. Berriman, or were described in a paper by someone called M. Berriman. Of these, 250328 were DNA/RNA sequences, 121075 were EST sequences (part of mRNAs), and 88649 were DNA sequences from genome sequencing projects (GSS or Genome Sequence Survey sequences).

Note that unfortunately the NCBI website does not allow us to search for "Berriman Matthew"[AU] so we cannot be sure 
that all of these sequences were submitted by Matthew Berriman. 

Q5. *How many nucleotide sequences from nematode worms are there in the NCBI Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov and select "Nucleotide" from the drop-down list, 
as you want to search for nucleotide sequences.

Then in the search box, type Nematoda[ORGN] and press 'Search'.

Here [ORGN] specifies the group of species that you want to search for sequences from. In Q4, [ORGN] was used to specify 
the name of one organism (*Chlamydia trachomatis*). However, you can also use [ORGN] to specify the name of a group of 
organisms, for example, Fungi[ORGN] would search for fungal sequences or Mammalia[ORGN] would search for mammalian 
sequences. The name of the group of species that you want to search for must be given in Latin, so to search for sequences
from nematode worms we use the Latin name Nematoda.

The search page should say at the top 'Found 2202458 nucleotide sequences.   Nucleotide (378255)   EST (1140454)   GSS (683749)' [as of 19-Feb-2011]. This means that 2,202,458 DNA or RNA sequences were found from nematode worm species in the database, of
which 378,255 are DNA/RNA sequences, 1,140,454 are ESTs, and 683,749 sequences are DNA sequences from genome sequencing
projects. These sequences are probably from a wide range of nematode worm species, including the model nematode worm
*Caenorhabditis elegans*.

Q6. *How many nucleotide sequences for collagen genes from nematode worms are there in the NCBI Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov and select "Nucleotide" from the drop-down list, 
as you want to search for nucleotide sequences.

Then in the search box, type Nematoda[ORGN] AND collagen.

Here [ORGN] specifies that we want sequences from nematode worms. The phrase "AND collagen" means that the word collagen 
must appear somewhere in the NCBI entries for those sequences, for example, in the sequence name, or in a description 
of the sequence, or in the title of a paper describing the sequence, etc.

On the results page, you should see 'Found 8341 nucleotide sequences.   Nucleotide (1546)   EST (6795)' [as of 19-Feb-2011].
This means that 8341 DNA or RNA sequences for collagen genes from nematode worms were found, of which 6795 are EST sequences
(parts of mRNAs). Note that these 8341 nucleotide sequences may not all necessarily be for collagen genes, as some of the
NCBI records found may be for other genes but contain the word 'collagen' somewhere in the NCBI record (for example, in
the title of a cited paper).

Q7. *How many mRNA sequences for collagen genes from nematode worms are there in the NCBI Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Nucleotide" from the drop-down sequences, as you want to search for nucleotide sequences (nucleotide sequences include DNA sequences and RNA sequences, such as mRNAs). 

Then in the search box, type Nematoda[ORGN] AND collagen AND "biomol mRNA"[PROP].

Here [ORGN] specifies the name of the group of species, collagen specifies that we want to find NCBI entries 
that include the word collagen, and [PROP] specifies a property of those sequences (that they are mRNAs, in this case).

The search page should say 'Found 7656 nucleotide sequences.   Nucleotide (861)   EST (6795)' [as of 19-Feb-2011].
This means that 7656 mRNA sequences were found that contain the word 'collagen' in the NCBI record. Of the
7656, 6795 are EST sequences (parts of mRNAs). 

Note that in Q7 we found 8341 nucleotide (DNA or RNA) sequences from nematode worms. In this question, we found out that 
only 7656 of those sequences are mRNA sequences. This means that the other (8341-7656=) 685 sequences must be DNA sequences, 
or other types of RNA sequences (not mRNAs) such as tRNAs or rRNAs.

Q8. *How many protein sequences for collagen proteins from nematode worms are there in the NCBI database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Protein" from the drop-down list, 
as you want to search for protein sequences.

Then type in the search box: Nematoda[ORGN] AND collagen and press 'Search'.

On the results page, you should see '1 to 20 of 1886'. This means that 1886 protein sequences from nematode
worms were found that include the word collagen in the NCBI sequence entries [as of 19-Feb-2011].

Q9. *What is the accession number for the Trypanosoma cruzi genome in NCBI?*

There are two ways that you can answer this.

The first method is to go to www.ncbi.nlm.nih.gov and select "Genome" from the drop-down list, 
as you want to search for genome sequences.

Then type in the search box: "Trypanosoma cruzi"[ORGN] and press 'Search'.

The results page says 'All:1', and lists just one NCBI record, the genome sequence for *Trypanosoma cruzi*
strain CL Brener, which has accession NZ\_AAHK00000000.

The second method of answering the question is to go to the NCBI Genomes webpage
http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome.

Click on the 'Eukaryota' link at the middle the page, as *Trypanosoma cruzi* is a eukaryotic species.

This will give you a complete list of all the eukaryotic genomes that have been sequenced.

Go to the 'Edit' menu of your web browser, and choose 'Find', and search for 'Trypanosoma cruzi'.  

You should find *Trypanosoma cruzi* strain CL Brener.
You will also find that there are several ongoing genome sequencing projects listed for other strains of
*Trypanosoma cruzi*: strains JR cl. 4, Sylvio X10/1, Y, and Esmeraldo Esmeraldo cl. 3.

The link 'GB' (in green) at the far right of the webpage gives a link to the NCBI record for the sequence.
In this case, the link for *Trypanosoma cruzi* strain CL Brener leads us to the NCBI record for accession
AAHK01000000. This is actually an accession for the *T. cruzi* strain CL Brener sequencing project, rather than
for the genome sequence itself. On the top right of the page, you will see a link "Genome", and if you click
on it, it will bring you to the NCBI accession NZ\_AAHK00000000, the genome sequence for *Trypanosoma cruzi* strain CL Brener.

Of the other *T. cruzi* strains listed, there is only a 'GB' link for one other strain, Sylvio X10/1.
Presumably there are no links for the other *Trypanosoma cruzi* strains, because the sequencing
projects are still in progress. If you click on the link for *Trypanosoma cruzi* strain Sylvio X10/1, it will bring you to the
NCBI record for accession ADWP01000000, the accession for the *T. cruzi* strain Sylvio X10/1 sequencing
project. At the top right of that page, there is no "Genome" link, which tells you that there is not yet
a genome assembly available for this strain. 

Note that the answer is slightly different for the answer from the first method above, which 
did not find the information on the genome projects for strains JR cl. 4, Sylvio X10/1, Y, and Esmeraldo Esmeraldo cl. 3,
because genome assemblies are not yet available for those strains.

Q10. *How many fully sequenced nematode worm species are represented in the NCBI Genome database?*

To answer this question, you need to go to the NCBI Genome webpage http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome. 

In the search box at the top of the page, type Nematoda[ORGN] to search for genome sequences from nematode   
worms, using the Latin name for the nematode worms. 

On the results page, you will see 'Items 1 - 20 of 62', indicating that 62 genome sequences from nematode worms
have been found. If you look down the page, you will see however that many of these are mitochondrial genome
sequences, rather than chromosomal genome sequences.

If you are just interested in chromosomal genome sequences, you can type 'Nematoda[ORGN] NOT mitochondrion' in the
search box, to search for non-mitochondrial sequences. This should give you 16 sequences, which are all chromosomal
genome sequences for nematode worms, including the species *Caenorhabditis elegans*, *Caenorhabditis remanei*,
*Caenorhabditis briggsae*, *Loa loa* (which causes subcutaneous filariasis), and 
*Brugia malayi* (which causes lymphatic filariasis). Thus, there are nematode genome sequences from five different
species that have been fully sequenced (as of 19-Feb-2011). Because nematode worms are multi-chromosomal species, 
there may be several chromosomal sequences for each species.

Note that when you search the NCBI Genome database at http://www.ncbi.nlm.nih.gov/sites/entrez?db=Genome, you will
find the NCBI records for completely sequenced genomes (completely sequenced nematode genomes in this case).

If you are interested in partially sequenced genomes, that is sequences from genome sequencing projects that are
still in progress, you can go to the NCBI Genome Projects website at http://www.ncbi.nlm.nih.gov/genomeprj. If you
search the NCBI Genome Projects database for Nematoda[ORGN], you will find that genome
sequencing projects for many other nematode species are ongoing, including for the species *Onchocerca volvulus*
(which causes onchocerciasis), *Wuchereria bancrofti* (which causes lymphatic filariasis), and 
*Necator americanus* (which causes soil-transmitted helminthiasis). 

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
.. |image7| image:: ../_static/P3_image7.png
.. |image8| image:: ../_static/P3_image8.png
            :width: 600
.. |image9| image:: ../_static/P3_image9.png
.. |image10| image:: ../_static/P3_image10.png
.. |image11| image:: ../_static/P3_image11.png

