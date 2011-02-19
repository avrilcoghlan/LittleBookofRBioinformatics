Answers to the exercises on Sequence Databases
==============================================   

Q1. *What information about the DEN-1 Dengue virus sequence (NCBI accession NC\_001477) can you obtain from its annotations in the NCBI Sequence Database?*

To do this, you need to go to the `www.ncbi.nlm.nih.gov <http://www.ncbi.nlm.nih.gov>`_ website 
and type the DEN-1 Dengue virus genome sequence accession (NC\_001477) in the search box, and press 'Search'. 

On the search results page, you should see '1' beside the word 'Nucleotide', meaning that there was one hit to a sequence record in the NCBI Nucleotide database, which contains DNA and RNA sequences. If you click on the word 'Nucleotide', it will bring you to the sequence record, which should be the NCBI sequence record for the DEN-1 Dengue virus' genome (ie. for accession NC\_001477). 
	
On the webpage for the NCBI sequence record for the DEN-1 Dengue virus' genome (webpage http://www.ncbi.nlm.nih.gov/nuccore/NC_001477) , you will see in the DEFINITION, ORGANISM and REFERENCE fields of its NCBI record: 

DEFINITION: Dengue virus type 1, complete genome.
ORGANISM: Dengue virus 1
REFERENCE: There are several papers (the first is):
AUTHORS: Puri,B., Nelson,W.M., Henchal,E.A., Hoke,C.H., Eckels,K.H., Dubois,D.R., Porter,K.R. and Hayes,C.G.
TITLE: Molecular analysis of dengue virus attenuation after serial passage in primary dog kidney cells
JOURNAL: J. Gen. Virol. 78 (PT 9), 2287-2291 (1997)

There are also some other references, for papers published about the Dengue virus genome sequence. 

Q2. *What were the nucleotide sequences published in Nature volume 460, page 352?*

To do this, you need to go to the NCBI website at www.ncbi.nlm.nih.gov and type in the search 
box on the top: "Nature"[JOUR] AND 460[VOL] AND 352[PAGE]

Here [JOUR] specifies the journal name, [VOL] the volume of the journal the paper is in, and [PAGE] the page number.

This should bring up a results page with "50890" beside the word "Nucleotide", and "1" beside the word
"Genome", and "25700" beside the word "Protein", indicating that there were 50890 hits to sequence records in the Nucleotide database, 
which contains DNA and RNA sequences, and 1 hit to the Genome database, which contains genome sequences, and 25700
hits to the Protein database, which contains protein sequences.

If you click on the word "Nucleotide", it will bring up a webpage with a list of links to the NCBI sequence 
records for those 50890 hits. The 50890 hits are all contigs from the schistosome worm *Schistosoma mansoni*.

Likewise, if you click on the word "Protein", it will bring up a webpage with a list of links to the NCBI
sequence records for the 25700 hits, and you will see that the hits are all predicted proteins for *Schistosoma
mansoni*.

If you click on the word "Genome", it will bring you to the NCBI record for the *Schistosoma mansoni* genome
sequence, which has NCBI accession NS\_00200. Note that the accession starts with "NS\_", which indicates that
it is a RefSeq accession. 

Therefore, in *Nature* volume 460, page 352, the *Schistosoma mansoni* genome sequence was published, along
with all the DNA sequence contigs that were sequenced for the genome project, and all the predicted proteins
for the gene predictions made in the genome sequence. You can view the original paper on the *Nature* website
at `http://www.nature.com/nature/journal/v460/n7253/abs/nature08160.html <http://www.nature.com/nature/journal/v460/n7253/abs/nature08160.html>`_.

Q3. *How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the NCBI Sequence Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Nucleotide" from the drop-down list at the top 
of the webpage, as you want to search for nucleotide (DNA or RNA) sequences.

Then in the search box, type "Chlamydia trachomatis"[ORGN] and press 'Search'.

Here [ORGN] specifies the organism you are interested in, that is, the species name in Latin.

The results page should give you a list of the hits to sequence records in the NCBI Nucleotide database. 
It will say "Found 35385 nucleotide sequences.   Nucleotide (35237)   GSS (148)". 
This means that 35,385 sequences were found, of which 35237 are DNA or RNA sequences, and 
148 are DNA sequences from the Genome Sequence Surveys (GSS), that is, from 
genome sequencing projects [as of 19-Feb-2011]. Note that there are new sequences 
being added to the database continuously, so if you check this again in a couple of months, you will 
probably find a higher number of sequences (eg. 36,000 sequences).

Note: if you just go to the www.ncbi.nlm.nih.gov database, and search for "Chlamydia trachomatis"[ORGN] 
(without choosing "Nucleotide" from the drop-down list), you will see 35237 hits to the Nucleotide 
database and 148 to the GSS (Genome Sequence Survey) database.

Note also that if you search for "Chlamydia trachomatis", without using [ORGN] to specify the organism, 
you will get 51046 hits to the Nucleotide database and 149 to the GSS database, but some of these might 
not be *Chlamydia trachomatis* sequences – they could just be sequences for which the NCBI sequence 
record contains the phrase "Chlamydia trachomatis" somewhere.

Q4. *How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the RefSeq part of the NCBI Sequence Database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov and select "Nucleotide" from the drop-down list 
at the top of the webpage, as you want to search for nucleotide sequences.

Then in the search box, type "Chlamydia trachomatis"[ORGN] AND srcdb_refseq[PROP] and press 'Search'.

Here [ORGN] specifies the organism, and [PROP] specifies a property of the sequences (in this case that 
they belong to the RefSeq subsection of the NCBI database).

At the top of the results page, it should say "Results: 1 to 20 of 29 sequences" [as of 19-Feb-2011]. 
As for Q3, if you try this again in a couple of months, the number will probably be higher, due to extra 
sequences added to the database. 

Note that the sequences in Q3 are all *Chlamydia trachomatis* DNA and RNA sequences in the NCBI database. 
The sequences in Q4 gives the *Chlamydia trachomatis* DNA and RNA sequences in the RefSeq part of the NCBI 
database, which is a subsection of the database for high-quality manually-curated data. 

The number of sequences in RefSeq is much fewer than the total number of *C. trachomatis* sequences, 
partly because low quality sequences are never added to RefSeq, but also because RefSeq curators have 
probably not had time to add all high-quality sequences to RefSeq (this is a time-consuming process, 
as the curators add additional information to the NCBI Sequence records in RefSeq, such as references to 
papers that discuss a particular sequence). 

Q5. *How many nucleotide sequences were submitted to NCBI by Matthew Berriman?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Nucleotide" from the drop-down list, 
as you want to search for nucleotide sequences.

Then in the search box, type "Berriman M"[AU] and press 'Search'.

Here [AU] specifies the name of the person who either submitted the sequence to the NCBI database, 
or wrote a paper describing the sequence. 

On the results page, it should say at the top: "Found 460052 nucleotide sequences.   Nucleotide (250328)   EST (121075)   GSS (88649)". This means that 460052 DNA/RNA sequences were either submitted to the NCBI database by someone called M. Berriman, or were described in a paper by someone called M. Berriman. Of these, 250328 were DNA/RNA sequences, 121075 were EST sequences (part of mRNAs), and 88649 were DNA sequences from genome sequencing projects (GSS or Genome Sequence Survey sequences).

Note that unfortunately the NCBI website does not allow us to search for "Berriman Matthew"[AU] so we cannot be sure 
that all of these sequences were submitted by Matthew Berriman. 

Q6. *How many nucleotide sequences from nematode worms are there in the NCBI Database?*

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

Q7. *How many nucleotide sequences for collagen genes from nematode worms are there in the NCBI Database?*

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

Q8. *How many mRNA sequences for collagen genes from nematode worms are there in the NCBI Database?*

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

Q9. *How many protein sequences for collagen proteins from nematode worms are there in the NCBI database?*

To answer this, you need to go to www.ncbi.nlm.nih.gov, and select "Protein" from the drop-down list, 
as you want to search for protein sequences.

Then type in the search box: Nematoda[ORGN] AND collagen and press 'Search'.

On the results page, you should see '1 to 20 of 1886'. This means that 1886 protein sequences from nematode
worms were found that include the word collagen in the NCBI sequence entries [as of 19-Feb-2011].

Q10. *What is the accession number for the Trypanosoma cruzi genome in NCBI?*

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

Q11. *How many fully sequenced nematode worm species are represented in the NCBI Genome database?*

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

