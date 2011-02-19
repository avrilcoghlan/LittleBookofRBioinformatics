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

Q4. *How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the RefSeq part of the 
    NCBI Sequence Database?*

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

