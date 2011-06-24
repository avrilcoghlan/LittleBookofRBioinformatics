Comparative Genomics
====================

.. highlight:: r

Introduction
------------

*Comparative genomics* is the field of bioinformatics that involves
comparing the genomes of two different species, or of two different
strains of the same species.

One of the first questions to ask when comparing the genomes of two
species is: do the two species have the same number of genes (ie.
the same *gene content*)? Since all life on earth shared a common
ancestor at some point, any two species, for example, human and a
fruitfly, must have descended from a common ancestor species. 

Since the time of the common ancestor of two species (eg. of human and
mouse), some of the genes that were present in the common ancestor
species may have been lost from either of the two descendant
lineages. Furthermore, the two descendant lineages may have gained
genes that were not present in the common ancestor species.

Using the biomaRt R Library to Query the Ensembl Database
---------------------------------------------------------

To carry out comparative genomic analyses of two animal species whose
genomes have been fully sequenced (eg. human and mouse), it is useful to 
analyse the data in the Ensembl database (`www.ensembl.org <http://www.ensembl.org>`_). 

The main Ensembl database which you can browse on the 
`main Ensembl webpage <http://www.ensembl.org>`_ 
contains genes from fully sequenced vertebrates, as
well as *Saccharomyces cerevisiae* (yeast) and a small number of
additional model organism animals (eg. the nematode worm *Caenorhabditis elegans* 
and the fruit-fly *Drosophila melanogaster*).

There are also Ensembl databases for other groups of organisms,
for example `Ensembl Protists <http://protists.ensembl.org/index.html>`_ for
Protists, `Ensembl Metazoa <http://metazoa.ensembl.org/index.html>`_ for 
Metazoans, `Ensembl Bacteria <http://bacteria.ensembl.org/index.html>`_ for Bacteria,
`Ensembl Plants <http://plants.ensembl.org/index.html>`_ for Plants, and
`Ensembl Fungi <http://fungi.ensembl.org/index.html>`_ for Fungi.

It is possible to carry out analyses of the Ensembl database using
R, with the "biomaRt" R package. The "biomaRt" package can connect
to the Ensembl database, and perform queries on the data. 

The "biomaRt" R package is part of the Bioconductor set of R packages,
and so can be installed as explained `here <./installr.html#how-to-install-a-bioconductor-r-package>`_.

Once you have installed the "biomaRt" package, you can get a list
of databases that can be queried using this package by typing:

::

    > library("biomaRt") # Load the biomaRt package in R
    > listMarts()        # List all databases that can be queried
                                       biomart
      1                                ensembl
      2                                    snp
      3                    functional_genomics
      4                                   vega
      5                       bacterial_mart_9
      6                          fungal_mart_9
      7                    fungal_variations_9
      8                         metazoa_mart_9
      9                   metazoa_variations_9
      10                          plant_mart_9
      11                    plant_variations_9
      12                        protist_mart_9
      13                  protist_variations_9
      14                                   msd
      15                                  htgt
      16                              REACTOME
      17                         WS220-testing
      ...    
                                                          version
      1                              ENSEMBL GENES 62 (SANGER UK)
      2                         ENSEMBL  VARIATION 62 (SANGER UK)
      3                ENSEMBL FUNCTIONAL GENOMICS 62 (SANGER UK)
      4                                      VEGA 42  (SANGER UK)
      5                               ENSEMBL BACTERIA 9 (EBI UK)
      6                                  ENSEMBL FUNGI 9 (EBI UK)
      7                        ENSEMBL FUNGI VARIATION 9 (EBI UK)
      8                                ENSEMBL METAZOA 9 (EBI UK)
      9                      ENSEMBL METAZOA VARIATION 9 (EBI UK)
      10                                ENSEMBL PLANTS 9 (EBI UK)
      11                      ENSEMBL PLANTS VARIATION 9 (EBI UK)
      12                              ENSEMBL PROTISTS 9 (EBI UK)
      13                    ENSEMBL PROTISTS VARIATION 9 (EBI UK)
      14                                             MSD (EBI UK)
      15                  WTSI MOUSE GENETICS PROJECT (SANGER UK)
      16                                       REACTOME (CSHL US)
      17                                   WORMBASE 220 (CSHL US)
      ...

The names of the databases are listed, and then an explanation of what each
database is, and what is the version of the database.

You will see that the "biomaRt" R package can actually be used to
query many different databases including WormBase, UniProt,
Ensembl, etc. 

Here, we will discuss using the
"biomaRt" package to query the Ensembl database, but it is worth
remembering that it also be used to perform queries on other
databases such as UniProt. 

You can see above that "biomaRt" tells
you which version of each database can be searched, for example,
the version of the main Ensembl database that can be searched is Ensembl 62 (the current
release), while the version of the Ensembl Protists database that can be searched is
Ensembl Protists 9.

If you want to perform a query on the Ensembl database using
"biomaRt", you first need to specify that this is the database that
you want to query. You can do this using the useMart() function
from the "biomaRt" package:

::

    > ensemblprotists <- useMart("protist_mart_9") # Specify that we want to query the Ensembl Protists database

This tells "biomaRt" that you want to query the Ensembl Protists database.
The Ensembl Protists database contains data sets of genomic information for
different protist species whose genomes have been fully sequenced.

To see which data sets you can
query in the database that you have selected (using useMart()), you
can type:

::

    > listDatasets(ensemblprotists)         # List the data sets in the Ensembl Protists database
                dataset                                   description
      1      pramorum_eg_gene         Phytophthora ramorum genes (Phyra1_1)
      2        pvivax_eg_gene                Plasmodium vivax genes (EPr 2)
      3   pfalciparum_eg_gene           Plasmodium falciparum genes (2.1.4)
      4  ptricornutum_eg_gene      Phaeodactylum tricornutum genes (Phatr2)
      5     pchabaudi_eg_gene          Plasmodium chabaudi genes (May_2010)
      6   ddiscoideum_eg_gene Dictyostelium discoideum genes (dictybase.01)
      7        lmajor_eg_gene    Leishmania major strain Friedlin genes (1)
      ...
        version
      1      Phyra1_1
      2         EPr 2
      3         2.1.4
      4        Phatr2
      5      May_2010
      6  dictybase.01
      7             1
      ...

You will see a long list of the organisms for which the Ensembl Protists
database has genome data, including *Plasmodium vivax* and *Plasmodium falciparium* (which cause malaria),
and *Leishmania major*, which causes `leishmaniasis <http://www.who.int/leishmaniasis/en/>`_, which is
classified by the WHO as a neglected tropical disease.

To perform a query on the Ensembl database using the "biomaRt" R
package, you first need to specify which Ensembl data set your
query relates to. You can do this using the useDataset() function
from the "biomaRt" package. For example, to specify that you want
to perform a query on the Ensembl Leishmania major data set, you would type:

::

    > ensemblleishmania <- useDataset("lmajor_eg_gene",mart=ensemblprotists) 

Note that the name of the *Leishmania major* Ensembl data set is
"lmajor\_eg\_gene"; this is the data set listed for *Leishmania major*
genomic information when we typed listDatasets(ensemblprotists) above.

Once you have specified the particular Ensembl data set that you
want to perform a query on, you can perform the query using the
getBM() function from the "biomaRt" package. 

Usually, you will want to perform a query to a particular set of features from the *Leishmania major*
Ensembl data set. What types of features can you search for? You
can find this out by using the listAttributes() function from the
"biomaRt" package:

::

    > leishmaniaattributes <- listAttributes(ensemblleishmania)

The listAttributes() function returns a list object, the first
element of which is a vector of all possible features that you can
select, and the second element of which is a vector containing
explanations of all those features:

::

    > attributenames <- leishmaniaattributes[[1]]
    > attributedescriptions <- leishmaniaattributes[[2]]
    > length(attributenames)                     # Find the length of vector "attributenames"
     [1] 292
    > attributenames[1:10]                       # Print out the first 10 entries in vector "attributenames"
     [1] "ensembl_gene_id"                "ensembl_transcript_id"         
     [3] "ensembl_peptide_id"             "canonical_transcript_stable_id"
     [5] "description"                    "chromosome_name"               
     [7] "start_position"                 "end_position"                  
     [9] "strand"                         "band"       
    > attributedescriptions[1:10]                # Print out the first 10 entries in vector "attributedescriptions"
    > attributedescriptions[1:10]   
     [1] "Ensembl Gene ID"                   "Ensembl Transcript ID"            
     [3] "Ensembl Protein ID"                "Canonical transcript stable ID(s)"
     [5] "Description"                       "Chromosome Name"                  
     [7] "Gene Start (bp)"                   "Gene End (bp)"                    
     [9] "Strand"                            "Band"     

This gives us a very long list of 292 features in the *Leishmania major* Ensembl
data set that we can search for by querying the database, such as
genes, transcripts (mRNAs), peptides (proteins), chromosomes, GO (Gene Ontology) terms, and so on.

When you are performing a query on the Ensembl *Leishmania major* data set using
getBM(), you have to specify which of these features you want to
retrieve. For example, you can see from the output of
listAttributes() (see above) that one possible type of feature we
can search for are *Leishmania major* genes. To retrieve a list of all *Leishmania major*
genes from the *Leishmania major* Ensembl data set, we just need to type:

::

    > leishmaniagenes <- getBM(attributes = c("ensembl_gene_id"), mart=ensemblleishmania)

This returns a list variable *leishmaniagenes*, the first element of which
is a vector containing the names of all *Leishmania major* genes. Thus, to find
the number of genes, and print out the names of the first ten genes
stored in the vector, we can type:

::

    > leishmaniagenenames <- leishmaniagenes[[1]] # Get the vector of the names of all L. major genes
    > length(leishmaniagenenames) 
    [1] 9379 
    > leishmaniagenenames[1:10]
    [1] "LmjF.01.0010" "LmjF.01.0020" "LmjF.01.0030" "LmjF.01.0040" "LmjF.01.0050"
    [6] "LmjF.01.0060" "LmjF.01.0070" "LmjF.01.0080" "LmjF.01.0090" "LmjF.01.0100"

This tells us that there are 9379 different *Leishmania major* genes in the
*L. major* Ensembl data set. Note that this includes various types of
genes including protein-coding genes (both "known" and "novel"
genes, where the "novel" genes are gene predictions that don't have
sequence similarity to any sequences in sequence databases), RNA
genes, and pseudogenes.

What if we are only interested in protein-coding genes? If you look at the output
of listAttributes(ensemblleishmania), you will see that one of the features
is "gene\_biotype", which is tells us what sort of gene each gene
is (eg. protein-coding, pseudogene, etc.):

::

    > leishmaniagenes2 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=ensemblleishmania)

In this case, the getBM() function will return a list variable
*leishmaniagenes2*, the first element of which is a vector containing the
names of all *Leishmania major* genes, and the second of which is a vector
containing the types of those genes:

::

    > leishmaniagenenames2 <- leishmaniagenes2[[1]] # Get the vector of the names of all L. major genes
    > leishmaniagenebiotypes2 <- leishmaniagenes2[[2]] # Get the vector of the biotypes of all genes

We can make a table of all the different types of genes using the
table() function:

::

    > table(leishmaniagenebiotypes2) 
      leishmaniagenebiotypes2
             ncRNA nontranslating_cds     protein_coding         pseudogene 
                84                  2               8310                 90 
              rRNA             snoRNA              snRNA               tRNA 
                63                741                  6                 83 

This tells us that there are 8310 protein-coding genes, 90
pseudogenes, and various types of RNA genes (tRNA genes, rRNA
genes, snRNA genes, etc.). Thus, there are 8310 human
protein-coding genes.

Comparing the number of genes in two species
--------------------------------------------

Ensembl is a very useful resource for comparing the gene content of
different species. For example, one simple question that we can ask
by analysing the Ensembl data is: how many protein-coding genes are
there in *Leishmania major*, and how many in *Plasmodium falciparum*? 

We know how many protein-coding genes are in *Leishmania major* (8310; see above), but what
about *Plasmodium falciparum*? To answer this question, we first need to tell the
"biomaRt" package that we want to make a query on the Ensembl *Plasmodium falciparum*
data set. 

We can do this using the useDataset() function to select
the *Plasmodium falciparum* Ensembl data set. 

::

    > ensemblpfalciparum <- useDataset("pfalciparum_eg_gene",mart=ensemblprotists) 

Note that the name of the *Plasmodium falciparum* Ensembl data set is
"pfalciparum\_eg\_gene"; this is the data set listed for *Plasmodium falciparum*
genomic information when we typed listDatasets(ensemblprotists) above.

We can then use getBM() as above to retrieve the names of all *Plasmodium falciparum*
protein-coding genes. This time we have to set the "mart" option in
the getBM() function to "ensemblpfalciparum", to specify that we want to
query the *Plasmodium falciparum* Ensembl data set rather than the *Leishmania major* Ensembl data
set:

::

    > pfalciparumgenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=ensemblpfalciparum)
    > pfalciparumgenenames <- pfalciparumgenes[[1]] # Get the names of the P. falciparum genes
    > length(pfalciparumgenenames)                  # Get the number of P. falciparum genes
    [1] 6213  
    > pfalciparumgenebiotypes <- pfalciparumgenes[[2]] # Get the types of the P. falciparum genes 
    > table(pfalciparumgenebiotypes)
      pfalciparumgenebiotypes
              ncRNA non_translating_cds      protein_coding                rRNA 
                712                   1                5428                  24 
              snRNA                tRNA 
                  3                  45 


This tells us that there are 5428 *Plasmodium falciparum* protein-coding genes in
Ensembl. That is, *Plasmodium falciparum* seems to have less protein-coding
genes than *Leishmania major* (8310 protein-coding genes; see above).

It is interesting to ask: why does *Plasmodium falciparum* have less protein-coding
genes than *Leishmania major*? There are several possible explanations: (i) that
there have been gene duplications in the *Leishmania major* lineage since *Leishmania*
and *Plasmodium* shared a common ancestor, which gave rise to new *Leishmania major*
genes; (ii) that completely new genes (that are not related to any
other *Leishmania major* gene) have arisen in the *Leishmania major* lineage since *Leishmania* and
*Plasmodium* shared a common ancestor; or (iii) that there have been genes
lost from the *Plasmodium falciparum* genome since *Leishmania* and *Plasmodium* shared a common
ancestor.

To investigate which of these explanations is most likely to be
correct, we need to figure out how the *Leishmania major* protein-coding genes
are related to *Plasmodium falciparum* protein-coding genes.

Identifying homologous genes between two species
------------------------------------------------

xxx
The Ensembl database groups homologous (related) genes together
into gene families. If a gene from human and a gene from mouse are
related, they should be placed together into the same Ensembl gene
family. In fact, if a human gene has any homologues (related
genes), it should be placed into some Ensembl gene family.

For all human and mouse genes that are placed together in a gene
family, Ensembl classifies the relationship between each pair of
human and mouse genes as *orthologues* (related genes that shared a
common ancestor in the ancestor of human and mouse, and arose due
to the human-mouse speciation event) or *paralogues* (related genes
that arose due to a duplication event within a species, for
example, due to a duplication event in mouse, or a duplication
event in the human-mouse ancestor).

If you type listAttributes(ensembl) again, you will see that one
possible feature that you can search for is "mouse\_ensembl\_gene",
which is the mouse orthologue of a human Ensembl gene. Another
possible feature that you can search for is
"mouse\_orthology\_type", which describes the type of orthology
relationship between a particular human gene and its mouse
orthologue. For example, if a particular human gene has two mouse
orthologues, the relationship between the human gene and each of
the mouse orthologues will be "ortholog\_one2many"
(one-human-to-many-mouse orthology). This can arise in the case
where there was a duplication in the mouse lineage after human and
mouse diverged, which means that two different mouse genes (which
are paralogues of each other) are both orthologues of the same
human gene.

Therefore, we can retrive the Ensembl identifiers of the mouse
orthologues of all human protein-coding genes by typing:

::

    > humgenes4 <- getBM(attributes = c("ensembl_gene_id", "mouse_ensembl_gene", "mouse_orthology_type"), mart=ensembl)

This will return an R list variable *humgenes4*, the first element
of which is a vector of Ensembl identifiers for all human
protein-coding genes, and the second element of which is a vector
of Ensembl identifiers for their mouse orthologues, and the third
element of which is a vector with information on the orthology
types.

We can print out the names of the first 10 human genes and their
mouse orthologues, and their orthology types, by typing:

::

    > humgenenames4 <- humgenes4[[1]]            # Get the names of all human genes
    > hummouseorthologues4 <- humgenes4[[2]]     # Get the names of the mouse orthologues of all human genes
    > hummouseorthologuetypes4 <- humgenes4[[3]] # Get the orthology relationship type
    > humgenenames4[1:10] 
     [1] "ENSG00000211890" "ENSG00000211892" "ENSG00000211892" "ENSG00000211892" "ENSG00000211892"
     [6] "ENSG00000211891" "ENSG00000211891" "ENSG00000211895" "ENSG00000211893" "ENSG00000211893"
    > hummouseorthologues4[1:10] 
     [1] "ENSMUSG00000076610" "ENSMUSG00000076614" "ENSMUSG00000076612" "ENSMUSG00000076613"
     [5] "ENSMUSG00000076615" "ENSMUSG00000076611" "ENSMUSG00000087642" "ENSMUSG00000076610"
     [9] "ENSMUSG00000076614" "ENSMUSG00000076612"
    > hummouseorthologuetypes4[1:10] 
     [1] "ortholog_one2many"  "ortholog_many2many" "ortholog_many2many" "ortholog_many2many"
     [5] "ortholog_many2many" "ortholog_one2many"  "ortholog_one2many"  "ortholog_one2many" 
     [9] "ortholog_many2many" "ortholog_many2many"

Not all human genes have mouse orthologues. To find out how many
human genes are orthologues, we can first find the indices of the
elements of the vector *hummouseorthologues4* that are empty:

::

    > myindex4 <- hummouseorthologues4=="" 

We can then find out the names of the human genes corresponding to
those indices:

::

    > humgenenames4b <- humgenenames4[myindex4]
    > length(unique(humgenenames4b))
    [1] 34115

This tells us that 34,115 human genes do not have mouse
orthologues. Note that we have to use the unique() function (which
removes duplicates from a vector) to count the number of human gene
names in vector *humgenenames4b*, as some human gene names appear
twice in that vector (because they have more than one mouse
orthologue listed in vector *hummouseorthologues4*).

How many of the 34,115 human genes that do not have mouse
orthologues are protein-coding genes? To answer this question, we
can merge together the information in the R list variable
*humgenes2* (which contains information on the name of each human
gene and its type), and the R list variable *humgenes4*. This can
be done using the merge() function in R, which can merge together
two list variables that contain some named elements in common (in
this case, both list variables contain a vector that has the names
of human genes):

::

    > humgenes5 <- merge(humgenes2, humgenes4)

The first element of the merged list variable *humgenes5* contains
a vector of the human gene names, the second has a vector of the
types of those genes (eg. protein-coding, pseudogene etc.), and the
third element has a vector of the mouse orthologues' names. We can
therefore find out how many protein-coding human genes lack mouse
orthologues by typing:

::

    > humgenenames5 <- humgenes5[[1]] 
    > humgenebiotypes5 <- humgenes5[[2]]
    > hummouseorthologues5 <- humgenes5[[3]]
    > myindex5 <- hummouseorthologues5=="" & humgenebiotypes5=="protein_coding"
    > humgenenames5b <- humgenenames5[myindex5] 
    > length(unique(humgenenames5b)) 
    [1] 4857 

This tells us that there are 4857 human protein-coding genes that
lack mouse orthologues.

Summary
-------

In this practical, you will have learnt to use the following R
functions:


#. useMart() to select a database to query (in the biomaRt package)
#. useDataset() to select a data set in a database to query (in the
   biomaRt package)
#. listDatasets() to get a list of all data sets in a database (in
   the biomaRt package)
#. listAttributes() to get a list of all features of a data set (in
   the biomaRt package)
#. getBM() to make a query on a database (in the biomaRt package)
#. unique() to remove duplicate elements from a vector
#. merge() to merge R list objects that contain some named elements
   in common

Links and Further Reading
-------------------------

Some links are included here for further reading, which will be
especially useful if you need to use the R package for your project
or assignments.

For background reading on comparative genomics, it is recommended
to read Chapter 8 of
*Introduction to Computational Genomics: a case studies approach*
by Cristianini and Hahn (Cambridge University Press;
`www.computational-genomics.net/book/ <http://www.computational-genomics.net/book/>`_).

Exercises
---------

Answer the following questions, using the R package. For each
question, please record your answer, and what you typed into R to
get this answer.

Q1. How many cow genes are there in the current version of the Ensembl database? 
    How many of the cow Ensembl genes are protein-coding genes?
Q2. How many cow protein-coding genes have human orthologues? 
    How many of the cow protein-coding genes have one-to-one
    orthologues in human?
Q3. How many cow genes have Pfam domains? Q4. What are the top 5 most common Pfam domains in cow genes, and how many copies of each are there in the cow protein set? Q5. How many of copies are there in the human protein set, of each of the top 5 cow protein domains? 
    Are the numbers of copies of some domains different in the two
    species?
    How would you check if this is a statistically significant
    difference?


