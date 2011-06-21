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
fruitfly, must have descended from a common ancestor species. Since
the time of the common ancestor of two species (eg. of human and
mouse), some of the genes that were present in the common ancestor
species may have been lost from either of the two descendant
lineages. Furthermore, the two descendant lineages may have gained
genes that were not present in the common ancestor species.

Using the biomaRt R Library to Query the Ensembl Database
---------------------------------------------------------

To carry out comparative genomic analyses of two vertebrate species
(eg. human and mouse), it is useful to analyse the data in the
Ensembl database (`www.ensembl.org <http://www.ensembl.org>`_). The
Ensembl database contains genes from fully sequenced vertebrate, as
well as *Saccharomyces cerevisiae* (yeast) and a small number of
additional model organism animals (eg. the nematode worm
*Caenorhabditis elegans* and the fruit-fly
*Drosophila melanogaster*).

It is possible to carry out analyses of the Ensembl database using
R, with the "biomaRt" R library. The "biomaRt" library can connect
to the Ensembl database, and perform queries on the data. The
"biomaRt" R library is part of the Bioconductor set of R libraries,
and so can be installed by typing:

::

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("biomaRt")

Once you have installed the "biomaRt" library, you can get a list
of databases that can be queried using this library by typing:

::

    > library("biomaRt") # Load the biomaRt library in R
    > listMarts()        # List all databases that can be queried
                         biomart                                                 version
    1                    ensembl                            ENSEMBL GENES 57 (SANGER UK)
    2                        snp                        ENSEMBL VARIATION 57 (SANGER UK)
    3        functional_genomics              ENSEMBL FUNCTIONAL GENOMICS 57 (SANGER UK)
    4                       vega                                     VEGA 37 (SANGER UK)
    5           bacterial_mart_4                             ENSEMBL BACTERIA 4 (EBI UK)
    6              fungal_mart_4                               ENSEMBL FUNGAL 4 (EBI UK)
    7             metazoa_mart_4                              ENSEMBL METAZOA 4 (EBI UK)
    8               plant_mart_4                                ENSEMBL PLANT 4 (EBI UK)
    9             protist_mart_4                             ENSEMBL PROTISTS 4 (EBI UK)
    10                       msd                                  MSD PROTOTYPE (EBI UK)
    11                      htgt HIGH THROUGHPUT GENE TARGETING AND TRAPPING (SANGER UK)
    12                  REACTOME                                      REACTOME (CSHL US)
    13          wormbase_current                                      WORMBASE (CSHL US)
    ...

You will see that the "biomaRt" R library can actually be used to
query many different databases including WormBase, UniProt,
Ensembl, etc. In today's practical, we will discuss using the
"biomaRt" library to query the Ensembl database, but it is worth
remembering that it also be used to perform queries on other
databases such as UniProt. You can see above that "biomaRt" tells
you which version of each database can be searched, for example,
the Ensembl version that can be searched is Ensembl 57 (the current
release).

If you want to perform a query on the Ensembl database using
"biomaRt", you first need to specify that this is the database that
you want to query. You can do this using the useMart() function
from the "biomaRt" library:

::

    > ensembl <- useMart("ensembl") # Specify that we want to query the Ensembl database

This tells "biomaRt" that you want to query the Ensembl database.
The Ensembl database contains data sets of genomic information for
many different vertebrate species. To see which data sets you can
query in the database that you have selected (using useMart()), you
can type:

::

    > listDatasets(ensembl)         # List the data sets in the Ensembl database
    1          oanatinus_gene_ensembl        Ornithorhynchus anatinus genes (OANA5)
    2           tguttata_gene_ensembl       Taeniopygia guttata genes (taeGut3.2.4)
    3         cporcellus_gene_ensembl               Cavia porcellus genes (cavPor3)
    4         gaculeatus_gene_ensembl        Gasterosteus aculeatus genes (BROADS1)
    5          lafricana_gene_ensembl            Loxodonta africana genes (loxAfr2)
    6         mlucifugus_gene_ensembl              Myotis lucifugus genes (myoLuc1)
    7           hsapiens_gene_ensembl                   Homo sapiens genes (GRCh37)
    ...

You will see a long list of the organisms for which the Ensembl
database has genome data, including *Ornithorhynchus anatinus*
(platypus), *Taeniopygia guttata* (zebra finch), *Cavia porcellus*
(guinea pig), *Gasterosteus aculeatus* (three-spined stickleback),
*Loxodonta africana* (African elephant), *Myotis lucifugus* (little
brown bat), *Homo sapiens* (human), and so on..

To perform a query on the Ensembl database using the "biomaRt" R
library, you first need to specify which Ensembl data set your
query relates to. You can do this using the useDataset() function
from the "biomaRt" library. For example, to specify that you want
to perform a query on the Ensembl human data set, you would type:

::

    > ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl) 

Note that the name of the human Ensembl data set is
"hsapiens\_gene\_ensembl"; this is the data set listed for human
genomic information when we typed listDatasets(ensembl) above.

Once you have specified the particular Ensembl data set that you
want to perform a query on, you can perform the query using the
getBM() function from the "biomaRt" library. Usually, you will want
to perform a query to a particular set of features from the human
Ensembl data set. What types of features can you search for? You
can find this out by using the listAttributes() function from the
"biomaRt" library:

::

    > attributes <- listAttributes(ensembl)

The listAttributes() function returns a list object, the first
element of which is a vector of all possible features that you can
select, and the second element of which is a vector containing
explanations of all those features:

::

    > attributenames <- attributes[[1]]
    > attributedescriptions <- attributes[[2]]
    > length(attributenames)                     # Find the length of vector "attributenames"
    [1] 961
    > attributenames[1:10]                       # Print out the first 10 entries in vector "attributenames"
     [1] "ensembl_gene_id"                "ensembl_transcript_id"          "ensembl_peptide_id"            
     [4] "canonical_transcript_stable_id" "description"                    "chromosome_name"               
     [7] "start_position"                 "end_position"                   "strand"                        
    [10] "band" 
    > attributedescriptions[1:10]                # Print out the first 10 entries in vector "attributedescriptions"
    > attributedescriptions[1:10]   
     [1] "Ensembl Gene ID"                   "Ensembl Transcript ID"            
     [3] "Ensembl Protein ID"                "Canonical transcript stable ID(s)"
     [5] "Description"                       "Chromosome Name"                  
     [7] "Gene Start (bp)"                   "Gene End (bp)"                    
     [9] "Strand"                            "Band"     

This gives us a very long list of 961 features in the human Ensembl
data set that we can search for by querying the database, such as
human genes, human transcripts (mRNAs), human peptides (proteins),
chromosomes, GO (Gene Ontology) terms, and so on.

When you are performing a query on the Ensembl human data set using
getBM(), you have to specify which of these features you want to
retrieve. For example, you can see from the output of
listAttributes() (see above) that one possible type of feature we
can search for are human genes. To retrieve a list of all human
genes from the human Ensembl data set, we just need to type:

::

    > humgenes <- getBM(attributes = c("ensembl_gene_id"), mart=ensembl)

This returns a list variable *humgenes*, the first element of which
is a vector containing the names of all human genes. Thus, to find
the number of genes, and print out the names of the first ten genes
stored in the vector, we can type:

::

    > humgenenames <- humgenes[[1]] # Get the vector of the names of all human genes
    > length(humgenenames) 
    [1] 51682
    > humgenenames[1:10]
     [1] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
     [5] "ENSG00000000460" "ENSG00000000938" "ENSG00000000971" "ENSG00000001036"
     [9] "ENSG00000001084" "ENSG00000001167"

This tells us that there are 51,682 different human genes in the
human Ensembl data set. Note that this includes various types of
genes including protein-coding genes (both "known" and "novel"
genes, where the "novel" genes are gene predictions that don't have
sequence similarity to any sequences in sequence databases), RNA
genes, and pseudogenes.

As mentioned above, the 51,682 different human genes in the human
Ensembl data set probably include various classes of genes, such as
protein-coding genes, RNA genes, or pseudogenes. What if we are
only interested in protein-coding genes? If you look at the output
of listAttributes(ensembl), you will see that one of the features
is "gene\_biotype", which is tells us what sort of gene each gene
is (eg. protein-coding, pseudogene, etc.):

::

    > humgenes2 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=ensembl)

In this case, the getBM() function will return a list variable
*humgenes2*, the first element of which is a vector containing the
names of all human genes, and the second of which is a vector
containing the types of those genes:

::

    > humgenenames2 <- humgenes2[[1]] # Get the vector of the names of all human genes
    > humgenebiotypes2 <- humgenes2[[2]] # Get the vector of the biotypes of all genes

We can make a table of all the different types of genes using the
table() function:

::

    > table(humgenebiotypes2) 
               IG_C_gene            IG_D_gene            IG_J_gene            IG_V_gene 
                      21                   30                   93                  226 
                 lincRNA                miRNA     miRNA_pseudogene             misc_RNA 
                    3517                 1698                   18                 1564 
     misc_RNA_pseudogene              Mt_rRNA              Mt_tRNA   Mt_tRNA_pseudogene 
                       7                    2                   22                  580 
    processed_transcript       protein_coding           pseudogene                 rRNA 
                    6762                22320                 9456                  461 
         rRNA_pseudogene     scRNA_pseudogene               snoRNA    snoRNA_pseudogene 
                     338                  834                 1217                  457 
                   snRNA     snRNA_pseudogene      tRNA_pseudogene 
                    1441                  490                  128 

This tells us that there are 22,320 protein-coding genes, 9456
pseudogenes, and various types of RNA genes (tRNA genes, rRNA
genes, snRNA genes, etc.). Thus, there are 22,320 human
protein-coding genes.

Comparing the number of genes in two vertebrate species
-------------------------------------------------------

Ensembl is a very useful resource for comparing the gene content of
different species. For example, one simple question that we can ask
by analysing the Ensembl data is: how many protein-coding genes are
there in mouse, and how many in human? We know how many
protein-coding genes are in humans (22,320; see above), but what
about mouse? To answer this question, we first need to tell the
"biomaRt" library that we want to make a query on the Ensembl mouse
data set. We can do this using the useDataset() function to select
the mouse (*Mus musculus*) Ensembl data set:

::

    > ensembl2 <- useDataset("mmusculus_gene_ensembl",mart=ensembl) 

We can then use getBM() as above to retrieve the names of all mouse
protein-coding genes. This time we have to set the "mart" option in
the getBM() function to "ensembl2", to specify that we want to
query the mouse Ensembl data set rather than the human Ensembl data
set:

::

    > mousegenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=ensembl2)
    > mousegenenames <- mousegenes[[1]]    # Get the names of the mouse genes
    > length(mousegenenames)               # Get the number of mouse genes
    [1] 34213
    > mousegenebiotypes <- mousegenes[[2]] # Get the types of the mouse genes 
    > table(mousegenebiotypes)
                 IG_C_gene              IG_D_gene              IG_J_gene              IG_V_gene 
                        20                     15                     87                    361 
                   lincRNA                  miRNA               misc_RNA                Mt_rRNA 
                       495                   1081                    148                      2 
                   Mt_tRNA polymorphic_pseudogene   processed_transcript         protein_coding 
                        22                      1                   2208                  23062 
                pseudogene                   rRNA                 snoRNA                  snRNA 
                      4677                    222                    949             

This tells us that there are 23,062 mouse protein-coding genes in
Ensembl. That is, mouse seems to have slightly more protein-coding
genes than human (23,062 mouse genes versus 22,320 human genes).

It is interesting to ask: why does mouse have more protein-coding
genes than human? There are several possible explanations: (i) that
there have been gene duplications in the mouse lineage since mouse
and human shared a common ancestor, which gave rise to new mouse
genes; (ii) that completely new genes (that are not related to any
other mouse gene) have arisen in the mouse lineage since mouse and
human shared a common ancestor; or (iii) that there have been genes
lost from the human genome since mouse and human shared a common
ancestor.

To investigate which of these explanations is most likely to be
correct, we need to figure out how the human protein-coding genes
are related to mouse protein-coding genes.

Identifying homologous genes between two species
------------------------------------------------

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


#. useMart() to select a database to query (in the biomaRt library)
#. useDataset() to select a data set in a database to query (in the
   biomaRt library)
#. listDatasets() to get a list of all data sets in a database (in
   the biomaRt library)
#. listAttributes() to get a list of all features of a data set (in
   the biomaRt library)
#. getBM() to make a query on a database (in the biomaRt library)
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


