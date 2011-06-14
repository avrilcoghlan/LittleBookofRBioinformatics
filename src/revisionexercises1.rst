Revision Exercises 1
====================

These are some revision exercises on sequence statistics and sequence
databases.

Exercises
---------

Answer the following questions. For each question, please record
your answer, and what you did/typed to get this answer.

Model answers to the exercises are given in the chapter entitled
`Answers to Revision Exercises 1 <./revisionexercises1_answers.html>`_.

Q1. What is the length of (total number of base-pairs in) the *Schistosoma mansoni* mitochondrial genome
(NCBI accession NC\_002545), and how many As, Cs, Gs and Ts does it contain?
    Note: *Schistmosoma mansoni* is a parasitic worm that is responsible for causing 
    `schistosomiasis <http://apps.who.int/tdr/svc/diseases/schistosomiasis>`_, 
    which is classified by the WHO as a neglected tropical disease.

Q2. What is the length of the *Brugia malayi* mitochondrial genome (NCBI accession NC\_004298),
and how many As, Cs, Gs and Ts does it contain?
    Note: *Brugia malayi* is a parasitic worm responsible for causing
    `lymphatic filariasis <http://apps.who.int/tdr/svc/diseases/lymphatic-filariasis>`_,
    which is classified by the WHO as a neglected tropical disease.

Q3. What is the probability of the *Brugia malayi* mitochondrial genome sequence (NCBI accession NC\_004298), 
   according to a multinomial model in which the probabilities of As, Cs, Gs and Ts (pA, pC, pG, and pT) 
    are set equal to the fraction of As, Cs, Gs and Ts in the *Schistosoma mansoni* mitochondrial genome?

Q4. What are the top three most frequent 4-bp words (4-mers) in the genome of the
    bacterium *Chlamydia trachomatis* strain D/UW-3/CX (NCBI accession NC\_000117), and
    how many times do they occur in its sequence?
    Note: *Chlamydia trachomatis* is a bacterium responsible for 
    `trachoma <http://www.who.int/blindness/causes/priority/en/index2.html>`_, which is
    classified by the WHO as a neglected tropical disease. 

Q5. Write an R function to generate a random DNA sequence that is *n* letters long (that is, 
    *n* bases long) using a multinomial model in which the probabilities *pA*, *pC*, *pG*, 
    and *pT* are set equal to the fraction of As, Cs, Gs and Ts in the *Schistosoma mansoni*
    mitochondrial genome (here *pA* stands for the probability of As, *pC* is the probability of Cs, etc.)
    Hint: look at the help page for the "sample()" function in R, as it might be useful to use within your R function.

Q6. Give an example of using your function from Q5 to calculate a random sequence that is 20 letters 
    long, using a multinomial model with *pA* =0.28, *pC* =0.21, *pG* =0.22, and *pT* =0.29 .

Q7. How many protein sequences from rabies virus are there in the NCBI Protein database?
    Note: rabies virus is the virus responsible for 
    `rabies <http://www.who.int/rabies/en/>`_, which is classified by the WHO as a neglected
    tropical disease.

Q8. What is the NCBI accession for the Mokola virus genome?
    Note: Mokola virus and rabies virus are closely related viruses that both belong to a group of 
    viruses called the Lyssaviruses. Mokola virus causes a rabies-like infection in mammals including humans.

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.


