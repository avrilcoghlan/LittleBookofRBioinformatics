REVISION EXERCISES 2
====================

These are some revision exercises on sequence alignment and phylogenetic trees.

Exercises
---------

Answer the following questions. For each question, please record
your answer, and what you did/typed to get this answer.

Model answers to the exercises are given in 
`Answers to Revision Exercises 2 <./revisionexercises_answers.html#revision-exercises-2>`_.

Q1.
---
One of the key proteins produced by rabies virus is the rabies phosphoprotein (also known as rabies virus protein P).  The UniProt accession for rabies virus phosphoprotein is P06747. The Mokola virus also produces a phosphoprotein, which has UniProt accession P0C569. Use the dotPlot() function in the SeqinR R library to make a dotplot of these two proteins, using a windowsize of 10 and threshold of 5. Are there any long regions of similarity between the two proteins (if so, where are they)?
    Note: rabies virus is the virus responsible for `rabies <http://www.who.int/rabies/en/>`_, which is classified by the WHO as a neglected tropical disease. Mokola virus and rabies virus are closely related viruses that both belong to a group of viruses called the Lyssaviruses. Mokola virus causes a rabies-like infection in mammals including humans.

Q2.
---
The function "makeDotPlot1()" below is an R function that makes a dotplot of two sequences by plotting a dot at every position where the two sequences share an identical letter.  Use this function to make a dotplot of the rabies virus phosphoprotein and the Mokola virus phosphoprotein, setting the argument "dotsize" to 0.1 (this determines the radius of each dot plotted). Are there any long regions of similarity between the two proteins (if so, where are they)? Do you find the same regions as found in Q1, and if not, can you explain why?

.. highlight:: r

::

    > makeDotPlot1 <- function(seq1,seq2,dotsize=1)
      {
         length1 <- length(seq1)
         length2 <- length(seq2)
         # make a plot:
         x <- 1
         y <- 1 
         plot(x,y,ylim=c(1,length2),xlim=c(1,length1),col="white") # make an empty plot
         # now plot dots at every position where the two sequences have the same letter:
         for (i in 1:length1)
         {
            letter1 <- seq1[i]
            for (j in 1:length2)
            {
               letter2 <- seq2[j]
               if (letter1 == letter2)
               {
                  # add a point to the plot
                  points(x=i,y=j,cex=dotsize,col="blue",pch=7)
               }   
            }
         }
      }

Q3.
---
Adapt the R code in Q2 to write a function that makes a dotplot using a window of size *x* letters, where a dot is plotted in the first  cell of the window if *y* or more letters compared in that window are identical in the two sequences.  

Q4.
---
Use the dotPlot() function in the SeqinR R library to make a dotplot of rabies virus phosphoprotein and Mokola virus phosphoprotein, using a window size of 3 and a threshold of 3. Use your own R function from Q3 to make a dotplot of rabies virus phosphoprotein and Mokola virus phosphoprotein, using a windowsize (*x*) of 3 and a threshold (*y*) of 3. Are the two plots similar or different, and can you explain why?

Contact
-------

I will be grateful if you will send me (`Avril Coghlan <http://www.ucc.ie/microbio/avrilcoghlan/>`_) corrections or suggestions for improvements to
my email address a.coghlan@ucc.ie 

License
-------

The content in this book is licensed under a `Creative Commons Attribution 3.0 License
<http://creativecommons.org/licenses/by/3.0/>`_.


