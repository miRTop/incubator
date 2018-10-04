project road map: https://github.com/miRTop/incubator/projects/2

## 10-04-2018

People: Lorena, Ioannis, Marc, Phillipe, Lauro, Gianvito

* We discussed the equimolar data compared to plasma data to look at difference between distribution of what is the most expressed and what is the percentage of the rest of the sequences.

Ideas that came up:

* maybe too much noise from miR like let7*
* Ioannis confirmed that the equimolar sample may no be purified. 
* Use cutoff tool from Phillipe group, discard low expressed reads
* equimolar from miRGe tool

### Ready to do points:

* Lorena has to send FASTA file of all sequences to Ioannis to get the features that will be used for the machine learning analysis once we can classify sequence into reproducible or not.
* Marc to contact tewari's paper authors: 
 * if we can know what are the extra sequences they added from another company
 * whether they have equimolar and plasma samples or all is gone
* Lorena to analyze other data set with spike-ints: working on this

### Next meeting 10-18-2018, 10am boston time.

* create a scheme from RNA material to results and spot all possibles places where artifacts can be created. It would be good for the paper and create new questions.

## 09-20-2018

People: Lorena, Arun, Ioannis, Marc, Phillipe

* We discussed the equimolar data after being analyzed by bcbio pipeline: https://github.com/miRTop/incubator/blob/master/projects/tewari_equimolar/commonality.md

### The discussion was focused on finding too many isomiRs from the tewari equimolar data

Ideas to explain this:

* non purified miRXplore sample
* illumina software generating 5' end that are not real
* trimming tool generating 3' end that are not real
* if isomiRs are errors from an unique sequence, and this happens after ligation, there should be a correlation of isomiRs and reference expression

### Ready to do points:

* Lorena has to send FASTA file of all sequences to Ioannis to get the features that will be used for the machine learning analysis once we can classify sequence into reproducible or not.
* Ioannis to contact Qiagen and miRXplore company
* Marc to contact tewari's paper authors
* Lorena to analyze other data set with spike-ints
* Lorena to share the fasta files to be analized with miRGe2 (maybe?)
* Lorena to compare distribution of plasma sample and equimolar to see if they are the same


## 08-30-2018

People: Lorena, Aarun, Ioannis

* We presented the variation in the output of IDR package when changing two of the main paramters. See [Methods](https://github.com/miRTop/incubator/tree/master/projects/tewari#methods)
* We presented an example of figure when visualizing the [commonality among replicates and labs](https://github.com/miRTop/incubator/blob/master/projects/tewari/figures/labs/bcbio.png)

### The discussion was focused on finding a method to flag sequences as reproducible or not:

* One strategy was to let Ioannis to contact IDR authors to ask for advice about how to adapt the code to our question
* Another strategy was to let Lorena find some statistician to join the collaboration to help in this particular step

### Ready to do points:

* Lorena has to send FASTA file of all sequences to Ioannis to get the features that will be used for the machine learning analysis once we can classify sequence into reproducible or not.
* Aarun will update plots and divided by low/medium/high expression


## 08-09-2018

We presented the data and results until now and discussed the following points:

* better analysis for variant calling: Ioannis
* chipseq-metrics-information to measure reproducibility: Ioannis will send information around
* mirnas matching the adapter or complemtary to it to measure the effect in expression due to the adapter used in the protocol
* measure the edit distance of the isomirs to all the mirnas in the database to spot problematic isomirs that we should remove
* Repeat plots separating between low/medium/high
* kmer analysis to explain difference among labs

### Papers related to reproducibility and small RNA bias by Ioannis

The metric is called IDR (irreproducible discovery rate) and you can find it here:

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/
 
http://cran.r-project.org/web/packages/idr/index.htm

These are the papers for the bias:
 
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0126049
 
https://www.ncbi.nlm.nih.gov/pubmed/22647250
 

### Questions for next meeting

* what to do with isomiRs that map equally to other mature miRNAs

### Final goal of the project

* How to filter False isomirs and True isomirs

Ideally:

* Machine learning to predict the real isomirs and remove the ones due to technical errors

features to look into that:

* PCR artifact
* Ligation
* Size selection
* Steps in protocols

