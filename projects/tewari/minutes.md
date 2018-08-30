project road map: https://github.com/miRTop/incubator/projects/2

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

