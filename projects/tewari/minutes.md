project road map: https://github.com/miRTop/incubator/projects/2

## 08-090-2018

We presented the data and results until now and discussed the following points:

* better analysis for variant calling: Ioannis
* chipseq-metrics-information to measure reproducibility: Ioannis will send information around
* mirnas matching the adapter or complemtary to it to measure the effect in expression due to the adapter used in the protocol
* measure the edit distance of the isomirs to all the mirnas in the database to spot problematic isomirs that we should remove
* Repeat plots separating between low/medium/high
* kmer analysis to explain difference among labs

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

