## Goals

* Optimize GFF format definiation and usability
* Detect methodology accuracy due to tools and some experimental step in the protocols.


Project road map: https://github.com/miRTop/incubator/projects/2

## Data

http://www.biorxiv.org/content/biorxiv/early/2017/05/17/113050.full.pdf

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=whipakmajrwprcv&acc=GSE94586

## Tools

* bcbio smallRNA-seq pipeline + isomiRs - On charge Lorena Pantano
* isomiR-SEA - On charge Gianvito Urgese
* ChimiRa, miRge - On charge Marck Halushka
* sRNAbench - On charge Michael Hackenberg
* Prost - Thomas Desvignes
* miRGe - Marc K. Halushka
* (Add your tool here and person will do it)

## Questions to address

* [To be added]

## Milestones:

### Set up

* [X] Select random public data
* [X] Run with all the tools listed above
* [X] Put data in common space
* [X] Adapt output tools to GFF format

### Random sample

Sample [SRR5756178](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5756178) is a whole blood small RNA-seq run from this manuscript https://academic.oup.com/nar/article/4080663 and is part of project PRJNA391912.  It has ~ 2.8 million reads, of which ~2.6 million are miRNAs.

### Synthetic data

Benchmark was done with synthetic isomiRs for one human miRNA, see [results](https://github.com/miRTop/incubator/tree/master/synthetic).

### Tewari data

* open an issue to ask for access if you are a new collaborator
* adapters for each sample is [Adaptors_to_remove.txt](Adaptors_to_remove.txt).
* [pilot config data](tewari_mini.csv) used to test adapters are ok
