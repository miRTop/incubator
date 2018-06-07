This folder contains the results of the re-analysis for tewari data.

# Limitations

* miRge did the trimming by itself
* bcbio did the trimming by itself and shared the trimmed files to be used by isomiRSEA and sRNAbench
* bcbio has a internal cutoff of a minimum of 2 counts to be annotated
* isomiRSEA considers only iso_5p:+/-1
* sRNAbench labels some sequences as `mv` for isomiRs, these are lost for now in the conversion to GFF3 format.


# pilot samples

To start the comparison, we decided to reduce the number of samples to analyze to give us an idea of the difference between labs and tools.

The samples selected are: https://github.com/miRTop/incubator/blob/master/projects/tewari/meta_pilot.csv

These include 4 labs, and 2 protocols, 4 replicates per lab and protocol.


Questions to answer here are:

* Inside the same tool, how replicates from each lab and protocol reproduce by isomiR type
* Inside the same tool, how replicates from each lab but different protocols reproduce by isomiR type
* Inisde each lab and protocol, how tools reproduce each other by isomiR type

## For each tool, similarity between replicates of each lab/protocol


###  Number of differents isomiRs for each lab for protocol TrueSeq:

![bcbio_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/bcbio_truseq_count.png)
![mirge_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/miRge_truseq_count.png)
![isomirsea_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/isomiR-SEA_truseq_count.png)
![srnabench_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/sRNAbench_truseq_count.png)

## Number of differents isomiRs for each lab for protocol NEBNext:

![bcbio_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/bcbio_nebnext_count.png)
![mirge_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/miRge_nebnext_count.png)
![isomirsea_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/isomiR-SEA_nebnext_count.png)
![srnabench_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/sRNAbench_nebnext_count.png)


