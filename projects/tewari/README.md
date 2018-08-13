This folder contains the results of the re-analysis for tewari data.

# Communication

* BOSC2018: https://f1000research.com/slides/7-953

# Goals

* How to filter False isomirs and True isomirs

Ideally:

* Machine learning to predict the real isomirs and remove the ones due to technical errors


# Limitations

* miRge did the trimming by itself
* bcbio did the trimming by itself and shared the trimmed files to be used by isomiRSEA and sRNAbench
* bcbio has a internal cutoff of a minimum of 2 counts to be annotated
* isomiRSEA considers only iso_5p:+/-1
* sRNAbench labels some sequences as `mv` for isomiRs, these are lost for now in the conversion to GFF3 format.

# reproducibility

The conversion to GFF3 is done by the script:  `scripts/pre_cmd.sh`. This will generate the `expression_cpunts.tsv.gz` files for each tool showed here and script in `R` are the responsible to reproduce the figures:

* `r_code/run_isomir_stats.R` generates the plot to show the amount of reads and sequences for all the samples
* `r_code/run_isomir_commonality.R` generates the plot to show the commonality among replicates for each lab and protocol

# pilot samples

To start the comparison, we decided to reduce the number of samples to analyze to give us an idea of the difference between labs and tools.

The samples selected are: https://github.com/miRTop/incubator/blob/master/projects/tewari/meta_pilot.csv

These include 4 labs, and 2 protocols, 4 replicates per lab and protocol.

# Interesting links

* https://www.ncbi.nlm.nih.gov/pubmed/22647250 (Reducing ligation bias of small RNAs in libraries for next generation sequencing)
* http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0126049 (Bias in Ligation-Based Small RNA Sequencing Library)
* https://projecteuclid.org/download/pdfview_1/euclid.aoas/1318514284 (IDR paper)

# Results

Questions to answer here are:

* Inside the same tool, how replicates from each lab and protocol reproduce by isomiR type
* Inside the same tool, how replicates from each lab but different protocols reproduce by isomiR type
* Inisde each lab and protocol, how tools reproduce each other by isomiR type

## For each tool, general stats of each lab/protocol

Look to individual stats figures at [stats](md/stats.md)


###  Number of miRNA reads

![summary_sum](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/summary_sum.png)

### Number of unique sequences

![summary_counts](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/stats/summary_counts.png)


## For each tool, commonality for replicates inside each lab and protocol

### bcbio

![bcbio_commonality](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/replicates/bcbio.png)

### miRge

![mirge_commonality](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/replicates/mirge.png)

### isomiR-SEA

![isomirsea_commonality](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/replicates/isomirsea.png)

### sRNAbench

![srnabench_commonality](https://github.com/miRTop/incubator/raw/master/projects/tewari/figures/replicates/srnabench.png)

