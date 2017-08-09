
As discussed here: https://github.com/miRTop/incubator/issues/10 we'll try to define a GFF3 format for output of small RNA pipelines

authors: @lpantano @gurgese @ThomasDesvignes @mhalushka @mlhack @keilbeck @BastianFromm @ivlachos @TJU-CMC 

GFF3 definition: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

**Goals**

We'll focus on miRNA first, then export to other databases

* define each column to be adapted to miRNA annotation
* define column 9 to support informatino like count data for samples, isomiRs labeling etc

Note: Keep in mind this is for the output of a pipeline, so we know there will be bias toward methodology, but the idea is to put enough information to be able to re-analyze or filter sequences using information described here. As well, it would be a proxy for downstream analysis or packages.

**Deadline**: https://github.com/miRTop/incubator/milestone/2

**Description**

Please add description for each columnd/attribute

* header:
  * database: `##source-ontology LINK TO DATABASE` include version and link
  * commands used to generate the file. At least information about adapter removal, filtering, aligner, mirna tool. All of them starting like: `## CMD: `
  * genome version used (maybe try to get from BAM file if GFF3 generated from it)
  * sample names used in attribute:Expression: `## colData:` separated by spaces
  * small RNA GFF version `## version: 0.9`
* column1: seqID: precursor name
* column2: source: databases (lower case) used for the annotation (miRBase, mirDBgene,tRNA...etc): https://github.com/miRTop/incubator/issues/13. With the version number after `_` character: `mirbase_21`
* column3: type: `ref_miRNA, isomiR`: https://github.com/miRTop/incubator/issues/13  (SO:0002166 ref_miRNA and SO:0002167 isomiR)
* column4/5: start/end: precursor start/end as indicated by alignment tool
* column6: score:
* column7: strand
* column8: phase: (For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame)
* column9: attributes: 
  * ID: unique ID based on sequence like mintmap has for tRNA: prefix-22-BZBZOS4Y1 (https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates). good way to use it as cross-mapper ID between different naming or future changes. The tool will implement this, so an API can be used to fill this field.
  * Name: mature name
  * Parent: hairpin precursor name
  * Variant: categorical types: iso_5p, iso_3p, iso_snp(_seed/_central_supp), iso_add (adapted from isomiR-SEA)
  * Cigar: CIGAR string as indicated here: []
  * Alias: get names from miRBase/miRgeneDB or other database separated by `,`
  * Genomic: positions on the genome in the following format: `chr:start-end,chr:start-end`
  * Expression: raw counts separated by `,`. It should be in the same order than `colData` in the header.
  * Filter: PASS or REJECT (this allow to keep all the data and select the one you really want to conside as valid features). PASS can have subclases: `PASS:te`: meaning the sequence pass but the tools consider variants showed here are not trusted. REJECT can go with any short word explaining why it was rejected: `REJECT:lowcounts`. In this case the sequence will be skipped for data mining of the file when quering counts or summarize miRNA expression.
  * Seed_fam: in the format of 2-8 nts and reference miRNA sharing the seed. Usefull to go for pre-computed target predictions: `ATGCTGT:mir34a_5p`
  
**API**

Developing an API to check the format and help with the following feature will help to get people using it. Features should have:

* command line and API
* from BAM file and miRNA annotation files create the GFF3 file
* parse GFF3 to filter, get stats, or convert to tabular format. something similar to samtools/vcftools, etc...
* get fasta sequence from GFF3 attribute and genome FASTA file
* (Add more)

I think python is a good way, installing via pipy or conda should get people using it pretty quickly. 

Repo: https://github.com/miRTop/mirtop (right now is adapted to produce tab format but the idea is to migrate to this)
