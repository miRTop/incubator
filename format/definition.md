

[![Project Status: Moved to http://example.com – The project has been moved to a new location, and the version at that location should be considered authoritative.](http://www.repostatus.org/badges/latest/moved.svg)](http://www.repostatus.org/#moved) to [https://github.com/miRTop/mirGFF3](https://github.com/miRTop/mirGFF3)


As discussed here: https://github.com/miRTop/incubator/issues/10 we'll try to define a GFF3 format for output of small RNA pipelines

authors: @lpantano @gurgese @ThomasDesvignes @mhalushka @mlhack @keilbeck @BastianFromm @ivlachos @TJU-CMC @sbb25 

GFF3 definition: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

**Goals**

We'll focus on miRNA first, then export to other databases

* define each column to be adapted to miRNA annotation
* define column 9 to support information like count data for samples, isomiRs labeling etc

Note: Keep in mind this is for the output of a pipeline, so we know there will be bias toward methodology, but the idea is to put enough information to be able to re-analyze or filter sequences using information described here. As well, it would be a proxy for downstream analysis or packages.

**Deadline**: https://github.com/miRTop/incubator/milestone/2

**Description**

Please add description for each columnd/attribute

* header:
  * database: `##source-ontology LINK TO DATABASE` include version and link (mandatory)
  * commands used to generate the file. At least information about adapter removal, filtering, aligner, mirna tool. All of them starting like: `## CMD: `. Can be multiple lines starting with this tag. (optional)
  * genome/database version used (maybe try to get from BAM file if GFF3 generated from it): `## REFERENCE:` 
  * sample names used in attribute:Expression: `## COLDATA:` separated by commas (mandatory)
  * small RNA GFF version `## VERSION: 0.9`
  * Filter tags meaning: See Filter attribute below. Different filter tags should be separated by `,` character. Example: `## FILTER: ` and example would be `## FILTER: PASS(is ok), REJECT(false positive), REJECT lowcount(rejected due to low count in data)`. 
* column1: seqID: precursor name
* column2: source: databases (lower case) used for the annotation (miRBase, mirgeneDB,tRNA...etc): https://github.com/miRTop/incubator/issues/13. With the version number after `_` character: `mirbase_21`
* column3: type: `ref_miRNA, isomiR`: https://github.com/miRTop/incubator/issues/13  (SO:0002166 ref_miRNA and SO:0002167 isomiR)
* column4/5: start/end: precursor start/end as indicated by alignment tool
* column6: score (Optional): It can be the mapping score or any other score the tool wants to assign to the sequence.
* column7: strand. In the case of mapping against precursor should be always `+`. It should accept mapping against the genome: `+/-` allowed.
* column8: phase: (For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame): Not relevant righ now. This can be: `.`
* column9: attributes: 
  * UID: unique ID based on sequence like mintmap has for tRNA: prefix-22-BZBZOS4Y1 (https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates). good way to use it as cross-mapper ID between different naming or future changes. The tool will implement this, so an API can be used to fill this field. Currently supported by [mirtop](https://github.com/miRTop/mirtop/blob/dev/mirtop/mirna/realign.py#L124) code.
  * Read: read name
  * Name: mature name
  * Parent: hairpin precursor name
  * Variant: (categorical types - adapted from isomiR-SEA)
    * iso_5p:+/-N. `+` indicates extra nucleotides not in the reference miRNA. `-` indicates removed nucleotides not in the sequence. `N` the number of nucleotides of difference. For instance, if the sequence starts 2 nts after the reference miRNA, the label will be: `iso_5p:-2`, but if it starts before, the label will be `iso_5p:+2`.
    * iso_3p:+/-N. Same explanation applied.
    * iso_add:+N. Same explanation applied.
    * iso_snp_seed: when affected nucleotides are betweem [2-7].
    * iso_snp_central_offset: when affected nucleotides is at position [8].
    * iso_snp_central: when affected nucleotides are betweem [9-12].
    * iso_snp_central_supp: when affected nucleotides are betweem [13-17].
    * iso_snp: anything else.
  * Changes (optional): similar to previous one but indicating the nucleotides being changed. 
    * additions are in capital case
    * deletions are in lower case
    * exaple: `Changes iso_5p:0,iso_3p:TT,iso_add:GTC,iso_snp:0` where `Variant iso_add:+3,iso_3p:-2`.
  * Cigar: CIGAR string as indicated [here](https://samtools.github.io/hts-specs/SAMv1.pdf). It is the standard CIGAR for aligners. With the restriction that `M` means exact match always. That's a difference with some aligners where `M` includes mismatches. In this case, if there is a mismatch, then it should be output like: `11MA7M` to indicates there is a mismatch at position 12, where `A` is the reference nucleotide.
  * Hits: number of hits in the database.
  * Alias (Optional): get names from miRBase/miRgeneDB or other database separated by `,`
  * Genomic (Optional): positions on the genome in the following format: `chr:start-end,chr:start-end`
  * Expression: raw counts separated by `,`. It should be in the same order than `colData` in the header.
  * Filter: PASS or REJECT (this allow to keep all the data and select the one you really want to conside as valid features). PASS can have subclases: `PASS:te`: meaning the sequence pass but the tools consider variants showed here are not trusted. REJECT can go with any short word explaining why it was rejected: `REJECT:lowcounts`. In this case the sequence will be skipped for data mining of the file when quering counts or summarize miRNA expression.
  * Seed_fam (Optional): in the format of 2-8 nts and reference miRNA sharing the seed. Usefull to go for pre-computed target predictions: `ATGCTGT:mir34a_5p`
  
**API**

Developing an API to check the format and help with the following feature will help to get people using it. Features should have:

* command line and API
* from BAM file (precursor/genome) and miRNA annotation files create the GFF3 file
* parse GFF3 to filter, get stats, or convert to tabular format. something similar to samtools/vcftools, etc...
* get fasta sequence from GFF3 attribute and genome FASTA file
* (Add more)

I think python is a good way, installing via pipy or conda should get people using it pretty quickly. 

Repo: https://github.com/miRTop/mirtop
