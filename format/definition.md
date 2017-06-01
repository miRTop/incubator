
As discussed here: https://github.com/miRTop/incubator/issues/10 we'll try to define a GFF3 format for output of small RNA pipelines

authors: @lpantano @gurgese @ThomasDesvignes @mhalushka @mlhack @keilbeck @BastianFromm @ivlachos @TJU-CMC 

GFF3 definition: http://gmod.org/wiki/GFF3

**Goals**

We'll focus on miRNA first, then export to other databases

* define each column to be adapted to miRNA annotation
* define column 9 to support informatino like count data for samples, isomiRs labeling etc

Note: Keep in mind this is for the output of a pipeline, so we know there will be bias toward methodology, but the idea is to put enough information to be able to re-analyze or filter sequences using information described here. As well, it would be a proxy for downstream analysis or packages.

**Deadline**: https://github.com/miRTop/incubator/milestone/2

**Description**

Please add description for each columnd/attribute

* header:
  * commands used to generate the file. At least information about adapter removal and filtering
  * genome version used (maybe try to get from BAM file if GFF3 generated from it)
  * sample names used in attribute:Expression
* column1: seqID:
* column2: source:
* column3: type:
* column4/5: start/end:
* column6: score:
* column7: strand:
* column8: phase: (For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame)
* column9: attributes
  * ID: unique ID based on sequence like mintmap has for tRNA: prefix-22-BZBZOS4Y1 (https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates). good way to use it as cross-mapper ID between different naming or future changes.
  * Name:
  * Parent: hairpin precursor name
  * Alias: get names from miRBase/miRgeneDB
  * Expression: raw counts separated by `,`
  * Filter: PASS or REJECT (this allow to keep all the data and select the one you really want to conside as valid features)
  
  
**API**

Developing an API to check the format and help with the following feature will help to get people using it. Features should have:

* command line and API
* from BAM file and miRNA annotation files create the GFF3 file
* parse GFF3 to filter, get stats, or convert to tabular format. something similar to samtools/vcftools, etc...
* get fasta sequence from GFF3 attribute and genome FASTA file
* (Add more)

I think python is a good way, installing via pipy or conda should get people using it pretty quickly. 

Repo: https://github.com/miRTop/mirtop (right now is adapted to produce tab format but the idea is to migrate to this)
