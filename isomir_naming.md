### Annotation agreement

This section will have the result of the discussion for the naming.

In miRNAs, the variations that can be observed are (highlighted the most common): 

-	3’ modifications:
  -	**Nucleotide deletions**
  - **Non-templated nucleotide additions**
  - **Templated nucleotide additions**
-	5’ modifications:
  - **Nucleotide deletions (-> Seed-shift)**
  - Non-templated nucleotide additions (-> Seed-shift) **can this happen**
  - Nucleotide deletions compensated by Non-templated additions (nb: we can only analyze this case, if there is deletion followed by templated addition, there is no visible difference with the non-modified 5’ miRNA)
  - **Templated nucleotide additions (-> Seed-shift)**
-	Editions (they are all the same but we can differentiate 3 types based on their putative effect on function):
  - **In the seed**
  - **In the 3’ complementary region**
- Outside the seed or 3’CR
  -	Insertions and Deletions: I personally never observed such modifications but technically they are possible.
  -	Combinations

Accepted rules:
- deletions and additions will be lower and capital case respectively.
- adding to miR name the modification one byt one. Need to agree in character separation.

** Naming in papers**

miR name plus modification added after that (need discussion an examples)

**File format**

tab format with following columns (need discussion):

- sequence
- mirna name
- precursor
- isomir TAGs in the sequence

## Discussion

### From Lorena

Right now, I have tags for each type of modification I found:

**reference_miRNA_name:subs_tag:addition_tag:trimming_start_tag:trimming_end_tag**

it would be like:

`hsa-let7a-5p:0:AA:0:at`

0 would mean no changes, capital letters will be insertions compared to reference, lower letters will be deletions. 
Nucleotides at the addition_tag DONOT match the precursor, the trimming_tag MATCH the precursor.


### From Thomas and John:

There could be many ways to refer to an isomiR compared to his RefSeq-miRNA.
Our goal when we proposed the unified nomenclature of miRNA was to match what is already in use for protein-coding genes. 
Furthermore, we made the decision to use protein nomenclature guidelines when referring to the mature miRNAs because 
the mature miRNAs are the final products, like proteins are the final products of protein-coding genes. 

Currently, there is basically only one rule for proteins mutations description and refers to precise modification at specific location. 
We can differentiate two types:
-	Missense mutations resulting from a point mutation creating a codon coding for s different amino acid: e.g. R527L, 
which means that an original Arginine (R) at position 527 of the protein is changed to a Glutamate (L).
-	Nonsense mutations resulting from a point mutation creating a stop codon truncating the protein: e.g. R527*,
which means that an original Arginine (R) at position 527 of the protein is changed to a stop codon (*)


The edited variants are the easy ones to describe because we can use the system in use for proteins, 
e.g A13G to refer to an isomiR that display an edited base (A to I) at position 13.

I like the idea of using lower case vs capital letters to differentiate deletions and additions, respectively.

Because 5’ modifications inducing a seed-shift are functionally important ones,
I wonder if we should include in the name something in reference to that, like a simple “S” or “s” symbol standing for “seed-shift”.
Also, this is important because by changing the start of the sequence all the positions along the miRNA get changed.
For example, if there is a common A to G edition at position 7 of the RefSeq-miRNA (edition case below) then in the isomiR
with a 5’ end non-templated addition, the edition site becomes position 8. 
A way to deal with that could be to keep fixed the positions based on the RefSeq-miRNA and a 5’ addition could be “AsT”
if there is an untemplated T addition, or “as” if there is a a deletion?

For 3’ end shortening, because they may be quite shorter, one could also use the protein denomination and for example,
if the isomiR is 18nt instead of 21, name it something like “miRA-5p-19*” which would make the naming simpler than “miRA-5p-a19-a20-a21”.

Here are some examples of a system that can be used to refer to a specific isomiR within a text or article.
I think this system is informative and not too heavy but we probably need to think about it a bit more: 

```
Description	Sequence	Name
DNA	                          AAAAACAAAAAAAAAAAAAAAAAAA	mirA
RefSeq-miRNA	                  AAACAAAAAAAAAAAAAAAAA	miRA-5p
3' end deletion	                AAACAAAAAAAAAAAAAAAA	miRA-5p-a21
3' end non-templated addition	  AAACAAAAAAAAAAAAAAAAAT	miRA-5p-A22T
3' end templated addition	      AAACAAAAAAAAAAAAAAAAAA	miRA-5p-A22
5' end deletion	                 AACAAAAAAAAAAAAAAAAA	miRA-5p-as
5' end non-templated addition	TAAACAAAAAAAAAAAAAAAAA	miRA-5p-AsT
5' end templated addition	    AAAACAAAAAAAAAAAAAAAAA	miRA-5p-As
Seed Edition	                 AAACAAGAAAAAAAAAAAAAA	miRA-5p-A7G
3' CR edition	                 AAACAAAAAAAAAGAAAAAAA	miRA-5p-A14G
Other edition	                 AAACAAAAAAGAAAAAAAAAA	miRA-5p-A11G
Combination Ex1	               AAACAAGAAAAAAAAAAAAA	miRA-5p-A7G-a21
Combination Ex2	               AAACAAAAAAAAAGAAAAAAAT	miRA-5p-A14G-A22T
Combination Ex3	              TAAACAAGAAAAAAAAAAAAAAT	miRA-5p-AsT-A7G-A22T
```
