# microRNA gene and products Nomenclature (created by Thomas Desvignes)

The situation:
Clear unambiguous microRNA naming is essential for accurate communication. Many ways have been used over the years to name microRNA genes, their products but also their intermediate forms (primary transcript and precursor). miRBase, the international reference microRNA database, Gene Nomenclature Consortia (HGNC, MGI, ZFIN, etc...) used to use different nomenclature guidelines.

A unified nomenclature system is needed to improve communication.

The optimal system should:
  - be applicable to all metazoan species 
  - be consistent with each species' gene nomenclature guidelines (HGNC for Human, MGI for Mouse, ZFIN for Fish, etc) 
  - be evolutionarily robust
  
A unified system for microRNA naming would permit more efficient communication between researchers within a species or across species because the naming would be consistent and follow familiar guidelines used for protein-coding genes.
  
*Note:*
*Plants miRNA research is facing the same problems. The different origin of plant and metazoan miRNA systems, however, make the creation of a unified “plant &amp; metazoan” pattern hazardous. Nevertheless, tentative standardization on specific aspects is a future goal of the group.*

*Budak, H., Bulut, R., Kantar, M., &amp; Alptekin, B. (2015). MicroRNA nomenclature and the need for a revised naming prescription. Briefings in Functional Genomics, elv026–. <a href="http://doi.org/10.1093/bfgp/elv026">http://doi.org/10.1093/bfgp/elv026</a>*

#Comments from Thomas Desvignes:
Two recent articles from two independent groups have simultaneously proposed novel nomenclature guidelines to improve microRNA naming:
Desvignes, T., Batzel, P., Berezikov, E., Eilbeck, K., Eppig, J. T., McAndrews, M. S., Singer, A., Postlethwait, J. H. (2015). miRNA Nomenclature: A View Incorporating Genetic Origins, Biosynthetic Pathways, and Sequence Variants. Trends in Genetics : TIG, 31(11), 613–626. <a href="http://doi.org/10.1016/j.tig.2015.09.002">http://doi.org/10.1016/j.tig.2015.09.002</a>

Fromm, B., Billipp, T., Peck, L. E., Johansen, M., Tarver, J. E., King, B. L., Newcomb J., Sempere L. F., Flatmark K., Hovig E., Peterson, K. J. (2015). A Uniform System for the Annotation of Vertebrate microRNA Genes and the Evolution of the Human microRNAome. Annual Review of Genetics. <a href="http://doi.org/10.1146/annurev-genet-120213-092023">http://doi.org/10.1146/annurev-genet-120213-092023</a> (<a href="http://mirgenedb.org">http://mirgenedb.org</a>)</li>

One article (Fromm et al, 2015) focused on evolutionary conservation and clear orthology establishment to provide phylogenetically informative names based on Human annotation, while the other (Desvignes et al, 2015) envisioned a system for sequence variation analysis following nomenclature guidelines in use for protein-coding genes. 

The first system (Fromm et al, 2015) has the advantage of establishing clear orthologies between genes across species and rename orthologs in a similar way in all species and harmonizes names within a miRNA family. Some problems, however, exist with this system: 
  1) a lot of names are modified. A lot of researchers have been now used to specific microRNA numbers for years. Changing many names, even though clarifying evolutionary links, can be confusing and is unlikely to be followed.
  2) the nomenclature system chosen doesn't follow any of the gene nomenclature guidelines established by Nomenclature consortia. miRNAs are also genes and should thus follow similar guidelines.
  3) the complete system is Human centric and poses problem because Human should not be envisioned as the end or start point in evolutionary studies, and also because this human-centered system is not evolutionarily robust with genome duplication events such as the two rounds of whole genome duplication that occurred at the stem of the vertebrate radiation, or the teleost specific whole genome duplication event.

The second system (Desvignes et al, 2015) on the other hand, doesn't changes the number constituting the name but encourages researchers, when identifying new miRNAs, to perform accurate phylogenetic verifications to provide phylogenetically conform names. In addition, this nomenclature system has been elaborated in collaboration with Mouse and Zebrafish Gene nomenclature consortia (MGI and ZFIN respectively), and is now the system in use in these gene nomenclature authorities. It is also in agreement with the guidelines in use for Human and other animal models. Finally, the phylogenetic relationship between orthologs or co-orthologs are dealt similarly to protein-coding gene and is thus a system already used by everyone, and can be phylogenetically accurate when appropriately used. 

Because Gene Nomenclature Consortia use of the nomenclature system proposed by Desvignes et al, we propose that this nomenclature system be the one used a basis.

#Star or not star?
The use of the start strand denomination '*' is not approved by any nomenclature consortium, including miRBase. Nonetheless, for cases of extreme difference in levels of expression, this additional symbol can provide potentially useful functional information. If people wants to use it informally, I would not see any objections.
But because of the risk of fluctuation, which was the motivation of its removal in first place, an official use, from my point of view, can't be promoted. Indeed, the denomination would rely too much on sequencing depth and the nature of the tissue studied. In addition, at what level of arm selection can we say that a strand is likely only a by-product? Fromm et al propose a two-fold change. This appears to me not strong enough of a difference to call the second strand star strand, given the non-representation of the complete expressed miRNAOme of an organism and the sequence bias known in miRNA-Seq library preparation. I would suggest at least a 10-fold change and a good representation of tissue types in the organism considered.

#The microRNA genes and products agreement:
<pre>
           Symbol	                 Signifies
Human	        Mouse	      Zebrafish	
MIRX	        MirX	      mirX	         microRNA gene X
MIRXy	        MirXy	      mirXy	         yth member of the mirX family (y can be a letter or a "dot-number")
MIRXOS          MirXos	      mirXos	     mirror microRNA gene of the microRNA gene X 
pri-MIRX	    pri-MiRX     pri-miRX	     primary microRNA transcript X
pre-MIRX	    pre-MiRX	  pre-miRX	     microRNA precursor X
MIRX-5p	        MIRX-5p	      MiRX-5p	5'     arm mature miRNA sequence 
MIRX-3p	        MIRX-3p	      MiRX-3p	3'     arm mature miRNA sequence 
MORX-5p	        MORX-5p	      MoRX-5p	5'     arm mature moR sequence 
MORX-3p	        MORX-3p	      MoRX-3p	3'     arm mature moR sequence 
LORX	        LORX	      LoRX	         loop fragment of the miRNA precursor
MIRCX	        MircX	      mircX	         mir gene cluster X
</pre>

If you want to add the nomenclature system for another species, please do so.
