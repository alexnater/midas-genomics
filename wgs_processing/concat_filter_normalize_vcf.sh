#!/bin/bash


### manually check breakpoints of vcf files and make sure that there are no overlapping variants or truncated file endings:
FILES=$(find all_Midas_withCenSiqTFA.raw.snps.indels.*.vcf | sort -k1V)

for FILE in $FILES; do echo $FILE; grep -v "^#" $FILE | head -1 | cut -f 1-5; tail -1 $FILE | cut -f 1-5; tail -1 $FILE | awk '{print NF}'; done \
	> vcf_breakpoints.txt 


### concatenate blocks into single vcf file:
vcf-concat $FILES > all_Midas_withCenSiqTFA.raw.snps.indels.vcf


### sort according to order in freebayes_all_Midas_withCenSiqTFA.poplist.txt file:
vcf-shuffle-cols -t header_VCF_sorted.txt all_Midas_withCenSiqTFA.raw.snps.indels.vcf  \
	> all_Midas_withCenSiqTFA.raw.snps.indels.sorted.vcf


### filter out poorly-supported or otherwise biased variant sites:
vcffilter -s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" all_Midas_withCenSiqTFA.raw.snps.indels.sorted.vcf \
	> all_Midas_withCenSiqTFA.sitefilt.snps.indels.vcf


### normalize variant representation (left aligned and no superfluous bases represented) and remove potentially duplicated variants (just in case):
vt normalize all_Midas_withCenSiqTFA.sitefilt.snps.indels.vcf \
-r ../../../reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.fasta \
	| vt uniq - > all_Midas_withCenSiqTFA.sitefilt.snps.indels.norm.vcf 

