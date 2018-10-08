
# Sequence data organization system

## "assemblies" directory contents

On server 128.32.146.137, there is a /media/walllab/zhanna/owl/assemblies directory.

* StrOccCau_1.0_nuc_finalMito.fa : This the Str_Occ_Cau_1.0 genome assembly.

*

## "filtered_reads" directory contents

On server 128.32.146.137, there is a /media/walllab/zhanna/owl/filtered_reads directory.
There are subdirectories for each sample, which contain the following:

* Adapter-trimmed .fastq.gz files.

* Aligned, merged singletons and paired reads, sorted, de-duplicated .bam files.

    - If the file name contains the string "alnStrOccCau1", then I aligned it to the Str_Occ_Cau_1.0 genome assembly.

    - If the file name contains the string "alnStrOccCau2". then I aligned it to the Str_Occ_Cau_2.0 genome assembly.

* BQSR alignment files.

    - SpBarecal.bam : This extension denotes that I performed BQSR using a set of filtered SNPs from all of the high-coverage spotted and barred owl sequences that I had available.

    - CAspowrecal.bam : This extension denotes that I performed BQSR using a set of filtered SNPs from just the high-coverage spotted owl sequences that I had available (4 samples from California populations).

* g.vcf files. To accompany the .bam files described above, there are various g.vcf files that I created using the GATK HaplotypeCaller. Unless you are re-doing an old analysis, the most recent and best file to use is the one ending in the string "alnStrOccCau2_mates_singles_sorted_dedup.bam_SpBarecal.g.vcf", which denotes that I aligned the sequences to the Str_Occ_Cau_2.0 genome assembly and performed BQSR using the set of all of the high-coverage spotted and barred owl sequences.
