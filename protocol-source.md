---
title: Supplementary methodology for fine scale recombination rates in barred owls (*Strix varia*) from eastern North America
author: Zachary R. Hanna
author: Jeffrey D. Wall
affiliation: Institute for Human Genetics, University of California, San Francisco, San Francisco, California, United States of America
---

## Repeat masking of genome assembly

We performed a homology-based annotation of repetitive regions using RepeatMasker version open-4.0.7 [@smitRepeatMaskerOpen4Accessed2013a] with the repeat databases of the DFAM library version 2.0 [@hubleyDfamDatabaseRepetitive2016] and the Repbase RepeatMasker libraries version 20170127 [@baoRepbaseUpdateDatabase2015; @jurkaRepbaseUpdateDatabase2005; @jurkaRepeatsGenomicDNA1998; @jurkaRepbaseUpdateDatabase2000]. Our installation of the RepeatMasker program relied on tandem repeats finder (TRF) version 4.09 [@bensonTandemRepeatsFinder1999] in addition to the NCBI BLAST+ version 2.7.1 [@camachoBLASTArchitectureApplications2009], RMBlast version 2.2.28 [@smitRMBlastAccessed20162015], and HMMER version 3.1b2 (@eddyProfileHiddenMarkov1998; http://hmmer.org) sequence search tools. We first performed a homology-based annotation of repetitive regions in the reference genome using RepeatMasker with the options '-gccalc -species aves'.

```
RepeatMasker -pa 12 -gccalc -species aves reference_genome.fa 1>reference_genome_RMask_aves.log 2>reference_genome_RMask_aves.err
```

* We created a de novo model of the repeats in the reference genome using RepeatModeler version 1.0.8 [@smitRepeatModelerOpen1Accessed2015a]. Our RepeatModeler installation used the RepeatMasker version open4.0.7 installation as detailed in the previous step with Repbase RepeatMasker libraries version 20170127, the RMBlast version 2.2.28 sequence search tool, TRF version 4.09, and two de novo repeat finding tools, RECON version 1.08 [@baoAutomatedNovoIdentification2002] and RepeatScout version 1.0.5 [@priceNovoIdentificationRepeat2005]. We first used the RepeatMasker BuildDatabase tool to build a sequence database from the reference nuclear genome, StrOccCau_2.0_nuc.fa (we didn't want the model to include repeats from the mitochondrial genome).

```
BuildDatabase -name nuclear_genome nuclear_genome.fa 1>nuclear_genome_RModbuilddb.log 2>nuclear_genome_RModbuilddb.err
```

* We then ran RepeatModeler with default options.

```
RepeatModeler -pa 6 -database nuclear_genome 1>nuclear_genome_RMod.log 2>nuclear_genome_RMod.err
```

* We performed a final round of repeat masking by supplying the masked genome output from the homology-based masking to RepeatMasker with the options "-gccalc -lib <RepeatModeler_output>", where the "RepeatModeler_output" was the "consensi.fa.classified" repeat library file output from our RepeatModeler run.

```
RepeatMasker -pa 12 -gccalc -lib consensi.fa.classified reference_genome.fa.masked.fa 1>RMask.log 2>RMask.err
```

* Since RepeatMasker masked all of the repetitive regions with N characters, we used seqtk version 1.2-r94 [@liSeqtkVersion2r942016] with options "cutN -n1 -p100000000 -g" to create a browser extensible data (BED) formatted file of the coordinates of all of the N regions (these included the gap regions in the original reference genome).

```
seqtk cutN -n1 -p100000000 -g reference_genome.fa.masked.fa.masked >reference_genome.fa.masked.fa.masked_Nregions.bed
```

* We then used the GNU core utility sort version 8.25 [@haertelSortGNUCoreutils2016] with options '-k1,1 -k2,2n' to sort the BED file.

```
sort -k1,1 -k2,2n reference_genome.fa.masked.fa.masked_Nregions.bed >reference_genome.fa.masked.fa.masked_Nregions_sorted.bed
```

## Sequence data

### Trimming
We used Trimmomatic version 0.36 [@bolgerTrimmomaticFlexibleTrimmer2014]
to trim adapter sequences.
In the Trimmomatic command, we employed the options below that included the use of an adapter file.
For all sequences apart from those for CAS:ORN:98821, we used the “TruSeq3-PE-2.fa” adapter file supplied in Trimmomatic version 0.36.
For the CAS:ORN:98821 sequences, we used the adapter sequences appropriate to each set of sequences [@hannaAdapterSequencesUsed2018; @hannaNorthernSpottedOwl2017].

```
java -jar trimmomatic-0.36.jar PE -threads 48 -trimlog <trimlog_file_path> <R1_input_file> <R2_input_file> <R1_paired_output_file> <R1_unpaired_output_file> <R2_paired_output_file> <R2_unpaired_output_file> ILLUMINACLIP:adapter_file.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
```

We trimmed the CAS:ORN:98821 sequences using the [trim_Sequoia.py](trim_Sequoia.py) wrapper script.

We trimmed the other sequences using [trim_high_cov_owls.py](trim_high_cov_owls.py).

### Reference genome

We combined the following using cat (GNU core utilities) version 8.25 [@granlundCatGNUCoreutils2017] in order to generate the reference genome:

1. StrOccCau_2.0_nuc.fa

2. mitochondrial genome [@hannaCompleteMitochondrialGenome2017]

```
cat StrOccCau_2.0_nuc.fa mitochondrial_genome.fa >reference_genome.fa
```

### Alignment, processing, and individual sample variant calls

For the sample CAS:ORN:98821 sequence sets, we aligned, processed the alignments, and called variants using the [aln_Sequoia.py](aln_Sequoia.py) wrapper script in this repository. The steps are very similar to those we followed for the other samples with the main difference being that in the paired-end alignment step we used a maximum insert size parameter (-w) appropriate to each sequence set.

For all other samples, we aligned, processed the alignments, and called variants using the [aln_high_cov_owls.py](aln_high_cov_owls.py) wrapper script in this repository.

The general commands we provide below in each section are specific to the samples other than CAS:ORN:98821. Please refer to the wrapper scripts we have provided for the exact manner in which we executed these commands for each sample.

#### Paired reads

Other than sample CAS:ORN:98821, for each sample, we used BWA-MEM version 0.7.17-r1188 [@liAligningSequenceReads2013] to align the trimmed paired reads to the reference genome.
We employed default options other than parameters 'bwa mem -M'. We used 420 nt as the insert size and set the maximum insert size (parameter '-w') equal to 1000.

As there were multiple read sets from the different SRA run accessions for sample CAS:ORN:98821, we aligned the paired reads for each of these independently.

```
bwa mem -M -t 12 -R '@RG\tID:<sample_run_ln>1\tSM:<sample>\tPL:illumina\tLB:<run_lib>\tPU:lane1\tPI:420' -w 1000 <ref_genome.fa> <R1_trimmed_file> <R2_trimmed_file> ><alignment_output_file 2><alignment_file.log>
```

#### Unpaired reads
We used BWA-MEM version 0.7.17-r1188 [@liAligningSequenceReads2013] to align the trimmed unpaired reads to the reference genome.
We used the default options other than parameters 'bwa mem -M'.

We independently aligned each of the multiple read sets from the different SRA run accessions for sample CAS:ORN:98821.

```
zcat <R1_unpaired_trimmed_file> <R2_unpaired_trimmed_file> | bwa mem -M -t 12 -R '@RG\tID:<sample_run_ln>1\tSM:<sample>\tPL:illumina\tLB:<run_lib>\tPU:lane1\tPI:420' <ref_genome.fa> - ><alignment_output_file 2><alignment_file.log>
```

### Alignment processing

#### Merge alignments

We merged the paired-end and unpaired sequence alignments for each sample using the Picard version 2.17.6 tool MergeSamFiles (http://broadinstitute.github.io/picard; Accessed 2018 Mar 22) with default settings.

We merged all of the independent read set alignments for sample CAS:ORN:98821 together at this step.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar picard.jar MergeSamFiles MAX_RECORDS_IN_RAM=200000 INPUT=<paired_alignment_file.sam> INPUT=<unpaired_alignment_file.sam> OUTPUT=<merged_alignment_file.sam> SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=</temporary/directory/path> 1><output_file.log> 2><output_file.err>
```

#### Sort alignments

For each sample, we then sorted the aligned sequences using the Picard version 2.17.6 function SortSam (http://broadinstitute.github.io/picard) with default settings.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar picard.jar SortSam MAX_RECORDS_IN_RAM=200000 INPUT=<merged_alignment_file.sam> OUTPUT=<merged_sorted_alignment_file.bam> SORT_ORDER=coordinate TMP_DIR=</temporary/directory/path> 1><output_file.log> 2><output_file.err>
```

#### Mark duplicates

We marked duplicate sequences (both PCR and optical) for each sample using the Picard version 2.17.6 MarkDuplicates tool (http://broadinstitute.github.io/picard) with default settings.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar picard.jar MarkDuplicates MAX_RECORDS_IN_RAM=200000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=</temporary/directory/path> INPUT=<merged_sorted_alignment_file.bam> OUTPUT=<merged_sorted_dupmarked_alignment_file.bam> METRICS_FILE=<output_metrics_file> 1><output_file.log> 2><output_file.err>
```

#### Build BAM index file

We indexed the duplicate-marked alignment file for each sample for use with downstream applications using the Picard version 2.17.6 MarkDuplicates tool.

```
java -Xmx10g -Djava.io.tmpdir=/temporary/directory/path -jar picard.jar BuildBamIndex MAX_RECORDS_IN_RAM=200000 TMP_DIR=/temporary/directory/path INPUT=merged_sorted_dupmarked_alignment_file.bam OUTPUT=merged_sorted_dupmarked_alignment_file.bai 1>output_file.log 2>output_file.err
```

#### Individual variant calls

We used the Genome Analysis Toolkit (GATK) version 3.8.0 HaplotypeCaller tool [@depristoFrameworkVariationDiscovery2011; @mckennaGenomeAnalysisToolkit2010; @vanderauweraFastQDataHigh2013] to call variants for each sample individually and create a genomic variant call format file (gVCF).
Our input into HaplotypeCaller for each sample was the corresponding bwa-aligned, sorted, duplicate-marked binary alignment map (BAM) file.
Other than setting “--emitRefConfidence GVCF', we used default options.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference_genome.fa -I <sorted_dedup.bam> --emitRefConfidence GVCF -o output.g.vcf 1>output.g.vcf.log 2>output.g.vcf.err
```

#### Creation of combined sample VCF

In order to combine the individual genomic variant call format file (gVCF) files into a variant call format (VCF) file, we used the GATK version 3.8.0 GenotypeGVCFs tool with default settings.

```
java -Xmx100g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R <reference_genome.fa> --variant <sample1.bam.g.vcf> --variant <sample2.bam.g.vcf> <...other variant files...> --variant <sample29.bam.g.vcf> -o Str2.0_SpBa.vcf 1>Str2.0_SpBa.vcf.log 2>Str2.0_SpBa.vcf.err
```

### SNP filtering

We filtered SNPs using the [filt_SpBa_StrOccCau2_snps.sh](filt_SpBa_StrOccCau2_snps.sh) script. Below we describe each of the steps in further detail.

We followed an adapted version of the GATK guidelines for hard filtering [https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set (Accessed 2018 Mar 15); https://software.broadinstitute.org/gatk/documentation/article.php?id=3225 (Accessed 2018 Mar 15); https://software.broadinstitute.org/gatk/documentation/article.php?id=6925 (Accessed 2018 Mar 15)] to generate a subset of high quality variants to use for base quality score recalibration (BQSR) (https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr; Accessed 2018 Mar 15).

* First, we used the GATK SelectVariants tool with option "-selectType SNP" to extract the single nucleotide polymorphisms (SNPs) from the raw VCF file.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V Str2.0_SpBa.vcf -selectType SNP -o Str2.0_SpBa_snps.vcf
```

* We filtered the SNPs using the GATK VariantFiltration tool with options '--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0"' to filter out variants where the QUAL field value divided by the unfiltered depth of non-homogygous reference samples was less than 2.0, the Phred-scaled probability of stand bias was greater than 60.0, the root mean square mapping quality of all the reads at the site was less than 40.0, the u-based z-approximation from a rank-sum test comparing mapping qualities of reads in favor of the reference versus alternate alleles was less than -12.5, the u-based z-approximation from a rank-sum test comparing the positions of the reference and alternate alleles within reads was less than -8.0, or the strand odds ratio (SOR), which estimates strand bias while accounting for the ratio of reads that cover each of the alleles at a site, was greater than 3.0 (https://software.broadinstitute.org/gatk/documentation/article.php?id=6925; Accessed 2018 Mar 15).

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T VariantFiltration -R <reference_genome.fa> -V Str2.0_SpBa_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName "my_snp_filter" -o Str2.0_SpBa_snps_filt1.vcf
```

* We used the GATK SelectVariants tool with the "--excludeFiltered" option to produce a final VCF file of high-confidence SNPs for BQSR.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V <Str2.0_SpBa_snps_filt1.vcf> --excludeFiltered -o <Str2.0_SpBa_snps_filt2.vcf>
```

### Indel filtering

We filtered indels using the [filt_SpBa_StrOccCau2_indels.sh](filt_SpBa_StrOccCau2_indels.sh) script. Below we describe each of the steps in further detail.

* Similar to the SNPs, we used the GATK SelectVariants tool with option "-selectType INDEL" to extract the indels from the raw VCF file.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V Str2.0_SpBa.vcf -selectType INDEL -o Str2.0_SpBa_indels.vcf
```

* We filtered the indels using the GATK VariantFiltration tools with options '--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"' to filter out variants where the QUAL field value, which measured variant confidence, divided by the unfiltered depth of non-homogygous reference samples was less than 2.0, the Phred-scaled probability of stand bias was greater than 200.0, the u-based z-approximation from a rank-sum test comparing the positions of the reference and alternate alleles within reads was less than -20.0, or the strand odds ratio was greater than 10.0 (https://software.broadinstitute.org/gatk/documentation/article.php?id=6925; Accessed 2018 Mar 15).

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference_genome.fa -V Str2.0_SpBa_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "my_indel_filter" -o Str2.0_SpBa_indels_filt1.vcf
```

* We then used the GATK SelectVariants tool with the "--excludeFiltered" option to produce a final VCF file of high-confidence indels for BQSR.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V Str2.0_SpBa_indels_filt1.vcf --excludeFiltered -o Str2.0_SpBa_indels_filt2.vcf
```

### BQSR

We performed base quality score recalibration (BQSR) and called variants for each sample using the [SpBa_StrOccCau2_recal.sh](SpBa_StrOccCau2_recal.sh) script. Below we provide further information on the relevant steps.

* Following an adapted version of the GATK guidelines for BQSR (https://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr; Accessed 2018 Mar 15), we first used the GATK BaseRecalibrator tool with default options other than supplying our high-confidence SNP and indel VCF files to the "-knownSites" option to produce a report file of covariation data for use in recalibrating the alignments.

```
java -Xmx4g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference_genome.fa -I sample_alignment.bam -knownSites Str2.0_SpBa_snps_filt2.vcf -knownSites Str2.0_SpBa_indels_filt2.vcf -o recal_data_table.tb 1>recal_data_table.log 2>recal_data_table.err
```

* We then used the GATK PrintReads tool supplying the covariation data report file to the "-BQSR" option to perform the actual recalibration of the sorted, duplicate-marked BAM files for all of our samples.

```
java -Xmx4g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T PrintReads -R reference_genome.fa -I sample_alignment.bam -BQSR recal_data_table.tb -o sample_alignment_recal.bam 1>sample_alignment_recal.log> 2>sample_alignment_recal.err
```

* We performed a new round of variant calling for each sample using the GATK HaplotypeCaller tool with the “--emitRefConfidence GVCF' option and suppling the recalibrated BAM file in order to produce sample-specific gVCF files.

```
java -Xmx4g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference_genome.fa -I sample_alignment_recal.bam --emitRefConfidence GVCF -o sample_alignment_recal.bam.gvcf 1>sample_alignment_recal.bam.gvcf.log 2>sample_alignment_recal.bam.gvcf.err
```

* We combined the individual sample gVCF files into a file of variants across samples by using the GATK GenotypeGVCFs tool with default settings to produce a new VCF file.

```
java -Xmx100g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference_genome.fa --variant sample{1..29}_recal.g.vcf -o Str2.0_SpBa_recal.vcf 1>Str2.0_SpBa_recal.vcf.log 2>Str2.0_SpBa_recal.vcf.err
```

### Variant filtering

We applied the filters up to our calculation of the mean site sequence coverage depth using the script [filt1_SpBa_StrOccCau2_recal.sh](filt1_SpBa_StrOccCau2_recal.sh).

* Similar to the pre-BQSR filtering, we filtered the new VCF file by first using the GATK SelectVariants tool with the "-selectType SNP" option to extract the SNPs.

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V SpBa_recal.vcf -selectType SNP -o Str2.0_SpBa_recal_snps.vcf 1>Str2.0_SpBa_recal_snps.vcf.log 2>Str2.0_SpBa_recal_snps.vcf.err
```

* We then filtered the SNPs using the GATK VariantFiltration tool with the options '--filterExpression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"' to filter out the variants where the QUAL field value (a measure of variant confidence) divided by the unfiltered depth of non-homogygous reference samples was less than 2.0; the Phred-scaled probability of strand bias was greater than 60.0; the strand odds ratio (SOR), which is an estimate of strand bias while accounting for the ratio of reads that cover each of the alleles at a site, was greater than 3.0; the root mean square mapping quality of all the reads at the site was less than 40.0; the u-based z-approximation from a rank-sum test comparing mapping qualities of reads in favor of the reference versus alternate alleles was less than -12.5; or the u-based z-approximation from a rank-sum test comparing the positions of the reference and alternate alleles within reads was less than -8.0 (https://software.broadinstitute.org/gatk/documentation/article.php?id=6925 ; Accessed 2018 Mar 15).

```
java -Xmx10g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference_genome.fa -V SpBa_recal.vcf --filterExpression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o Str2.0_SpBa_recal_snps_filt1.vcf 1>Str2.0_SpBa_recal_snps_filt1.log 2>Str2.0_SpBa_recal_snps_filt1.err
```

* We removed the filtered variants with the GATK SelectVariants tool.

```
java -Xmx10g -Djava.io.tmpdir=/temporary/directory/path -jar $GATK -T SelectVariants -R reference_genome.fa -V SpBa_recal_snps_filt1.vcf --excludeFiltered -o Str2.0_SpBa_recal_snps_filt2.vcf 1>Str2.0_SpBa_recal_snps_filt2.log 2>Str2.0_SpBa_recal_snps_filt2.err
```

* We removed any variants that fell within repetitive or low complexity regions using BEDTools version 2.25.0 [@quinlanBEDToolsFlexibleSuite2010] with options 'intersect -v -a <file.vcf> -b <masked_regions.bed> -header -wa' where 'file.vcf' was our filtered VCF file and 'masked_regions.bed' was the BED file of N-masked regions that we produced following our annotation of repetitive regions in the reference genome.

```
bedtools intersect -v -a Str2.0_SpBa_recal_snps_filt1.vcf -b masked_regions.bed -header -wa >Str2.0_SpBarecal_snps_filt2.vcf
```



* We next removed the mitochondrial genome variants from the VCF file using GNU Awk (GAWK) version 4.2.0 [@freesoftwarefoundationGNUAwk2017] with the expression 'awk '$1 !~ /^Sequoia_complete_mtGenome$/''.

```
awk '$1 !~ /^Sequoia_complete_mtGenome$/' Str2.0_SpBa_recal_snps.vcf >Str2.0_SpBa_recal_snps_filt1.vcf
```
* We used the GATK SelectVariants tool with the " --restrictAllelesTo BIALLELIC --excludeFiltered" options to remove the variants that failed the above filters from the VCF file and retain only biallelic sites.

```
java -Xmx4g -Djava.io.tmpdir=</temporary/directory/path> -jar GenomeAnalysisTK.jar -T SelectVariants -R reference_genome.fa -V SpBa_recal_snps_filt3.vcf --restrictAllelesTo BIALLELIC --excludeFiltered -o Str2.0_SpBa_recal_snps_filt4.vcf 1>Str2.0_SpBa_recal_snps_filt4.log 2>Str2.0_SpBa_recal_snps_filt4.err
```

* We used the dp_cov_script.sh tool from SPOW-BDOW-introgression-scripts version 1.1.1 [@hannaSPOWBDOWintrogressionscriptsVersion2017] to calculate the mean and standard deviation of the total unfiltered read depth across all samples per site.

```
./dp_cov_script.sh Str2.0_SpBa_recal_snps_filt4.vcf
```
output:

meanDP = 868.805,stdevDP = 679.386,number of sites = 17792804

* We then used the vcf_filter_highDP.sh script from genetics-tools version 1.0.1 [@hannaGeneticstoolsVersion2018a] to only retain sites in our VCF file with an unfiltered read depth less than 4,266X, which removed any sites exceeding the mean plus fives the standard deviation, as suggested by the GATK documentation (https://software.broadinstitute.org/gatk/documentation/article.php?id=3225 ; Accessed 2018 Mar 16).

```
vcf_filter_highDP.sh SpBa_recal_snps_filt4.vcf 4266 >Str2.0_SpBa_recal_snps_filt5.vcf
```

* We then compressed the VCF using the bgzip tool from HTSlib version 1.8 [@daviesHTSlib2018].

```
bgzip -c SpBa_recal_snps_filtfinal.vcf >SpBa_recal_snps_filtfinal.vcf.bgz
```

* We indexed the compressed VCF using the Tabix tool from HTSlib version 1.8 [@daviesHTSlib2018; @liTabixFastRetrieval2011].

```
tabix -p vcf SpBa_recal_snps_filtfinal.vcf.bgz
```

### Sequence coverage per individual

* We calculated mean and standard deviation of read depth across the final set of variant sites for each individual using the DP_sample_calc.sh script from genetics-tools version 1.0.1 [@hannaGeneticstoolsVersion2018a].

```
DP_sample_calc.sh SpBa_recal_snps_filtfinal.vcf | head -1 >SpBa_recal_snps_filtfinal.vcf_dp_means_stdev.txt
```





## References
