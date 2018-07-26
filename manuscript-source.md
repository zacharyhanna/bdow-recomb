---
title: Fine scale recombination rates in barred owls (*Strix varia*) from eastern North America
author: Zachary R. Hanna
author: Jeffrey D. Wall
affiliation: Institute for Human Genetics, University of California, San Francisco, San Francisco, California, United States of America
---

## Methods

### Sequence data

We utilized high-coverage genomic sequence data from twelve eastern and thirteen western *Strix varia* [XXX] as well as from four spotted owls (*S. occidentalis*) [@hannaWholegenomeSequencesSuggest2018; @hannaNorthernSpottedOwl2017; XXX]. We downloaded the raw genomic sequences from the NCBI Sequence Read Archive (SRA) accessions SRR4011595-SRR4011597, SRR4011614-SRR4011617, SRR6048826-SRR6048853 (TableSampInfo). We used Trimmomatic version 0.36 [@bolgerTrimmomaticFlexibleTrimmer2014] to remove adapter sequences as well as low quality leading and trailing bases. We used BWA-MEM version 0.7.17-r1188 [@liAligningSequenceReads2013] to align the processed sequences to our new reference *Strix occidentalis caurina* nuclear genome assembly together with the Hanna et al. [-@hannaCompleteMitochondrialGenome2017] *S. o. caurina* mitochondrial genome sequence. We merged the alignments, sorted the alignments, and marked duplicate sequences with Picard version 2.17.6 (http://broadinstitute.github.io/picard).

### Variant calling and filtering

We called variants for each sample using the Genome Analysis Toolkit (GATK) version 3.8.0 HaplotypeCaller tool and then combined the variants across samples using the GATK GenotypeGVCFs tool [@depristoFrameworkVariationDiscovery2011; @mckennaGenomeAnalysisToolkit2010; @vanderauweraFastQDataHigh2013; @poplinScalingAccurateGenetic2017]. Using the GATK VariantFiltration and SelectVariants tools, we implemented a hard filtering of the variants to generate a set of high quality single nucleotide polymorphisms (SNPs) and indels by following an adapted version of the GATK guidelines for hard filtering of variants. We utilized these filtered sets SNPs and indels for base quality score recalibration (BQSR). We first recalibrated the alignment files for each sample using the GATK BaseRecalibrator and PrintReads tools and then performed a new round of individual sample variant calling with the GATK HaplotypeCaller tool followed by combining variants across samples with the GATK GenotypeGVCFs tool. We extracted SNPs from the resulting variant call format (VCF) file using the GATK SelectVariants tool.

We used GNU AWK (GAWK) version 4.2.0 [@freesoftwarefoundationGNUAwk2017] to exclude SNPs that fell on the mitochondrial genome and then removed SNPs in repetitive or low complexity regions using BEDTools version 2.25.0 [@quinlanBEDToolsFlexibleSuite2010] along with the file of repetitive regions that we generated following repeat annotation of the new nuclear genome. We excluded low quality SNPs using the GATK VariantFiltration tool and extracted only biallelic SNPs by using the GATK SelectVariants tool to omit multiallelic sites. We calculated the mean and standard deviation of the site-level coverage of the remaining SNPs and then removed sites with excessive coverage with the vcf_filter_highDP.sh script from genetics-tools version 1.0.1 [@hannaGeneticstoolsVersion2018a]. We determined the mean and standard deviation of the sequence coverage depth for each individual across the final set of SNPs using the DP_sample_calc.sh script from genetics-tools version 1.0.1 [@hannaGeneticstoolsVersion2018a]. For use in downstream analyses, we compressed the final set of SNPs using the bzip tool from HTSlib version 1.8 [@daviesHTSlib2018] and then indexed the file using the Tabix tool from HTSlib version 1.8 [@daviesHTSlib2018; @liTabixFastRetrieval2011].

## Results



## Discussion

## References
