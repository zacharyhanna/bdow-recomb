#!/bin/bash
reference="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa"
raw_vcf="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal.vcf"
snps_out="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps.vcf"
select_err="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps.err"
select_log="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps.log"
filt1_snps="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt1.vcf"
filt1_err="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt1.err"
filt1_log="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt1.log"
filt2_snps="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt2.vcf"
filt2_err="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt2.err"
filt2_log="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt2.log"
masked_regions="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa.masked.fa.masked_Nregions_sorted.bed"
filt3_snps="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt3.vcf"
filt4_snps="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt4.vcf"
filt4_err="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt4.err"
filt4_log="/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt4.log"
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $raw_vcf -selectType SNP -o $snps_out 1>$select_log 2>$select_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T VariantFiltration -R $reference -V $snps_out --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName "my_snp_filter" -o $filt1_snps 1>$filt1_log 2>$filt1_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $filt1_snps --excludeFiltered -o $filt2_snps 1>$filt2_log 2>$filt2_err
bedtools intersect -v -a $filt2_snps -b $masked_regions -header -wa >$filt3_snps
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $filt3_snps --restrictAllelesTo BIALLELIC -XL Sequoia_complete_mtGenome --excludeFiltered -o $filt4_snps 1>$filt4_log 2>$filt4_err
