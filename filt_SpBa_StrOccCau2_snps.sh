#!/bin/bash
reference="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa"
raw_vcf="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2.vcf"
snps_out="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps.vcf"
select_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps.err"
select_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps.log"
filt1_snps="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt1.vcf"
filt1_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt1.err"
filt1_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt1.log"
filt2_snps="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt2.vcf"
filt2_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt2.err"
filt2_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt2.log"
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $raw_vcf -selectType SNP -o $snps_out 1>$select_log 2>$select_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T VariantFiltration -R $reference -V $snps_out --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName "my_snp_filter" -o $filt1_snps 1>$filt1_log 2>$filt1_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $filt1_snps --excludeFiltered -o $filt2_snps 1>$filt2_log 2>$filt2_err
