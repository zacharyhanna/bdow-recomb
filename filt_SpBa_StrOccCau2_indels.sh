#!/bin/bash
reference="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa"
raw_vcf="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2.vcf"
indels_out="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels.vcf"
select_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels.err"
select_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels.log"
filt1_indels="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt1.vcf"
filt1_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt1.err"
filt1_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt1.log"
filt2_indels="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt2.vcf"
filt2_err="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt2.err"
filt2_log="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt2.log"
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $raw_vcf -selectType INDEL -o $indels_out 1>$select_log 2>$select_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T VariantFiltration -R $reference -V $indels_out --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "my_indel_filter" -o $filt1_indels 1>$filt1_log 2>$filt1_err
java -Xmx10g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T SelectVariants -R $reference -V $filt1_indels --excludeFiltered -o $filt2_indels 1>$filt2_log 2>$filt2_err
