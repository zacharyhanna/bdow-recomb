#!/bin/bash
reference="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa"
filt2_snps="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_snps_filt2.vcf"
filt2_indels="/media/walllab/zhanna/owl/20180628_SpBa_StrOccCau2_indels_filt2.vcf"
mybam=$1
recal_data_table="${mybam}_SpBarecaldata.tb"
recal_data_tablelog="${mybam}_SpBarecaldata.log"
recal_data_tableerr="${mybam}_SpBarecaldata.err"
set -x
java -Xmx4g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T BaseRecalibrator -R $reference -I $mybam -knownSites $filt2_snps -knownSites $filt2_indels -o $recal_data_table 1>$recal_data_tablelog 2>$recal_data_tableerr
set +x
post_recal_table="${mybam}_SpBapostrecaldata.tb"
post_recal_tablelog="${mybam}_SpBapostrecaldata.log"
post_recal_tableerr="${mybam}_SpBapostrecaldata.err"
set -x
java -Xmx4g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T BaseRecalibrator -R $reference -I $mybam -knownSites $filt2_snps -knownSites $filt2_indels -BQSR $recal_data_table -o $post_recal_table 1>$post_recal_tablelog 2>$post_recal_tableerr
set +x
plots="${mybam}_SpBarecal_plots.pdf"
plotslog="${mybam}_SpBarecal_plots.log"
plotserr="${mybam}_SpBarecal_plots.err"
set -x
java -Xmx4g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T AnalyzeCovariates -R $reference -before $recal_data_table -after $post_recal_table -plots $plots 1>$plotslog 2>$plotserr
set +x
recal_reads="${mybam}_SpBarecal.bam"
recal_readslog="${mybam}_SpBarecal.log"
recal_readserr="${mybam}_SpBarecal.err"
set -x
java -Xmx4g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T PrintReads -R $reference -I $mybam -BQSR $recal_data_table -o $recal_reads 1>$recal_readslog 2>$recal_readserr
set +x
out_gvcf="${mybam}_SpBarecal.g.vcf"
out_log="${mybam}_SpBarecal.g.vcf.log"
out_err="${mybam}_SpBarecal.g.vcf.err"
set -x
java -Xmx4g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T HaplotypeCaller -R $reference -I $recal_reads --emitRefConfidence GVCF -o $out_gvcf 1>$out_log 2>$out_err
set +x
