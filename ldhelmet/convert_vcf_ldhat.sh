#!/bin/bash
while read scaf; do
  set -x
  vcftools --gzvcf /media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt5_BADOeastGrEq1Mb.gt.vcf.gz \
  --out /media/walllab/zhanna/owl/ldhat/${scaf}_BADOeastGrEq1Mb \
  --ldhat --phased --chr $scaf
  set +x
done <$1
