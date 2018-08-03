#!/bin/bash
while read scaf; do
  set -x
  ~/bin/LDhat/convert -seq ${scaf}_BADOeastGrEq1Mb.ldhat.sites -loc ${scaf}_BADOeastGrEq1Mb.ldhat.locs -prefix ${scaf}_BADOeastGrEq1Mb_cvrt
  vcftools --gzvcf /media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt5_BADOeastGrEq1Mb.gt.vcf.gz \
  --out /media/walllab/zhanna/owl/ldhat/${scaf}_BADOeastGrEq1Mb \
  --ldhat --phased --chr $scaf
  set +x
done <$1
