#!/bin/bash
cat /media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal.vcf | \
awk '$1 ~ /^##contig=<ID=/ {gsub("^##contig=<ID=",""); gsub(",length=","\t"); gsub(">$",""); print $0}' | \
sort -k2,2n | \
awk '$2 >= 1000000 {print $1}' \
>/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_contigs_wLeng_sortGrEq1Mb_scafs.intervals
