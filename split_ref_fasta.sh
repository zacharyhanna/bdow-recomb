#!/bin/bash
reference="/media/walllab/zhanna/owl/assemblies/StrOccCau_2.0_nuc_finalMito.fa"
while read scaf; do
  set -x
  bioawk -v myscaf=$scaf -c fastx '$name == myscaf {print ">"$name"\n"$seq}' $reference >${scaf}_StrOccCau_2.0.fa
  set +x
done <$1
