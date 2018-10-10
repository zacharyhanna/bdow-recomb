# Supplementary methodology for fine scale recombination rates in barred owls (*Strix varia*) from eastern North America

## Phasing

* We phased the genotypes of the eastern barred owl individuals using Beagle version 5.0 (03Jul18.40b) (Browning & Browning, 2007; Browning, 2018) with default options.

```
java -Xmx100g -Djava.io.tmpdir=/path/to/tmp -jar beagle.03Jul18.40b.jar gt=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.vcf out=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt 1>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.log 2>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.err
```

## Sequence coverage per individual

* We calculated mean and standard deviation of read depth across the final set of variant sites for each individual using the DP_sample_calc.sh script from genetics-tools version 1.0.1 (Hanna, 2018).

```
DP_sample_calc.sh SpBa_recal_snps_filtfinal.vcf | head -1 >SpBa_recal_snps_filtfinal.vcf_dp_means_stdev.txt
```





## References

Browning B. Beagle [Internet]. 2018. Available: https://faculty.washington.edu/browning/beagle/beagle.html

Browning SR, Browning BL. Rapid and Accurate Haplotype Phasing and Missing-Data Inference for Whole-Genome Association Studies By Use of Localized Haplotype Clustering. *The American Journal of Human Genetics*. 2007;81: 1084–1097. doi:10.1086/521987

Hanna ZR. genetics-tools. Version 1.0.1. *Zenodo*. 2018; doi:10.5281/zenodo.1257508
