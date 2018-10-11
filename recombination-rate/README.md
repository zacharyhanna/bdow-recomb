# Supplementary methodology for fine scale recombination rates in barred owls (*Strix varia*) from eastern North America

## Phasing

* We phased the genotypes of the eastern barred owl individuals using Beagle version 5.0 (03Jul18.40b) (Browning & Browning, 2007; Browning, 2018) with default options.

```
java -Xmx100g -Djava.io.tmpdir=/path/to/tmp -jar beagle.03Jul18.40b.jar gt=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.vcf out=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt 1>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.log 2>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.err
```

## References

Browning B. Beagle [Internet]. 2018. Available: https://faculty.washington.edu/browning/beagle/beagle.html

Browning SR, Browning BL. Rapid and Accurate Haplotype Phasing and Missing-Data Inference for Whole-Genome Association Studies By Use of Localized Haplotype Clustering. *The American Journal of Human Genetics*. 2007;81: 1084–1097. doi:10.1086/521987
