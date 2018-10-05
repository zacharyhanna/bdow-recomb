---
title: Supplementary methodology for fine scale recombination rates in barred owls (*Strix varia*) from eastern North America
author: Zachary R. Hanna
author: Jeffrey D. Wall
affiliation: Institute for Human Genetics, University of California, San Francisco, San Francisco, California, United States of America
---







### Phasing

* We phased the genotypes of the eastern barred owl individuals using Beagle version 5.0 (03Jul18.40b) [@browningRapidAccurateHaplotype2007; @browningBeagle2018] with default options.

```
java -Xmx100g -Djava.io.tmpdir=/path/to/tmp -jar beagle.03Jul18.40b.jar gt=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.vcf out=Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt 1>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.log 2>Str2.0_SpBa_recal_snps_EBdowGrEq1Mb.gt.err
```



### Sequence coverage per individual

* We calculated mean and standard deviation of read depth across the final set of variant sites for each individual using the DP_sample_calc.sh script from genetics-tools version 1.0.1 [@hannaGeneticstoolsVersion2018a].

```
DP_sample_calc.sh SpBa_recal_snps_filtfinal.vcf | head -1 >SpBa_recal_snps_filtfinal.vcf_dp_means_stdev.txt
```





## References
