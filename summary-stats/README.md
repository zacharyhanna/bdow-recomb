# Supplementary methodology for summary statistics

## Sequence coverage per individual

* We calculated mean and standard deviation of read depth across the final set of variant sites for each individual using the DP_sample_calc.sh script from genetics-tools version 1.0.1 (Hanna, 2018).

```
DP_sample_calc.sh SpBa_recal_snps_filtfinal.vcf | head -1 >SpBa_recal_snps_filtfinal.vcf_dp_means_stdev.txt
```

## References

Hanna ZR. genetics-tools. Version 1.0.1. *Zenodo*. 2018; doi:10.5281/zenodo.1257508
