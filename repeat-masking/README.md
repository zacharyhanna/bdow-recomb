# Repeat masking of reference genome

## Reference genome

We combined the following using cat (GNU core utilities) version 8.25 (Granlund & Stallman, 2017) in order to generate the reference genome:

1. StrOccCau_2.0_nuc.fa

2. mitochondrial genome (Hanna et al., 2017)

```
cat StrOccCau_2.0_nuc.fa mitochondrial_genome.fa >reference_genome.fa
```

## Repeat masking of genome assembly

We performed a homology-based annotation of repetitive regions using RepeatMasker version open-4.0.7 (Smit et al., 2013) with the repeat databases of the DFAM library version 2.0 (Hubley et al., 2016) and the Repbase RepeatMasker libraries version 20170127 (Bao et al., 2005; Jurka et al., 2005; Jurka, 1998; Jurka, 2000). Our installation of the RepeatMasker program relied on tandem repeats finder (TRF) version 4.09 (Benson, 1999) in addition to the NCBI BLAST+ version 2.7.1 (Camacho et al., 2009), RMBlast version 2.2.28 (Smit & Hubley, 2015), and HMMER version 3.1b2 (Eddy, 1998; http://hmmer.org) sequence search tools. We first performed a homology-based annotation of repetitive regions in the reference genome using RepeatMasker with the options '-gccalc -species aves'.

```
RepeatMasker -pa 12 -gccalc -species aves reference_genome.fa 1>reference_genome_RMask_aves.log 2>reference_genome_RMask_aves.err
```

* We created a de novo model of the repeats in the reference genome using RepeatModeler version 1.0.8 (Smit & Hubley, 2015). Our RepeatModeler installation used the RepeatMasker version open4.0.7 installation as detailed in the previous step with Repbase RepeatMasker libraries version 20170127, the RMBlast version 2.2.28 sequence search tool, TRF version 4.09, and two de novo repeat finding tools, RECON version 1.08 (Bao & Eddy, 2002) and RepeatScout version 1.0.5 (Price et al., 2005). We first used the RepeatMasker BuildDatabase tool to build a sequence database from the reference nuclear genome, StrOccCau_2.0_nuc.fa (we didn't want the model to include repeats from the mitochondrial genome).

```
BuildDatabase -name nuclear_genome nuclear_genome.fa 1>nuclear_genome_RModbuilddb.log 2>nuclear_genome_RModbuilddb.err
```

* We then ran RepeatModeler with default options.

```
RepeatModeler -pa 6 -database nuclear_genome 1>nuclear_genome_RMod.log 2>nuclear_genome_RMod.err
```

* We performed a final round of repeat masking by supplying the masked genome output from the homology-based masking to RepeatMasker with the options "-gccalc -lib <RepeatModeler_output>", where the "RepeatModeler_output" was the "consensi.fa.classified" repeat library file output from our RepeatModeler run.

```
RepeatMasker -pa 12 -gccalc -lib consensi.fa.classified reference_genome.fa.masked.fa 1>RMask.log 2>RMask.err
```

* Since RepeatMasker masked all of the repetitive regions with N characters, we used seqtk version 1.2-r94 (Li, 2016) with options "cutN -n1 -p100000000 -g" to create a browser extensible data (BED) formatted file of the coordinates of all of the N regions (these included the gap regions in the original reference genome).

```
seqtk cutN -n1 -p100000000 -g reference_genome.fa.masked.fa.masked >reference_genome.fa.masked.fa.masked_Nregions.bed
```

* We then used the GNU core utility sort version 8.25 (Haertel & Eggert, 2016) with options '-k1,1 -k2,2n' to sort the BED file.

```
sort -k1,1 -k2,2n reference_genome.fa.masked.fa.masked_Nregions.bed >reference_genome.fa.masked.fa.masked_Nregions_sorted.bed
```

## Final bed file of masked regions

The output of the above was our final bed file of masked regions.

```
mv reference_genome.fa.masked.fa.masked_Nregions_sorted.bed masked_regions.bed
```

## References

Bao Z, Eddy SR. Automated De Novo Identification of Repeat Sequence Families in Sequenced Genomes. *Genome Res*. 2002;12: 1269–1276. doi:10.1101/gr.88502

Bao W, Kojima KK, Kohany O. Repbase Update, a database of repetitive elements in eukaryotic genomes. *Mob DNA*. 2015;6: 1–6. doi:10.1186/s13100-015-0041-9

Benson G. Tandem repeats finder: a program to analyze DNA sequences. *Nucl Acids Res*. 1999;27: 573–580. doi:10.1093/nar/27.2.573

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10: 421. doi:10.1186/1471-2105-10-421

Eddy SR. Profile hidden Markov models. *Bioinformatic*. 1998;14: 755–763.

Granlund T, Stallman RM. cat (GNU coreutils) [Internet]. 2017. Available: http://www.gnu.org/software/coreutils/coreutils.html

Hanna ZR, Henderson JB, Sellas AB, Fuchs J, Bowie RCK, Dumbacher JP. Complete mitochondrial genome sequences of the northern spotted owl (*Strix occidentalis caurina*) and the barred owl (*Strix varia*; Aves: Strigiformes: Strigidae) confirm the presence of a duplicated control region. *PeerJ*. 2017;5: e3901. doi:10.7717/peerj.3901

Haertel M, Eggert P. sort (GNU coreutils) [Internet]. 2016. Available: http://www.gnu.org/software/coreutils/coreutils.html

Hubley R, Finn RD, Clements J, Eddy SR, Jones TA, Bao W, et al. The Dfam database of repetitive DNA families. *Nucleic Acids Res*. 2016;44: D81–D89. doi:10.1093/nar/gkv1272

Jurka J. Repeats in genomic DNA: mining and meaning. *Current Opinion in Structural Biology*. 1998;8: 333–337. doi:10.1016/S0959-440X(98)80067-5

Jurka J. Repbase Update: a database and an electronic journal of repetitive elements. *Trends in Genetics*. 2000;16: 418–420. doi:10.1016/S0168-9525(00)02093-X

Jurka J, Kapitonov VV, Pavlicek A, Klonowski P, Kohany O, Walichiewicz J. Repbase Update, a database of eukaryotic repetitive elements. *Cytogenetic and Genome Research*. 2005;110: 462–467. doi:10.1159/000084979

Li H. Seqtk [Internet]. 2016. Available: https://github.com/lh3/seqtk

Price AL, Jones NC, Pevzner PA. De novo identification of repeat families in large genomes. *Bioinformatics*. 2005;21: i351–i358. doi:10.1093/bioinformatics/bti1018

Smit A, Hubley R, Green P. RepeatMasker Open-4.0 [Internet]. 2013. Available: http://www.repeatmasker.org

Smit AFA, Hubley R. RepeatModeler Open-1.0 [Internet]. 2015. Available: http://www.repeatmasker.org

Smit A, Hubley R, National Center for Biotechnology Information. RMBlast [Internet]. 2015. Available: http://www.repeatmasker.org/RMBlast.html
