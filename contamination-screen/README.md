# StrOccCau2.0 assembly contamination screen

## Introduction
We followed the advice of the National Center for Biotechnology Information (NCBI) GenBank Submissions Staff to replicate the foreign contamination screen used by NCBI to screen genome assemblies for foreign contaminants (Jianli Dai, personal communication, 9 February 2018).  

## Common contaminants
1. Obtained NCBI database of common contaminants in eukaryotic sequences.  
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa
```
**Tools**  
GNU Wget version 1.17.1 (Rühsen et al. 2015)  

---
2. Decompressed file.  
```
gunzip contam_in_euks.fa.gz
```
**Tools**  
Gzip version 8.25 (Gailly et al. 2016)  

---
3. Made BLAST database.  
```
makeblastdb -in contam_in_euks.fa -parse_seqids -dbtype nucl
```
**Tools**  
NCBI BLAST+ version 2.7.1 (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST; Camacho et al. 2009)  

---
4. Performed blast search.
```
blastn -query assembly_0.1.fa -db contam_in_euks.fa -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' >contam_in_euks_blast.out
```
**Tools**  
NCBI BLAST+ version 2.7.1  
GNU Awk (GAWK) version 4.2.0 (Free Software Foundation 2017)  

---
5. Instructions from NCBI GenBank Submissions Staff (Jianli Dai, personal communication, 9 February 2018):  
> 'Process contaminant matches:  
> Contaminant matches from (1) are merged if they are from the same class of sequence (VECTOR, E.coli, IS, PHG) and they overlap or are separated by 50 bases or less.  
> If the total coverage of contaminant matches from (1) is >75% of the sequence length then flag the sequence as a contaminant to be excluded.  
> If the contaminant is classed as VECTOR, E.coli, IS:.*, PERM:.* or PHG:* and the contaminant location is within 100 bases of the the start or end of the sequence (or gap is the sequence is not contiguous), or within 100 bases of another contaminant match that is at an end, flag the contaminant span for trimming.  
> If the contaminant is one of the above, and the match is longer than 700 bases flag the contaminant span for trimming.  
> Other matches may be false alarms. Treat them as suspect spans and reBLAST the hit span plus 10 Kbp of flanking sequence on each side against nr, HTGS, related and unrelated chromosomes (as described below).'  

---
6. Filtered output.  
```
cat contam_in_euks_blast.out | awk '$0!="# 0 hits found"' >contam_in_euks_blast_filt2.out
```
**Tools**  
cat (GNU coreutils) version 8.25 (Granlund & Stallman 2017)  
GNU Awk (GAWK), 4.2.0  

---
7. Contig6318 and Contig5778 appeared to be contaminants, so we removed them.
```
bioawk -
c fastx '$name !~ /^Contig6318$/ && $name !~ /^Contig5778$/ {print ">"$name" "$comment"\n"$seq}' genome_assembly_0.1.fa >genome_assembly_0.2.fa
```
**Tools**  
bioawk version 1.0 (Li 2013)  


## Adaptor sequence screening
1. Obtained NCBI database of adaptor sequences in eukaryotic sequences.
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa
```
**Tools**  
GNU Wget version 1.17.1

---
2. Formatted database.
```
makeblastdb -in adaptors_for_screening_euks.fa -parse_seqids -dbtype nucl
```
**Tools**  
NCBI BLAST+ version 2.7.1

---
3. Obtained VecScreen.
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
```
**Tools**  
GNU Wget version 1.17.1

---
4. Made VecScreen executable.
```
chmod +x vecscreen
```
**Tools**  
chmod (GNU coreutils) version 8.25 (MacKenzie & Meyering 2017)

---
5. Ran VecScreen.
```
vecscreen -d adaptors_for_screening_euks.fa -f3 -i genome_assembly_0.1.fa -o euks_adapters_vecscreen.out
```
**Tools**  
VecScreen (National Center for Biotechnology Information 2012)

---
6. Obtained script to filter the VecScreen results.
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
```
**Tools**  
GNU Wget version 1.17.1

---
7. Made script executable.
```
chmod +x VSlistTo1HitPerLine.awk
```
**Tools**  
chmod (GNU coreutils) version 8.25

---
8. Filtered out the "Weak" and "Suspect Origin" hits.
```
./VSlistTo1HitPerLine.awk suspect=0 weak=0 euks_adapters_vecscreen.out >euks_adapters_vecscreen_filt.out
```

---
9. Filtered output.
```
cat euks_adapters_vecscreen_filt.out | awk '$1 !~ /^VecScreen_No_Hits/' >euks_adapters_vecscreen_filt2.out
```
**Tools**  
GAWK version 4.2.0  
cat (GNU coreutils) version 8.25  

---
10. We inspected all remaining adapter spans in the filtered output file and replaced those intervals with N's. If the adapter spans were near N gap regions already, then we replaced any intervening nucleotides with N's as well.
```
$ cat euks_adapters_vecscreen_filt2.out  

VecScreen_Strong Super-Scaffold_13 1887886 1887921  
VecScreen_Moderate Super-Scaffold_13 2825049 2825074  
VecScreen_Moderate Super-Scaffold_34 11359624 11359656  
VecScreen_Strong Super-Scaffold_50 3450380 3450416  
VecScreen_Moderate Super-Scaffold_51 5663227 5663251  
VecScreen_Moderate Super-Scaffold_1_obj 32542469 32542494  
VecScreen_Moderate Super-Scaffold_11_obj 33927095 33927121  
VecScreen_Strong Super-Scaffold_34_obj 2644603 2644646  
VecScreen_Strong Super-Scaffold_12_obj 5340268 5340298  
VecScreen_Strong Contig10033 5 47  
```
**Tools**  
cat (GNU coreutils) version 8.25

---
11. Removal 1.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Strong Super-Scaffold_13 1887886 1887921
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_13$/) print ">"$name" "$comment"\n"substr($seq,1,1887881)"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,1887922,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.1.fa
```
**Tools**  
bioawk version 1.0

---
12. Removal 2.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Moderate Super-Scaffold_13 2825049 2825074  
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_13$/) print ">"$name" "$comment"\n"substr($seq,1,2825048)"NNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,2825075,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.1.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.2.fa
```
**Tools**  
bioawk version 1.0  

---
13. Removal 3.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Moderate Super-Scaffold_34 11359624 11359656  
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_34$/) print ">"$name" "$comment"\n"substr($seq,1,11359594)"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,11359657,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.2.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.3.fa
```
**Tools**  
bioawk version 1.0  

---
14. Removal 4.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Strong Super-Scaffold_50 3450380 3450416  
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_50$/) print ">"$name" "$comment"\n"substr($seq,1,3450379)"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,3450417,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.3.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.4.fa
```
**Tools**  
bioawk version 1.0  

---   
15. Removal 5.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Moderate Super-Scaffold_51 5663227 5663251
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_51$/) print ">"$name" "$comment"\n"substr($seq,1,5663226)"NNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,5663252,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.4.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.5.fa
```
**Tools**  
bioawk version 1.0

---
16. Removal 6.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Moderate Super-Scaffold_1_obj 32542469 32542494
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_1_obj$/) print ">"$name" "$comment"\n"substr($seq,1,32542468)"NNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,32542495,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.5.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.6.fa
```
**Tools**  
bioawk version 1.0

---
17. Removal 7.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Moderate Super-Scaffold_11_obj 33927095 33927121
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_11_obj$/) print ">"$name" "$comment"\n"substr($seq,1,33927094)"NNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,33927122,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.6.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.7.fa
```
**Tools**  
bioawk version 1.0

---
18. Removal 8.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Strong Super-Scaffold_34_obj 2644603 2644646

```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_34_obj$/) print ">"$name" "$comment"\n"substr($seq,1,2644602)"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,2644647,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.7.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.8.fa
```
**Tools**  
bioawk version 1.0  

---
19. Removal 9.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Strong Super-Scaffold_12_obj 5340268 5340298  
```
bioawk -c fastx '{if ($name ~ /^Super-Scaffold_12_obj$/) print ">"$name" "$comment"\n"substr($seq,1,5340267)"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"substr($seq,5340299,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.8.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.9.fa
```
**Tools**  
bioawk version 1.0  

---
20. Removal 10.  
Used bioawk version 1.0 to replace the following adapter span with N's:  
VecScreen_Strong Contig10033 5 47  
```
bioawk -c fastx '{if ($name ~ /^Contig10033$/) print ">"$name" "$comment"\n"substr($seq,48,length($seq)); else print ">"$name" "$comment"\n"$seq}' /home/zhanna/owl/assemblies/StrOccCau_hybrid_0.2.9.fa >/home/zhanna/owl/assemblies/StrOccCau_hybrid_0.3.fa
```
**Tools**  
bioawk version 1.0  

---
## Ribosomal RNA screen  
Instructions from NCBI GenBank Submissions Staff (Jianli Dai, personal communication, 9 February 2018):  
>'Ribosomal RNA genes are the cause of many false positives because the include some segments that align to distantly related organisms. Segments that match rRNA genes are identified so that such segments are not reported as being foreign.  
>BLAST is used to screen the input sequences against a database of the rRNA gene sequences .'  

1. Obtained database.  
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz  
```
**Tools**  
GNU Wget version 1.17.1  

---
2. Decompressed database.  
```
gunzip rrna.gz  
```
**Tools**  
Gzip (GNU coreutils) version 8.25  

---
3. Renamed database.  
```
mv rrna rrna.fa  
```
**Tools**  
mv (GNU coreutils) version 8.25  

---
4. Made BLAST database.
```
makeblastdb -in rrna.fa -parse_seqids -dbtype nucl  
```
**Tools**  
NCBI BLAST+ version 2.7.1  

---
5. Ran BLAST command.  
The NCBI GenBank Submissions Staff recommended requiring hits to have >95% identity and be a minimum of 100 bases long (Jianli Dai, personal communication, 9 February 2018).  
```
blastn -query genome_assembly_0.1.fa -db rrna.fa -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 7 | awk '$4>=100' >genome_assembly_rrnablast_filt.out  
```
**Tools**  
NCBI BLAST+ version 2.7.1  

---
6. Filter output file.  
```
awk '$1 !~ /#/' genome_assembly_rrnablast_filt.out >genome_assembly_rrnablast_filt2.out
```
**Tools**  
GAWK version 4.2.0  
---
7. Sort output.  
```
cat genome_assembly_rrnablast_filt2.out | sort -k1,1 -k7 >genome_assembly_rrnablast_filt3.out
```
**Tools**
cat (GNU coreutils) version 8.25  
sort (GNU coreutils) version 8.25  

We used the resulting file as a reference for which regions were likely rrna.  

## Mitochondrial contamination
Instructions from NCBI GenBank Submissions Staff (Jianli Dai, personal communication, 9 February 2018):  
>'BLAST is used to screen the input sequences against a database of the mitochondrial genome sequences in the NCBI Reference Sequences (RefSeq) collection.'  

1. Obtained database.  
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz
```
**Tools**  
GNU Wget version 1.17.1  

---
2. Decompressed database.  
```
gunzip mito.nt.gz
```
**Tools**  
Gzip (GNU coreutils) version 8.25  

---
3. Renamed database.
```
mv mito.nt mito.nt.fa
```
**Tools**  
mv (GNU coreutils), 8.25  

---
4. Made BLAST database.  
```
makeblastdb -in mito.nt.fa -parse_seqids -dbtype nucl
```
**Tools**  
NCBI BLAST+ version 2.7.1

---
5. Ran BLASTN against database. We removed a '-perc_identity 98.6' originally filter suggested by NCBI in order to obtain more divergent hits.  
```
blastn -query genome_assembly_0.1.fa -db mito.nt.fa -out genome_assembly_mitoblast_filt.out -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -soft_masking true -outfmt 7  
```
**Tools**  
NCBI BLAST+ version 2.7.1  

---
6. Filtering instructions.  
Instructions from NCBI GenBank Submissions Staff (Jianli Dai, personal communication, 9 February 2018):  
>'If the total coverage of mitochondrial matches from (3) is >75% of the
sequence length then flag the sequence as being mitochindrial sequence to be
excluded.'

7. Filtered output.  
```
cat genome_assembly_mitoblast_filt.out | awk '$1 !~ /#/' >genome_assembly_mitoblast_filt2.out
```
**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

---
8. Sorted output.  
```
cat genome_assembly_mitoblast_filt2.out | sort -k1,1 -k7n >genome_assembly_mitoblast_filt2.txt
```
**Tools**  
cat (GNU coreutils) version 8.25  
sort (GNU coreutils) version 8.25  

---
9. Process output to combine overlapping blast result regions and output a file with the overlapping regions combined.  
```
python process_mito_matches.py genome_assembly_mitoblast_filt2.txt output_mito.txt
```
**Tools**
process_mito_matches.py (script included in this repository)  

---
10. Sort output and output unique scaffold names.  
```
cat output_mito.txt | awk '{print $1}' | sort -u >mito_poss_scafs.txt
```
**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  
sort (GNU coreutils) version 8.25  

---
11. Retrieved lengths of contigs with possible mitochondrial contamination.  
```
bioawk -c fastx '$name ~ /^<contig>$/ {print $name"\t"length($seq)}\' genome_assembly.fa
```
**Tools**  
bioawk version 1.0  

---
12. We calculated the percentage that the length of the span from the mitochondrial blast result represented of the total length of the scaffold / contig and, using the NCBI's guidelines, determined that the following contigs were possible mitochondrial sequences:  
Contig2141  
Contig3409  
Contig4223 - this also matched rrna  
Contig5159  
Contig5776  
Contig6157  
Contig7914  
Contig796  

---
13. We extracted these contigs.  
```
bioawk -c fastx '$name ~ /^Contig2141$/ || $name ~ /^Contig3409$/ || $name ~ /^Contig4223$/ || $name ~ /^Contig5159$/ || $name ~ /^Contig5776$/ || $name ~ /^Contig6157$/ || $name ~ /^Contig7914$/ || $name ~ /^Contig796$/ {print ">"$name" "$comment"\n"$seq}' genome_assembly_0.3.1.fa >genome_assembly_0.3.1_poss_mitos.fa  
```
**Tools**  
bioawk version 1.0  

---
14. Web BLAST  
We used the web version of the NCBI BLAST+ version 2.8.0 tool BLASTN (Altschul et al., 1997; Camacho et al., 2009) with default parameters to search the NCBI nucleotide collection (Johnson et al., 2008; Boratyn et al., 2013; Benson et al., 2015; NCBI Resource Coordinators, 2015) (NCBI-nt) for the top hits to these contigs.  

Contig2141 - matched mitochondrial DNA  
Contig3409 - matched bacteria  
Contig4223 - this matched pigeon rrna  
Contig5159 - matched bacteria  
Contig5776 - matched mitochondrial DNA  
Contig6157 - matched bacteria  
Contig7914 - matched bird genomic DNA  
Contig796 - matched bacteria  

---
15. Based on these BLAST results, we removed the following contigs from the assembly:  
Contig6157  
Contig5776  
Contig3409  
Contig2141  
Contig796  
Contig5159  
```
bioawk -c fastx '$name !~ /^Contig2141$/ && $name !~ /^Contig3409$/ && $name !~ /^Contig5159$/ && $name !~ /^Contig5776$/ && $name !~ /^Contig6157$/ && $name !~ /^Contig796$/ {print ">"$name" "$comment"\n"$seq}' genome_assembly_0.3.1.fa >genome_assembly_0.3.2.fa  
```
**Tools**  
bioawk version 1.0  

## Foreign chromosome screen  
We obtained FASTA sequence databases of eight taxonomic groups (we did not need the chordata group, as our genome fell in that group) by using query terms suggested by the NCBI GenBank Submissions Staff (Jianli Dai, personal communication, 9 February 2018). The query terms are included in the Perl scripts provided in this repository that we used to obtain the FASTA sequences. The Perl scripts are modifications of the Perl script model in the "Application 3: Retrieving large datasets" section of Sayers (2017).  

The taxonomic groups were the following: archaea, arthropoda, bacteria, fungi, other_eukaryota, other_metazoa, viridiplantae, and viruses_and_viroids.  

### Archaea  
```
cat 20180212_hybrid_archaea_blast_fm7.out | awk '$1 !~ /#/' >20180212_hybrid_archaea_blast_fm7_filt.out  
```
```
$ cat 20180212_hybrid_archaea_blast_fm7_filt.out  
```
> 0 results (empty file)  

**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Arthropoda  
```
cat 20180212_hybrid_arthropoda_blast_fm7.out | awk '$1 !~ /#/' >20180212_hybrid_arthropoda_blast_fm7_filt.out  
```
```
$ cat 20180212_hybrid_arthropoda_blast_fm7_filt.out  
```
> 0 results (empty file)  

**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Fungi  
```
cat 20180212_hybrid_fungi_blast_fm7.out | awk '$1 !~ /#/' >20180212_hybrid_fungi_blast_fm7_filt.out
```
```
$ cat 20180212_hybrid_fungi_blast_fm7_filt.out
```
> 0 results (empty file)  

**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Other eukaryota  
```
cat 20180212_hybrid_other_eukaryota_blast_fm7.out | awk '$1 !~ /#/' >20180212_hybrid_other_eukaryota_blast_fm7_filt.out
```
```
$ cat 20180212_hybrid_other_eukaryota_blast_fm7_filt.out  
Contig3553      NC_031180.2     98.596  285     4       0       2385    2669    212131  212415  4.54e-146       523  
```
> 1 result  

**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Other metazoa  
```
cat genome_assembly_other_metazoa_blast_fm7.out | awk '$1 !~ /#/' >genome_assembly_other_metazoa_blast_fm7_filt.out  
```
```
$ cat 20180212_hybrid_other_metazoa_blast_fm7_filt.out  
Contig8503      NC_031499.1     98.352  182     3       0       1578    1759    3014058 3014239 8.06e-87        329  
Contig8503      NC_031499.1     99.091  110     1       0       1658    1767    3052823 3052932 1.34e-50        208  
Contig4969      NC_031499.1     98.352  182     3       0       475     656     3014058 3014239 3.37e-87        329  
Contig4969      NC_031499.1     99.091  110     1       0       555     664     3052823 3052932 5.59e-51        208  
Contig6727      NC_031499.1     98.352  182     3       0       3800    3981    3014058 3014239 5.18e-87        329  
Contig6727      NC_031499.1     99.091  110     1       0       3880    3989    3052823 3052932 8.61e-51        208  
```
**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Viridiplantae  
```
cat genome_assembly_viridiplantae_blast_fm7.out | awk '$1 !~ /#/' >genome_assembly_viridiplantae_blast_fm7_filt.out  
```
```
cat genome_assembly_viridiplantae_blast_fm7_filt.out  
```
> Produces a large output (many hits).  

**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Viruses
```
cat genome_assembly_viruses_blast_fm7.out | awk '$1 !~ /#/' >genome_assembly_viruses_blast_fm7_filt.out  
```
```
$ cat 20180212_hybrid_viruses_blast_fm7_filt.out  
Super-Scaffold_28       NC_001407.1     100.000 101     0       0       11573931        11574031        7570    7470    3.70e-46        202  
Super-Scaffold_28       NC_008094.1     100.000 101     0       0       11573931        11574031        1148    1048    3.70e-46        202  
Contig6318      NC_001422.1     100.000 1982    0       0       489     2470    5386    3405    0.0     3975  
Contig6318      NC_001422.1     100.000 488     0       0       1       488     488     1       0.0     979  
Contig5778      NC_001422.1     99.828  2914    5       0       1       2914    490     3403    0.0     5783  
```
**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  

### Bacteria  
We saved processing this for last as it had the greatest number of results that we needed to investigate.  
* Perform BLAST to database.  

---
* Filter output.  
```
cat genome_assembly_bacteria_blast_fm7.out | awk '$1 !~ /#/ {print $1,$7,$8}' |sort -k1,1 -k2n|uniq >genome_assembly_bacteria_blast_fm7_filt1.out  
```
**Tools**  
cat (GNU coreutils) version 8.25  
GAWK version 4.2.0  
sort (GNU coreutils) version 8.25  
uniq (GNU coreutils) version 8.25 (Stallman & MacKenzie 2017)  

---
* Combine overlapping records.  
```
python
process_bacteria_matches.py genome_assembly_bacteria_blast_fm7_filt1.out genome_assembly_bacteria_blast_fm7_filt2.out
```
```
$ cat 20180212_hybrid_bacteria_blast_fm7_filt2.out  
Contig805	1	314  
Contig805	415	739  
Contig805	840	1153  
Contig2449	373	1416  
Contig796	4	335  
Contig796	788	1119  
Contig2319	1	699  
Contig2319	800	1181  
Contig3553	1	2870  
Contig3409	1	3066  
Contig6157	1	1895  
Contig1521	1	968  
Contig1521	1069	1523  
Contig5778	1	2914  
Contig6318	1	1516  
Contig6318	1517	2470  
Contig5159	1	1231  
Contig5159	1281	1815  
Contig3981	1	1511  
Contig1909	1	481  
Contig1909	582	1062  
```
Contigs already gone from mitochondrial screen:  
Contig6157  
Contig5776  
Contig3409  
Contig2141  
Contig796  
Contig5159  

Contigs already gone from common contaminant screen:  
Contig6318 and Contig5778  

Removed from other_eukaryota blast result:  
Contig3553  

New, reduced list:  
Contig805	1	314  
Contig805	415	739  
Contig805	840	1153  
Contig2449	373	1416  
Contig2319	1	699  
Contig2319	800	1181  
Contig1521	1	968  
Contig1521	1069	1523  
Contig3981	1	1511  
Contig1909	1	481  
Contig1909	582	1062  

**Tools**  
process_bacteria_matches.py (provided in this repository)  
cat (GNU coreutils) version 8.25  

---
* We extracted the sequences of the suspect contigs.  
```
bioawk -c fastx '$name ~ /^Contig805$/ || $name ~ /^Contig2449$/ || $name ~ /^Contig2319$/ || $name ~ /^Contig3553$/ || $name ~ /^Contig1521$/ || $name ~ /^Contig3981$/ || $name ~ /^Contig1909$/ {print ">"$name" "$comment"\n"$seq}' genome_assembly_0.3.2.fa >genome_assembly_0.3.2_bacteria_suspects.fa
```
**Tools**  
bioawk version 1.0  

---
* We used the web version of the NCBI BLAST+ version 2.8.0 tool BLASTN (Altschul et al., 1997; Camacho et al., 2009) with default parameters to search the NCBI nucleotide collection (Johnson et al., 2008; Boratyn et al., 2013; Benson et al., 2015; NCBI Resource Coordinators, 2015) (NCBI-nt) for the top hits to these contigs.

They were all bacteria.  

---
* We removed those contigs by running the below:

```
bioawk -c fastx '$name !~ /^Contig805$/ && $name !~ /^Contig2449$/ && $name !~ /^Contig2319$/ && $name !~ /^Contig3553$/ && $name !~ /^Contig1521$/ && $name !~ /^Contig3981$/ && $name !~ /^Contig1909$/ {print ">"$name" "$comment"\n"$seq}' genome_assembly_0.3.2.fa >genome_assembly_0.3.3.fa
```

**Tools**  
bioawk version 1.0  

### Final assembly version

genome_assembly_0.3.3.fa became our cleaned assembly, which we renamed "StrOccCau_2.0_nuc.fa".

```
mv genome_assembly_0.3.3.fa StrOccCau_2.0_nuc.fa
```

## References

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, et al. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucl Acids Res. 1997;25: 3389–3402. doi:[10.1093/nar/25.17.3389](https://doi.org/10.1093/nar/25.17.3389)  

Benson DA, Cavanaugh M, Clark K, Karsch-Mizrachi I, Ostell J, Pruitt KD, et al. GenBank. Nucleic Acids Res. 2018;46: D41–D47. doi:[10.1093/nar/gkx1094](https://doi.org/10.1093/nar/gkx1094)  

Boratyn GM, Camacho C, Cooper PS, Coulouris G, Fong A, Ma N, et al. BLAST: a more efficient report with usability improvements. Nucl Acids Res. 2013;41: W29–W33. doi:[10.1093/nar/gkt282](https://doi.org/10.1093/nar/gkt282)  

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009;10: 421. doi:10.1186/1471-2105-10-421  

Free Software Foundation. GNU Awk [Internet]. 2017. Available: https://www.gnu.org/software/gawk  

Gailly J, Adler M, Meyering J, Eggert P. Gzip (GNU coreutils) [Internet]. 2016. Available: http://www.gnu.org/software/coreutils/coreutils.html  

Granlund T, Stallman RM. cat (GNU coreutils) [Internet]. 2017. Available: https://www.gnu.org/software/coreutils/coreutils.html  

Haertel M, Eggert P. sort (GNU coreutils) [Internet]. 2016. Available: https://www.gnu.org/software/coreutils/coreutils.html  

Johnson M, Zaretskaya I, Raytselis Y, Merezhuk Y, McGinnis S, Madden TL. NCBI BLAST: a better web interface. Nucl Acids Res. 2008;36: W5–W9. doi:[10.1093/nar/gkn201](https://doi.org/10.1093/nar/gkn201)  

Li H. bioawk [Internet]. 2013. Available: https://github.com/lh3/bioawk  

Parker M, MacKenzie D, Meyering J. mv (GNU coreutils) [Internet]. 2017. Available: https://www.gnu.org/software/coreutils/coreutils.html  

MacKenzie D, Meyering J. chmod (GNU coreutils) [Internet]. 2017. Available: https://www.gnu.org/software/coreutils/coreutils.html  

National Center for Biotechnology Information. VecScreen [Internet]. 2012. Available: https://www.ncbi.nlm.nih.gov/tools/vecscreen/  

Rühsen T, Shah D, Scrivano G, Nikšić H. GNU Wget [Internet]. 2015. Available: https://www.gnu.org/software/wget/  

Sayers E. Sample Applications of the E-utilities. Entrez Programming Utilities Help [Internet]. Bethesda, Maryland, United States of America: National Center for Biotechnology Information, United States of America; 2017. Available: https://www.ncbi.nlm.nih.gov/books/NBK25498  

Stallman RM, MacKenzie D. uniq (GNU coreutils) [Internet]. 2017. Available: http://www.gnu.org/software/coreutils/coreutils.html  
