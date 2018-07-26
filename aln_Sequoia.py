import sys
import mimetypes
import gzip
import subprocess
import time

def Run(command): #run the command in the bash shell
    #print "not actually running"
    try:
        proc = subprocess.check_output(command, shell=True)
        return proc
    except subprocess.CalledProcessError:
        pass
    except:
        pass
'''
def Run(command):
    pass
'''
run_dict = {
"SRR4011615":"Sequoia_550bpnoPCR_2014Sep10_CAS",
"SRR4011616":"Sequoia_900bpPCR_2014Sep10_CAS",
"SRR4011596":"Sequoia_2011Aug_Ill",
"SRR4011595":"Sequoia_2011Jul_UCSF",
"SRR4011597":"Sequoia_2013Oct07_UCB",
"SRR4011614":"Sequoia_2014May1224_CAS",
"SRR4011617":"Sequoia_Hydro6_2014Sep10_CAS"
}

samp1 = "Sequoia"
sampfull1 = "SRR4011595"
runID1 = "2011Jul_UCSF_ln1"
lib1 = "2011NexFragSequoia"
lane1 = "lane1"
insert1 = "247"
maxinsert1 = "719"
pair1 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011595_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011595_2_2017Feb07AdptTrim_prd.fastq"
unpair1 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011595_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011595_2_2017Feb07AdptTrim_unprd.fastq"

samp2 = "Sequoia"
sampfull2 = "SRR4011596"
runID2 = "2011Aug_Ill_ln1"
lib2 = "2011NexFragSequoia"
lane2 = "lane1"
insert2 = "247"
maxinsert2 = "719"
pair2 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011596_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011596_2_2017Feb07AdptTrim_prd.fastq"
unpair2 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011596_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011596_2_2017Feb07AdptTrim_unprd.fastq"

samp3 = "Sequoia"
sampfull3 = "SRR4011597"
runID3 = "2013Oct07_UCB_ln1"
lib3 = "2013Nex809Sequoia"
lane3 = "lane1"
insert3 = "566"
maxinsert3 = "1342"
pair3 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011597_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011597_2_2017Feb07AdptTrim_prd.fastq"
unpair3 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011597_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011597_2_2017Feb07AdptTrim_unprd.fastq"

samp7 = "Sequoia"
sampfull7 = "SRR4011615"
runID7 = "2014Sep10_ln1"
lib7 = "2014550bpnoPCR"
lane7 = "lane1"
insert7 = "619"
maxinsert7 = "1147"
pair7 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011615_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011615_2_2017Feb07AdptTrim_prd.fastq"
unpair7 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011615_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011615_2_2017Feb07AdptTrim_unprd.fastq"

samp8 = "Sequoia"
sampfull8 = "SRR4011616"
runID8 = "2014Sep10_ln1"
lib8 = "20149000bpPCR"
lane8 = "lane1"
insert8 = "687"
maxinsert8 = "919"
pair8 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011616_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011616_2_2017Feb07AdptTrim_prd.fastq"
unpair8 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011616_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011616_2_2017Feb07AdptTrim_unprd.fastq"

samp9 = "Sequoia"
sampfull9 = "SRR4011614"
runID9 = "2014May1224_ln1"
lib9 = "2014Nex560Sequoia"
lane9 = "lane1"
insert9 = "560"
maxinsert9 = "660"
pair9 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011614_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011614_2_2017Feb07AdptTrim_prd.fastq"
unpair9 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011614_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011614_2_2017Feb07AdptTrim_unprd.fastq"

samp11 = "Sequoia"
sampfull11 = "SRR4011617"
runID11 = "2014Sep10_ln1"
lib11 = "2011Hydro6Sequoia"
lane11 = "lane1"
insert11 = "500"
maxinsert11 = "708"
pair11 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011617_1_2017Feb07AdptTrim_prd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011617_2_2017Feb07AdptTrim_prd.fastq"
unpair11 = "/media/walllab/zhanna/owl/Sequoia_filt/SRR4011617_1_2017Feb07AdptTrim_unprd.fastq /media/walllab/zhanna/owl/Sequoia_filt/SRR4011617_2_2017Feb07AdptTrim_unprd.fastq"

def merge_sams(sam1, sam2):
    out_parse = sam1.split("_")
    out_ls = out_parse[:-1]
    outfile_bs = ("_").join(out_ls)
    outfile = outfile_bs + "_mates_singles.sam"
    out_log = outfile_bs + "_mates_singles.log"
    out_err = outfile_bs + "_mates_singles.err"
    merge_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $PICARD MergeSamFiles MAX_RECORDS_IN_RAM=200000 INPUT=" + sam1 + " INPUT=" + sam2 + " OUTPUT=" + outfile + " SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=/media/walllab/zhanna/tmp 1>" + out_log + " 2>" + out_err
    return outfile, merge_comm

def sort_sam(in_sam):
    out_parse = in_sam[:-4] #take off the .sam
    out_bam = out_parse + "_sorted.bam"
    out_log = out_parse + "_sorted.log"
    out_err = out_parse + "_sorted.err"
    sort_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $PICARD SortSam MAX_RECORDS_IN_RAM=200000 INPUT=" + in_sam + " OUTPUT=" + out_bam + " SORT_ORDER=coordinate TMP_DIR=/media/walllab/zhanna/tmp 1>" + out_log + " 2>" + out_err
    return out_bam, sort_comm

def mark_dups(in_bam):
    out_parse = in_bam[:-4] #take off the .bam
    out_bam = out_parse + "_dedup.bam"
    out_log = out_parse + "_dedup.log"
    out_err = out_parse + "_dedup.err"
    out_metrics = out_parse + "_dedup_metrics.txt"
    mark_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $PICARD MarkDuplicates MAX_RECORDS_IN_RAM=200000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=/media/walllab/zhanna/tmp INPUT=" + in_bam + " OUTPUT=" + out_bam + " METRICS_FILE=" + out_metrics + " 1>" + out_log + " 2>" + out_err
    return out_bam, mark_comm

# copied make_aln_comm so that I can just have one that returns the file names and doesn't run anything
def make_aln_comm2(samp, fullsamp, runID, library, lane, insert, maxinsert, paired_ins, unpaired_ins): #receives a list of parameters as input
    global date
    global ref_genome
    aln_file = "/media/walllab/zhanna/owl/Sequoia_filt/" + fullsamp + "_" + date + "alnStrOccCau1_mates.sam"
    aln_log = "/media/walllab/zhanna/owl/Sequoia_filt/" + fullsamp + "_" + date + "alnStrOccCau1_mates.log"
    aln_comm = str("bwa mem -M -t 42 -R '@RG\\tID:" + runID + "\\tSM:" + samp + "\\tPL:illumina\\tLB:" + library + "\\tPU:" + lane + "\\tPI:" + insert + "' -w " + maxinsert + " " + ref_genome + " " + paired_ins + " >" + aln_file + " 2>" + aln_log)
    print aln_comm
    Run(aln_comm)
    aln_file_sing = "/media/walllab/zhanna/owl/Sequoia_filt/" + fullsamp + "_" + date + "alnStrOccCau1_singles.sam"
    aln_log_sing = "/media/walllab/zhanna/owl/Sequoia_filt/" + fullsamp + "_" + date + "alnStrOccCau1_singles.log"
    aln_comm_sing = str("zcat " + unpaired_ins + " | bwa mem -M -t 42 -R '@RG\\tID:" + runID + "\\tSM:" + samp + "\\tPL:illumina\\tLB:" + library + "\\tPU:" + lane + "\\tPI:" + insert + "' " + ref_genome + " - >" + aln_file_sing + " 2>" + aln_log_sing)
    print aln_comm_sing
    Run(aln_comm_sing)
    merged_two = merge_sams(aln_file, aln_file_sing)
    print merged_two[1]
    Run(merged_two[1])
    sorted = sort_sam(merged_two[0])
    print sorted[1]
    Run(sorted[1])
    dedup = mark_dups(sorted[0])
    print dedup[1]
    Run(dedup[1])
    #wgs_met = run_wgs_metrics(dedup[0])
    #print wgs_met
    #Run(wgs_met)
    return dedup[0]

def merge_ls_bams(ls_bams):
    global date
    outfile = "/media/walllab/zhanna/owl/Sequoia_filt/Sequoia_" + date + "alnStrOccCau1_merged.bam"
    out_log = "/media/walllab/zhanna/owl/Sequoia_filt/Sequoia_" + date + "alnStrOccCau1_merged.log"
    out_err = "/media/walllab/zhanna/owl/Sequoia_filt/Sequoia_" + date + "alnStrOccCau1_merged.err"
    ls_comm_ins = []
    for bam in ls_bams:
        in_format = "INPUT=" + bam
        ls_comm_ins.append(in_format)
    all_ins = ' '.join(ls_comm_ins)
    merge_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $PICARD MergeSamFiles MAX_RECORDS_IN_RAM=400000 " + all_ins + " OUTPUT=" + outfile + " SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=/media/walllab/zhanna/tmp 1>" + out_log + " 2>" + out_err
    return outfile, merge_comm

def run_wgs_metrics(in_bam):
    out_parse = in_bam[:-4] #take off the .bam
    out_metrics = out_parse + "_wgs_metrics.txt"
    out_log = out_parse + "_wgs_metrics.log"
    out_err = out_parse + "_wgs_metrics.err"
    wgs_metrics_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar /data/commonOwl/bin/picard-tools_archive/picard-tools-1.141/picard.jar CollectWgsMetrics TMP_DIR=/media/walllab/zhanna/tmp INPUT=" + in_bam + " OUTPUT=" + out_metrics + " REFERENCE_SEQUENCE=" + ref_genome + " COVERAGE_CAP=1000" + " INCLUDE_BQ_HISTOGRAM=true" + " COUNT_UNPAIRED=true" + " 1>" + out_log + " 2>" + out_err
    return wgs_metrics_comm

def build_index(in_bam):
    out_parse = in_bam[:-1] #take off the m
    outfile = out_parse + "i"
    out_log = out_parse + "i.log"
    out_err = out_parse + "i.err"
    build_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $PICARD BuildBamIndex MAX_RECORDS_IN_RAM=200000 TMP_DIR=/media/walllab/zhanna/tmp INPUT=" + in_bam + " OUTPUT=" + outfile + " 1>" + out_log + " 2>" + out_err
    return build_comm

def first_gatk(in_bam):
    out_parse = in_bam[:-4] #take off the .bam
    out_gvcf = out_parse + ".raw.snps.indels.g.vcf"
    out_log = out_parse + ".g.vcf.log"
    out_err = out_parse + ".g.vcf.err"
    first_gvcf_comm = "java -Xmx40g -Djava.io.tmpdir=/media/walllab/zhanna/tmp -jar $GATK -T HaplotypeCaller -R "+ ref_genome + " -I " + in_bam + " --emitRefConfidence GVCF -o " + out_gvcf + " 1>" + out_log + " 2>" + out_err
    return out_gvcf, first_gvcf_comm

ref_genome = '/media/walllab/zhanna/owl/assemblies/StrOccCau_1.0_nuc_finalMito.fa'
date = "2018Feb02"
aln1 = make_aln_comm2(samp1, sampfull1, runID1, lib1, lane1, insert1, maxinsert1, pair1, unpair1)
aln2 = make_aln_comm2(samp2, sampfull2, runID2, lib2, lane2, insert2, maxinsert2, pair2, unpair2)
aln3 = make_aln_comm2(samp3, sampfull3, runID3, lib3, lane3, insert3, maxinsert3, pair3, unpair3)
aln7 = make_aln_comm2(samp7, sampfull7, runID7, lib7, lane7, insert7, maxinsert7, pair7, unpair7)
aln8 = make_aln_comm2(samp8, sampfull8, runID8, lib8, lane8, insert8, maxinsert8, pair8, unpair8)
aln9 = make_aln_comm2(samp9, sampfull9, runID9, lib9, lane9, insert9, maxinsert9, pair9, unpair9)
aln11 = make_aln_comm2(samp11, sampfull11, runID11, lib11, lane11, insert11, maxinsert11, pair11, unpair11)

dedup_ls = [aln1,aln2,aln3,aln7,aln8,aln9,aln11]

merged_dedups = merge_ls_bams(dedup_ls)
print merged_dedups[1]
Run(merged_dedups[1])

sorted_merged = sort_sam(merged_dedups[0])
print sorted_merged[1]
Run(sorted_merged[1])
dedup_merged = mark_dups(sorted_merged[0])
print dedup_merged[1]
Run(dedup_merged[1])
build = build_index(dedup_merged[0])
print build
Run(build)
gvcf_1 = first_gatk(dedup_merged[0])
print gvcf_1[1]
Run(gvcf_1[1])
