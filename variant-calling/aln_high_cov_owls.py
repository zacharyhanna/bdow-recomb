import sys
import mimetypes
import gzip
import subprocess
# index Assembly
# install GATK

def Run(command): #run the command in the bash shell
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
date = '2018Feb02' #time.strftime("%Y%b%d") date of trimming to add to end of trim files - usually when running script

raw_path = "/home/zhanna/owl/raw_reads/" # parent directory of all raw read files
filt_path = "/home/zhanna/owl/filtered_reads/" # path to parent directory for all filtered reads
paired = "prd"
unpaired = "unprd"

run_ls = ["owl_2017Aug01_MedGenome"]

def getTrimDate(run):
    if run == "owl_2017Aug01_MedGenome":
        trimdate = '2017Sep13'
    elif run == "owl_2013Oct24_UCSF":
        trimdate = '2015Mar19'
    return trimdate


def getDataDir(Fullrun):
    DirName = raw_path + Fullrun + "/"
    return DirName

def getRunName(Fullrun):
    runName = Fullrun.split("_")[1:]
    run = ("_").join(runName)
    return run

def getRunDate(Fullrun):
    runDate = Fullrun.split("_")[1]
    return runDate


def getSampList(Fullrun):
    sampFile = getDataDir(Fullrun) + Fullrun + "_samplesRun2"
    sampList = []
    with open(sampFile, 'r') as file:
        for line in file:
            line = line.strip()
            sampList.append(line)
    return sampList

def getPhred(Fullrun):
    phredFile = getDataDir(Fullrun) + Fullrun + "_phred"
    with open(phredFile, 'r') as file:
        phred = file.readline()
        phred = str(phred)
        return phred

def getAdapPath(samp, Fullrun):
    adapterDir = "adapters/"
    adapter_ext = "_adapters.fa"
    adapter_path = "/home/zhanna/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa" #getDataDir(Fullrun) + adapterDir + samp + "_" + getRunName(Fullrun) + adapter_ext
    return adapter_path

def runTrim(samp, Fullrun): #takes a list of samples and runs a paired-end trim for each
    read1 = "R1" # what follows the sample name in the full name of all raw read1 files
    read2 = "R2" # what follows the sample name in the full name of all raw read2 files
    final_ext = ".fastq.gz"
    trimmed_samps = []
    filt_out_paired1 = filt_path + samp + "_" + getRunName(Fullrun) + "_" + read1 + "_" + getTrimDate(Fullrun) + trimtype + "_" + paired + final_ext
    filt_out_unpaired1 = filt_path + samp + "_" + getRunName(Fullrun) + "_" + read1 + "_" + getTrimDate(Fullrun) + trimtype + "_" + unpaired + final_ext
    filt_out_paired2 = filt_path + samp + "_" + getRunName(Fullrun) + "_" + read2 + "_" + getTrimDate(Fullrun) + trimtype + "_" + paired + final_ext
    filt_out_unpaired2 = filt_path + samp + "_" + getRunName(Fullrun) + "_" + read2 + "_" + getTrimDate(Fullrun) + trimtype + "_" + unpaired + final_ext
    trimmed_samps.append(filt_out_paired1)
    trimmed_samps.append(filt_out_unpaired1)
    trimmed_samps.append(filt_out_paired2)
    trimmed_samps.append(filt_out_unpaired2)
    return trimmed_samps

def ListTrim(run):
    sampList = getSampList(run) # retrieving the list of samples put on the run, using the samples file in the run raw_reads directory to get the list of samples. output is a list, assigned to "sampList"
    trim_out = runTrim(sampList, run)
    return trim_out

def SampLsAdd(trim_out, list):
    for samp in trim_out[0]:
        list.append(samp)
    return

def RawLsAdd(trim_out):
    for raw in trim_out[1]:
        raw_samps_ls.append(raw)
    return

def get_insert(samp, run):
    insertFile = getDataDir(run) + run + "_insert_sizes"
    with open(insertFile, 'r') as file:
        for line in file:
            line = line.strip()
            splitline = line.split('\t')
            if splitline[0] == samp:
                insert = str(splitline[1])
                return insert

def aln_comm(samp, run, paired):
    maxinsert = "1000"
    global date
    aln_dir_comm = "mkdir -p /home/zhanna/owl/filtered_reads/" + samp
    print aln_dir_comm
    Run(aln_dir_comm)
    infile_ls = runTrim(samp, run)
    if paired == 1: #paired-end files used
        aln_file = "/home/zhanna/owl/filtered_reads/" + samp + "/" + samp + "_" + getRunName(run) + "_" + getTrimDate(run) + trimtype + "_" + date + "alnStrOccCau1_mates.sam"
        aln_log = "/home/zhanna/owl/filtered_reads/" + samp + "/" + samp + "_" + getRunName(run) + "_" + getTrimDate(run) + trimtype + "_" + date + "alnStrOccCau1_mates.log"
        aln_comm = "bwa mem -M -t 36 -R '@RG\\tID:" + samp + run + "ln1\\tSM:" + samp + "\\tPL:illumina\\tLB:" + run + "lib\\tPU:lane1\\tPI:420' -w " + maxinsert + " " + ref_genome + " " + infile_ls[0] + " " + infile_ls[2] + " > " + aln_file + " 2> " + aln_log
    elif paired == 0: #single reads being aligned
        aln_file = "/home/zhanna/owl/filtered_reads/" + samp + "/" + samp + "_" + getRunName(run) + "_" + getTrimDate(run) + trimtype + "_" + date + "alnStrOccCau1_singles.sam"
        aln_log = "/home/zhanna/owl/filtered_reads/" + samp + "/" + samp + "_" + getRunName(run) + "_" + getTrimDate(run) + trimtype + "_" + date + "alnStrOccCau1_singles.log"
        aln_comm = "zcat " + infile_ls[1] + " " + infile_ls[3] + " | bwa mem -M -t 36 -R '@RG\\tID:" + samp + run + "ln1\\tSM:" + samp + "\\tPL:illumina\\tLB:" + run + "lib\\tPU:lane1\\tPI:420' " + ref_genome + " - > " + aln_file + " 2> " + aln_log
    return aln_file, aln_comm

def merge_sams(sam1, sam2):
    out_parse = sam1.split("_")
    out_ls = out_parse[:-1]
    outfile_bs = ("_").join(out_ls)
    outfile = outfile_bs + "_mates_singles.sam"
    out_log = outfile_bs + "_mates_singles.log"
    out_err = outfile_bs + "_mates_singles.err"
    merge_comm = "java -Xmx100g -Djava.io.tmpdir=/home/zhanna/tmp -jar $PICARD MergeSamFiles MAX_RECORDS_IN_RAM=200000 INPUT=" + sam1 + " INPUT=" + sam2 + " OUTPUT=" + outfile + " SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=/home/zhanna/tmp 1>" + out_log + " 2>" + out_err
    return outfile, merge_comm

def sort_sam(in_sam):
    out_parse = in_sam[:-4] #take off the .sam
    out_bam = out_parse + "_sorted.bam"
    out_log = out_parse + "_sorted.log"
    out_err = out_parse + "_sorted.err"
    sort_comm = "java -Xmx100g -Djava.io.tmpdir=/home/zhanna/tmp -jar $PICARD SortSam MAX_RECORDS_IN_RAM=200000 INPUT=" + in_sam + " OUTPUT=" + out_bam + " SORT_ORDER=coordinate TMP_DIR=/home/zhanna/tmp 1>" + out_log + " 2>" + out_err
    return out_bam, sort_comm

def mark_dups(in_bam):
    out_parse = in_bam[:-4] #take off the .bam
    out_bam = out_parse + "_dedup.bam"
    out_log = out_parse + "_dedup.log"
    out_err = out_parse + "_dedup.err"
    out_metrics = out_parse + "_dedup_metrics.txt"
    mark_comm = "java -Xmx100g -Djava.io.tmpdir=/home/zhanna/tmp -jar $PICARD MarkDuplicates MAX_RECORDS_IN_RAM=200000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=/home/zhanna/tmp INPUT=" + in_bam + " OUTPUT=" + out_bam + " METRICS_FILE=" + out_metrics + " 1>" + out_log + " 2>" + out_err
    return out_bam, mark_comm

def build_index(in_bam):
    out_parse = in_bam[:-1] #take off the m
    outfile = out_parse + "i"
    out_log = out_parse + "i.log"
    out_err = out_parse + "i.err"
    build_comm = "java -Xmx100g -Djava.io.tmpdir=/home/zhanna/tmp -jar $PICARD BuildBamIndex MAX_RECORDS_IN_RAM=200000 TMP_DIR=/home/zhanna/tmp INPUT=" + in_bam + " OUTPUT=" + outfile + " 1>" + out_log + " 2>" + out_err
    return build_comm

def first_gatk(in_bam):
    out_parse = in_bam[:-4] #take off the .bam
    out_gvcf = out_parse + ".raw.snps.indels.g.vcf"
    out_log = out_parse + ".g.vcf.log"
    out_err = out_parse + ".g.vcf.err"
    first_gvcf_comm = "java -Xmx100g -Djava.io.tmpdir=/home/zhanna/tmp -jar $GATK -T HaplotypeCaller -R "+ ref_genome + " -I " + in_bam + " --emitRefConfidence GVCF -o " + out_gvcf + " 1>" + out_log + " 2>" + out_err
    return out_gvcf, first_gvcf_comm

def aln_and_filter(run):
    sampList = ["ZRHG101","ZRHG102","ZRHG103","ZRHG104","ZRHG105","ZRHG106","ZRHG107","ZRHG108","ZRHG109","ZRHG110","ZRHG111","ZRHG112","ZRHG113","ZRHG114","ZRHG115","ZRHG116","ZRHG117","ZRHG118","ZRHG119","ZRHG120","ZRHG121","ZRHG122","ZRHG123","ZRHG124","ZRHG125","ZRHG126","ZRHG127","ZRHG128"] # retrieving the list of samples put on the run, using the samples file in the run raw_reads directory to get the list of samples. output is a list, assigned to "sampList"
    for samp in sampList:
        aln_pairs = aln_comm(samp, run, 1) #paired reads alignment
        print aln_pairs[1]
        Run(aln_pairs[1])
        aln_singles = aln_comm(samp, run, 0) #single reads alignment
        print aln_singles[1]
        Run(aln_singles[1])
        merge_two = merge_sams(aln_pairs[0], aln_singles[0])
        print merge_two[1]
        Run(merge_two[1])
        sorted = sort_sam(merge_two[0])
        print sorted[1]
        Run(sorted[1])
        dedup = mark_dups(sorted[0])
        print dedup[1]
        Run(dedup[1])
        build = build_index(dedup[0])
        print build
        Run(build)
        gvcf_1 = first_gatk(dedup[0])
        print gvcf_1[1]
        #print "--variant " + str(gvcf_1[0])
        Run(gvcf_1[1])

sorted_aln_list = []

ref_genome = '/home/zhanna/owl/assemblies/StrOccCau_1.0_nuc_finalMito.fa'
trimtype = "AdptTrim" #type of trimming are performing AdptTrim = only adapter trimming, trim = adapter trimming and quality trimming
quality = ''
for run in run_ls:
    aln_and_filter(run)
