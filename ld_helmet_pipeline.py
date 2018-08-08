import subprocess
import sys
import glob
"""
Example run:
python ld_helmet_pipeline.py <run || test> scaffold_list_file
"""

def Run(command): #run the command in the bash shell
    print(command)
    if realrun:
        try:
            proc = subprocess.check_output(command, shell=True)
            return proc
        except subprocess.CalledProcessError:
            pass
        except:
            pass

def get_vcf_data(scaf):
    vcf_input = "/media/walllab/zhanna/owl/20180716_StrOccCau2_SpBarecal_snps_filt5_BADOeastGrEq1Mb.gt.vcf.gz"
    out_prefix = "/media/walllab/zhanna/owl/ldhelmet/"
    output = out_prefix + scaf + "_BADOeastGrEq1Mb"
    cmd = "vcftools --gzvcf " + vcf_input + " --chr " + scaf + " --recode --out " + output
    Run(cmd)

def get_scafs(file_in):
    scaf_ls = []
    with open(file_in, 'r') as infile:
        for line in infile:
            splitline = line.strip().split()
            scaf_ls.append(splitline[0])
    return scaf_ls

def rename_sample_seq(fasta, samp, haplo, new_fasta):
    cmd = "bioawk -v mysamp=" + samp + " -v myhaplo=" + haplo + \
    " -c fastx \'{print \">\"$name\"_$\{mysamp\}_$\{myhaplo\}\n\"$seq}\'" + " " \
    + fasta + " >" + new_fasta
    Run(cmd)

def rename_fastas(scaf):
    fasta_dir = "/media/walllab/zhanna/owl/ldhelmet/"
    fasta_pattern = fasta_dir + "*_" + scaf + ":*.fasta"
    fasta_ls = glob.glob(fasta_pattern)
    for fasta in fasta_ls:
        # split full file path
        splitfasta = fasta.split("/")
        # get just file name without path
        base_str = splitfasta[-1]
        # file name example: ZRHG117_Super-Scaffold_47_obj:1.fasta
        splitbase = base_str.split("_")
        samp = splitbase[0]
        haplosplit = base_str.split(":")
        haploside = haplosplit[1]
        haplosidesplit = haploside.split(".")
        haplo = haplosidesplit[0]
        new_fasta_dir = fasta_dir + scaf
        new_fasta = fasta_dir + scaf + "/" + scaf + "_" + samp + "_" + haplo + ".fasta"
        mkdir_cmd = "mkdir -p " + new_fasta_dir
        Run(mkdir_cmd)
        rename_sample_seq(fasta, samp, haplo, new_fasta)

    pass

def vcf2fasta(scaf):
    reference_dir = "/media/walllab/zhanna/owl/assemblies/StrOccCau2_scafs/"
    reference = reference_dir + scaf + "_StrOccCau_2.0.fa"
    input = "/media/walllab/zhanna/owl/ldhelmet/" + scaf + "_BADOeastGrEq1Mb.recode.vcf"
    output = "/media/walllab/zhanna/owl/ldhelmet/"
    combined_output = scaf + "_BADOeastGrEq1Mb.fa"
    cmd = "vcf2fasta --reference " + reference + " --prefix " + output + " " + input
    Run(cmd)
    #cat output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles/* > {output}
    #rm -rf output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles/

def run_pipeline(scaf_file):
    scaf_ls = get_scafs(scaf_file)
    for scaf in scaf_ls:
        #get_vcf_data(scaf)
        #vcf2fasta(scaf)
        rename_fastas(scaf)


if sys.argv[1] == "run":
    realrun = True
else:
    realrun = False
run_pipeline(sys.argv[2])
