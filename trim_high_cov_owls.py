import os
import subprocess
import time


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
date = time.strftime("%Y%b%d") #date of trimming to add to end of trim files - usually when running script


trimtype = "AdptTrim" #type of trimming are performing AdptTrim = only adapter trimming, trim = adapter trimming and quality trimming
#quality = ''

#from here down should be fine to leave unedited
raw_path = "/media/walllab/zhanna/owl/raw_reads/" # parent directory of all raw read files
filt_path = "/media/walllab/zhanna/owl/filtered_reads/" # path to parent directory for all filtered reads
paired = "prd"
unpaired = "unprd"

run_ls = ["owl_2017Aug01_MedGenome"]


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
    sampFile = getDataDir(Fullrun) + Fullrun + "_samples"
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
    adapter_path = getDataDir(Fullrun) + adapterDir + samp + "_" + getRunName(Fullrun) + adapter_ext
    return adapter_path


def runTrim(sampLs, Fullrun): #takes a list of samples and runs a paired-end trim for each
    read1 = "R1" # what follows the sample name in the full name of all raw read1 files
    read2 = "R2" # what follows the sample name in the full name of all raw read2 files
    final_ext = ".fastq.gz"
    for samp in sampLs:
        in_file1 = getDataDir(Fullrun) + samp + "_" + read1 + final_ext
        in_file2 = getDataDir(Fullrun) + samp + "_" + read2 + final_ext

        filt_out_paired1 = filt_path + samp + "/" + samp + "_" + getRunName(Fullrun) + "_" + read1 + "_" + date + trimtype + "_" + paired + final_ext
        filt_out_unpaired1 = filt_path + samp + "/" + samp + "_" + getRunName(Fullrun) + "_" + read1 + "_" + date + trimtype + "_" + unpaired + final_ext
        filt_out_paired2 = filt_path + samp + "/" + samp + "_" + getRunName(Fullrun) + "_" + read2 + "_" + date + trimtype + "_" + paired + final_ext
        filt_out_unpaired2 = filt_path + samp + "/" + samp + "_" + getRunName(Fullrun) + "_" + read2 + "_" + date + trimtype + "_" + unpaired + final_ext

        #adapPath = getAdapPath(samp, Fullrun)

        trimlog_path = filt_path + samp + "/" + samp + "_" + getRunName(Fullrun) + "_" + date + trimtype + "_trimlog"
        mkdir_com = "mkdir -p " + filt_path + samp
        print mkdir_com
        Run(mkdir_com)
    	trim_com = "java -jar /home/zhanna/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 48 " + "-trimlog "+trimlog_path+" "+in_file1+" "+in_file2+" "+filt_out_paired1+" "+filt_out_unpaired1+" "+filt_out_paired2+" "+filt_out_unpaired2+" ILLUMINACLIP:/home/zhanna/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36"
    	print trim_com
        Run(trim_com)
    return

for run in run_ls:
    sampList = getSampList(run)
    runTrim(sampList, run)
