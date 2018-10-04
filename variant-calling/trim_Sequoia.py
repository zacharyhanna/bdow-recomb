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

def changeRun(Fullrun):
    diffrun = run_dict[Fullrun]
    return diffrun

def getAdapPath(Fullrun):
    adapter_ext = "_adapters.fa"
    changed_run = changeRun(Fullrun)
    adapter_path = raw_path + changed_run + adapter_ext
    return adapter_path


def runTrim(run): #takes a list of samples and runs a paired-end trim for each
    read1 = "1" # what follows the sample name in the full name of all raw read1 files
    read2 = "2" # what follows the sample name in the full name of all raw read2 files
    final_ext = ".fastq"
    in_file1 = raw_path + run + "_" + read1 + final_ext
    in_file2 = raw_path + run + "_" + read2 + final_ext

    filt_out_paired1 = filt_path + run + "_" + read1 + "_" + date + trimtype + "_" + paired + final_ext
    filt_out_unpaired1 = filt_path + run + "_" + read1 + "_" + date + trimtype + "_" + unpaired + final_ext
    filt_out_paired2 = filt_path + run + "_" + read2 + "_" + date + trimtype + "_" + paired + final_ext
    filt_out_unpaired2 = filt_path + run + "_" + read2 + "_" + date + trimtype + "_" + unpaired + final_ext

    adapPath = getAdapPath(run)

    trimlog_path = filt_path + run + "_" + date + trimtype + "_trimlog"
    trim_com = "java -jar /home/zhanna/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 42 " + "-trimlog "+trimlog_path+" "+in_file1+" "+in_file2+" "+filt_out_paired1+" "+filt_out_unpaired1+" "+filt_out_paired2+" "+filt_out_unpaired2+" ILLUMINACLIP:" + adapPath + ":2:30:10 LEADING:3 TRAILING:3 MINLEN:36"
    print trim_com
    Run(trim_com)
    return
date = "2017Feb07"#time.strftime("%Y%b%d") #date of trimming to add to end of trim files - usually when running script
trimtype = "AdptTrim" #type of trimming are performing AdptTrim = only adapter trimming, trim = adapter trimming and quality trimming
raw_path = "/media/walllab/zhanna/owl/Sequoia_raw/" # parent directory of all raw read files
filt_path = "/media/walllab/zhanna/owl/Sequoia_filt/" # path to parent directory for all filtered reads
paired = "prd"
unpaired = "unprd"

run_dict = {
"SRR4011615":"Sequoia_550bpnoPCR_2014Sep10_CAS",
"SRR4011616":"Sequoia_900bpPCR_2014Sep10_CAS",
"SRR4011596":"Sequoia_2011Aug_Ill",
"SRR4011595":"Sequoia_2011Jul_UCSF",
"SRR4011597":"Sequoia_2013Oct07_UCB",
"SRR4011614":"Sequoia_2014May1224_CAS",
"SRR4011617":"Sequoia_Hydro6_2014Sep10_CAS"
}

run_ls = ["SRR4011595","SRR4011614","SRR4011615","SRR4011616","SRR4011617","SRR4011596","SRR4011597"]

for run in run_ls:
    runTrim(run)
