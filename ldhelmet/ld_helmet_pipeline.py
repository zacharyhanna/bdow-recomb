import subprocess
import sys
import glob
"""
Example run:
python ld_helmet_pipeline.py <run || test> scaffold_list_file theta_list_file
"""

def Run(command): #run the command in the bash shell
    print(command)
    # realrun is a global variable specified by sys.argv[2]
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
    " -c fastx \'{print \">\"$name\"_\"mysamp\"_\"myhaplo\"\\n\"$seq}\'" + " " \
    + fasta + " >" + new_fasta
    Run(cmd)
    # remove old fasta
    rm_fasta_cmd = "rm " + fasta
    Run(rm_fasta_cmd)

def rename_fastas(scaf):
    fasta_dir = "/media/walllab/zhanna/owl/ldhelmet/"
    fasta_pattern = fasta_dir + "*_" + scaf + ":*.fasta"
    fasta_ls = glob.glob(fasta_pattern)
    new_fasta_ls = []
    combo_fasta = fasta_dir + scaf + "/" + scaf + "_BADOeast_haplos.fasta"
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
        new_fasta_ls.append(new_fasta)
        mkdir_cmd = "mkdir -p " + new_fasta_dir
        Run(mkdir_cmd)
        rename_sample_seq(fasta, samp, haplo, new_fasta)
    combine_fastas_cmd = "cat " + fasta_dir + scaf + "/*.fasta >" + combo_fasta
    Run(combine_fastas_cmd)
    # remove the individual sample fastas after combining them into the total fasta
    for mynew_fasta in new_fasta_ls:
        rm_fasta_cmd = "rm " + mynew_fasta
        Run(rm_fasta_cmd)
    return combo_fasta

def vcf2fasta(scaf):
    reference_dir = "/media/walllab/zhanna/owl/assemblies/StrOccCau2_scafs/"
    reference = reference_dir + scaf + "_StrOccCau_2.0.fa"
    input = "/media/walllab/zhanna/owl/ldhelmet/" + scaf + "_BADOeastGrEq1Mb.recode.vcf"
    output = "/media/walllab/zhanna/owl/ldhelmet/"
    combined_output = scaf + "_BADOeastGrEq1Mb.fa"
    cmd = "vcf2fasta --reference " + reference + " --prefix " + output + " " + input
    Run(cmd)

def ldhelmet_find_confs(scaf, combined_fasta):
    input = combined_fasta
    # take off "fasta" extension
    output = combined_fasta[:-5] + "conf"
    find_confs_cmd = "~/bin/LDhelmet_v1.10/ldhelmet find_confs --num_threads 12 -w 50 -o " + output + " " + input
    Run(find_confs_cmd)
    return output

def get_thetas(file_in):
    theta_dict = {}
    with open(file_in, 'r') as infile:
        for line in infile:
            splitline = line.strip().split()
            # first column is the scaffold
            # second column is the theta value
            theta_dict[splitline[0]] = str(splitline[1])
    return theta_dict

def ldpop(scaf, combined_fasta, scaf_theta):
    activate_ldpop = "source activate LDpopPy3"
    deactivate_ldpop = "source deactivate"
    Run(activate_ldpop)
    # take off "fasta" extension
    output = combined_fasta[:-5] + "ldpop"
    ldpop_cmd = "~/bin/ldpop/run/ldtable.py -n 24 -th " + scaf_theta \
    + " -rh 101,100 --approx --cores 12 >" + output
    Run(ldpop_cmd)
    Run(deactivate_ldpop)
    return output

def ldhelmet_convert_table(scaf, input_ldpop, input_conf):
    # take off "ldpop" extension
    output = input_ldpop[:-5] + "lk"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet convert_table --input_table " \
    + input_ldpop + " --output_table " + output + " --config_file " + input_conf
    Run(cmd)
    return output

def ldhelmet_pade(scaf, scaf_theta, input_conf):
    pade_coeff = "11"
    # take off "conf" extension
    output = input_conf[:-4] + "pade"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet pade --num_threads 24 -t " \
    + scaf_theta + " -x " + pade_coeff + " -c " + input_conf + " -o " + output
    Run(cmd)
    return output

def ldhelmet_rjmcmc(scaf, input_lk, input_pade, input_fasta):
    # block penalty is a parameter that can be optimized
    block_penalty = "50.0"
    # chain steps in burn-in
    burn_in = "100000"
    # total iterations of chain after burn-in
    total_iter = "1000000"
    # take off "fasta" extension
    post_output = input_fasta[:-5] + "post"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet rjmcmc --num_threads 24 -l " \
    + input_lk + " -p " + input_pade + " -s "+ input_fasta + " -b " \
    + block_penalty + " --burn_in " + burn_in + " -n " + total_iter + " -o " + post_output
    Run(cmd)
    return post_output

def ldhelmet_post_to_text(scaf, input_post):
    low_percentile = "0.025"
    mid_percentile = "0.5"
    high_percentile = "0.975"
    # take off "post" extension
    output = input_post[:-4] + "txt"
    # -m option outputs the mean
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet post_to_text -m -p " + low_percentile \
    + " -p " + mid_percentile + " -p " + high_percentile + " -o " + output \
    + " " + input_post
    Run(cmd)

def run_pipeline(scaf_file):
    scaf_ls = get_scafs(scaf_file)
    theta_dict = get_thetas(sys.argv[3])
    for scaf in scaf_ls:
        scaf_theta = theta_dict[scaf]
        get_vcf_data(scaf)
        vcf2fasta(scaf)
        combo_fasta = rename_fastas(scaf)
        conf_file = ldhelmet_find_confs(scaf, combo_fasta)
        ldpop_file = ldpop(scaf, combo_fasta, scaf_theta)
        lk_file = ldhelmet_convert_table(scaf, ldpop_file, conf_file)
        pade_file = ldhelmet_pade(scaf, scaf_theta, conf_file)
        post_file = ldhelmet_rjmcmc(scaf, lk_file, pade_file, combo_fasta)
        ldhelmet_post_to_text(scaf, post_file)

if sys.argv[1] == "run":
    realrun = True
else:
    realrun = False
run_pipeline(sys.argv[2])
