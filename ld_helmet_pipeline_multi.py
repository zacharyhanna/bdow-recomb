import subprocess
import sys
import glob
from functools import partial
from multiprocessing.dummy import Pool
"""
Example run:
python ld_helmet_pipeline.py <run || test> scaffold_list_file theta_list_file number_haplotypes total_threads
"""
total_threads = int(sys.argv[5])
number_haplotypes = int(sys.argv[4])
output_dir_stem = "/media/walllab/zhanna/owl/ldhelmet/"
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
    out_prefix = output_dir_stem
    output = out_prefix + scaf + "_BADOeastGrEq1Mb"
    cmd = "vcftools --gzvcf " + vcf_input + " --chr " + scaf + " --recode --out " + output
    return cmd, 1

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
    # remove old fasta
    rm_fasta_cmd = "rm " + fasta
    return cmd, rm_fasta_cmd

def rename_fastas(scaf):
    comm_dict = {"rename_sample_com": [], "rm_fasta_cmd1": [], "rm_fasta_cmd2": []}
    fasta_dir = output_dir_stem
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
        comm_dict["mkdir_cmd"] = mkdir_cmd
        rename_sample_com, rm_fasta_cmd1 = rename_sample_seq(fasta, samp, haplo, new_fasta)
        comm_dict["rename_sample_com"].append(rename_sample_com)
        comm_dict["rm_fasta_cmd1"].append(rm_fasta_cmd1)
    combine_fastas_cmd = "cat " + fasta_dir + scaf + "/*.fasta >" + combo_fasta
    comm_dict["combine_fastas_cmd"] = combine_fastas_cmd
    # remove the individual sample fastas after combining them into the total fasta
    for mynew_fasta in new_fasta_ls:
        rm_fasta_cmd2 = "rm " + mynew_fasta
        comm_dict["rm_fasta_cmd2"].append(rm_fasta_cmd2)
    return combo_fasta, comm_dict, 1

def vcf2fasta(scaf):
    reference_dir = "/media/walllab/zhanna/owl/assemblies/StrOccCau2_scafs/"
    reference = reference_dir + scaf + "_StrOccCau_2.0.fa"
    input = output_dir_stem + scaf + "_BADOeastGrEq1Mb.recode.vcf"
    output = output_dir_stem
    combined_output = scaf + "_BADOeastGrEq1Mb.fa"
    cmd = "vcf2fasta --reference " + reference + " --prefix " + output + " " + input
    return cmd

def ldhelmet_find_confs(scaf, combined_fasta):
    input = combined_fasta
    # take off "fasta" extension
    output = combined_fasta[:-5] + "conf"
    find_confs_cmd = "~/bin/LDhelmet_v1.10/ldhelmet find_confs --num_threads 6 -w 50 -o " + output + " " + input
    return output, find_confs_cmd, 6

def get_thetas(file_in):
    theta_dict = {}
    with open(file_in, 'r') as infile:
        for line in infile:
            splitline = line.strip().split()
            # first column is the scaffold
            # second column is the theta value
            theta_dict[splitline[0]] = str(splitline[1])
    return theta_dict

def activate_ldpop_envir():
    activate_ldpop = "source activate LDpopPy3"
    return activate_ldpop

def deactivate_ldpop_envir():
    deactivate_ldpop = "source deactivate"
    return deactivate_ldpop

def ldpop(scaf, combined_fasta, scaf_theta):
    # take off "fasta" extension
    output = combined_fasta[:-5] + "ldpop"
    ldpop_cmd = "~/bin/ldpop/run/ldtable.py -n " + number_haplotypes + " -th " + scaf_theta \
    + " -rh 101,100 --approx --cores 6 >" + output
    return output, ldpop_cmd, 6

def ldhelmet_convert_table(scaf, input_ldpop, input_conf):
    # take off "ldpop" extension
    output = input_ldpop[:-5] + "lk"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet convert_table --input_table " \
    + input_ldpop + " --output_table " + output + " --config_file " + input_conf
    return output, cmd

def ldhelmet_pade(scaf, scaf_theta, input_conf):
    pade_coeff = "11"
    # take off "conf" extension
    output = input_conf[:-4] + "pade"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet pade --num_threads 6 -t " \
    + scaf_theta + " -x " + pade_coeff + " -c " + input_conf + " -o " + output
    return output, cmd, 6

def ldhelmet_rjmcmc(scaf, input_lk, input_pade, input_fasta):
    # block penalty is a parameter that can be optimized
    block_penalty = "50.0"
    # chain steps in burn-in
    burn_in = "100000"
    # total iterations of chain after burn-in
    total_iter = "1000000"
    # take off "fasta" extension
    post_output = input_fasta[:-5] + "post"
    cmd = "~/bin/LDhelmet_v1.10/ldhelmet rjmcmc --num_threads 6 -l " \
    + input_lk + " -p " + input_pade + " -s "+ input_fasta + " -b " \
    + block_penalty + " --burn_in " + burn_in + " -n " + total_iter + " -o " + post_output
    return post_output, cmd, 6

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
    return cmd

def run_ls_comm(comm_list, threads):
    num_comm = total_threads // threads
    print(num_comm)
    print(comm_list)
    if realrun:
        pool = Pool(num_comm) # two concurrent commands at a time
        for i, returncode in enumerate(pool.imap(partial(subprocess.call, shell=True), comm_list)):
            if returncode != 0:
                print("%d command failed: %d" % (i, returncode))

def run_pipeline(scaf_file):
    scaf_ls = get_scafs(scaf_file)
    theta_dict = get_thetas(sys.argv[3])
    get_vcf_data_dict = {"com_ls":[], "threads":1}
    vcf2fasta_dict = {"com_ls":[], "threads":1}
    rename_fastas_dict1 = {"com_ls":[], "threads":1}
    rename_fastas_dict2 = {"com_ls":[], "threads":1}
    rename_fastas_dict3 = {"com_ls":[], "threads":1}
    rename_fastas_dict4 = {"com_ls":[], "threads":1}
    rename_fastas_dict5 = {"com_ls":[], "threads":1}
    ldhelmet_find_confs_dict = {"com_ls":[], "threads":1}
    ldpop_dict = {"com_ls":[], "threads":1}
    ldpop_dict["com_ls"].append(activate_ldpop_envir())
    ldpop_deactivate = {"com_ls":[deactivate_ldpop_envir()], "threads":1}
    ldhelmet_convert_table_dict = {"com_ls":[], "threads":1}
    ldhelmet_pade_dict = {"com_ls":[], "threads":1}
    ldhelmet_rjmcmc_dict = {"com_ls":[], "threads":1}
    ldhelmet_post_to_text_dict = {"com_ls":[], "threads":1}
    for scaf in scaf_ls:
        scaf_theta = theta_dict[scaf]
        vcfcom, vcfthreads = get_vcf_data(scaf)
        get_vcf_data_dict["com_ls"].append(vcfcom)
        get_vcf_data_dict["threads"] = vcfthreads
        vcf2fasta_dict["com_ls"].append(vcf2fasta(scaf))
        combo_fasta, rename_fastas_comm_dict, rename_fastas_threads = rename_fastas(scaf)
        rename_fastas_dict1["com_ls"].append(rename_fastas_comm_dict["mkdir_cmd"])
        rename_fastas_dict2["com_ls"] = rename_fastas_dict2["com_ls"] + rename_fastas_comm_dict["rename_sample_com"]
        rename_fastas_dict3["com_ls"] = rename_fastas_dict3["com_ls"] + rename_fastas_comm_dict["rm_fasta_cmd1"]
        rename_fastas_dict4["com_ls"].append(rename_fastas_comm_dict["combine_fastas_cmd"])
        rename_fastas_dict5["com_ls"] = rename_fastas_dict5["com_ls"] + rename_fastas_comm_dict["rm_fasta_cmd2"]
        conf_file, find_confs_com, find_confs_com_threads = ldhelmet_find_confs(scaf, combo_fasta)
        ldhelmet_find_confs_dict["com_ls"].append(find_confs_com)
        ldhelmet_find_confs_dict["threads"] = find_confs_com_threads
        ldpop_file, ldpop_com, ldpop_threads = ldpop(scaf, combo_fasta, scaf_theta)
        ldpop_dict["com_ls"].append(ldpop_com)
        ldpop_dict["threads"] = ldpop_threads
        lk_file, ld_convert_com = ldhelmet_convert_table(scaf, ldpop_file, conf_file)
        ldhelmet_convert_table_dict["com_ls"].append(ld_convert_com)
        pade_file, pade_com, pade_threads = ldhelmet_pade(scaf, scaf_theta, conf_file)
        ldhelmet_pade_dict["com_ls"].append(pade_com)
        ldhelmet_pade_dict["threads"] = pade_threads
        post_file, ldhelmet_rjmcmc_com, ldhelmet_rjmcmc_threads = ldhelmet_rjmcmc(scaf, lk_file, pade_file, combo_fasta)
        ldhelmet_rjmcmc_dict["com_ls"].append(ldhelmet_rjmcmc_com)
        ldhelmet_rjmcmc_dict["threads"] = ldhelmet_rjmcmc_threads
        ldhelmet_post_to_text_com = ldhelmet_post_to_text(scaf, post_file)
        ldhelmet_post_to_text_dict["com_ls"].append(ldhelmet_post_to_text_com)
    ls_filled_dicts = [\
    get_vcf_data_dict, \
    vcf2fasta_dict, \
    rename_fastas_dict1, \
    rename_fastas_dict2, \
    rename_fastas_dict3, \
    rename_fastas_dict4, \
    rename_fastas_dict5, \
    ldhelmet_find_confs_dict, \
    ldpop_dict, \
    ldpop_deactivate, \
    ldhelmet_convert_table_dict, \
    ldhelmet_pade_dict, \
    ldhelmet_rjmcmc_dict, \
    ldhelmet_post_to_text_dict]
    for filled_comm_dict in ls_filled_dicts:
        run_ls_comm(filled_comm_dict["com_ls"], filled_comm_dict["threads"])

if sys.argv[1] == "run":
    realrun = True
else:
    realrun = False
run_pipeline(sys.argv[2])
