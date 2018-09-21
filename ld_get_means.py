import sys
import glob
"""
Example run:
python get_means_ld.py scaffold_list_file means_output_file
"""

res_dir = "/media/walllab/zhanna/owl/ldhelmet/"

def get_scafs(file_in):
    scaf_ls = []
    with open(file_in, 'r') as infile:
        for line in infile:
            splitline = line.strip().split()
            scaf_ls.append(splitline[0])
    return scaf_ls

def get_ld_txt(scaf, res_path):
    res_pattern = res_path + scaf + "/" + scaf + "*haplos.txt"
    return res_pattern

def wtd_vals(splitline):
    dist = splitline[1]-splitline[0]
    wtd_mean = splitline[2]*dist
    return dist, wtd_mean


def get_mean(ld_res_file):
    tot_dist = 0
    tot_wtd_mean = 0.0
    with open(ld_res_file) as ld_file:
        line_num = 0
        for line in ld_file:
            line_num += 1
            if line_num > 3:
                splitline = line.strip().split()
                line_dist, line_wtd_mean = wtd_vals(splitline)
                tot_dist += line_dist
                tot_wtd_mean += line_wtd_mean
    return str(tot_wtd_mean / tot_dist)

def write_results(res_ls, out_file):
    with open(out_file, 'w') as outfile:
        for item_ls in res_ls:
            outfile.write('\t'.join(item_ls))
            outfile.write('\n')

def run_pipeline(scaf_file, out_file):
    scaf_ls = get_scafs(scaf_file)
    res_ls = []
    for scaf in scaf_ls:
        res_file = get_ld_txt(scaf, res_dir)
        scaf_mean = get_mean(res_file)
        res_ls.append([scaf, scaf_mean])
    write_results(res_ls, out_file)

run_pipeline(sys.argv[1], sys.argv[2])
