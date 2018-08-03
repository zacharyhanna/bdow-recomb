import sys
"""
Script to estimate Watterson's Theta from ldhat input files.
Example:
python theta_watterson.py scaffold_list_file output_file_name
"""

def wattersons_theta_seq(seg_sites, num_seqs):
    """
    Return Watterson's Theta (per sequence)
    """
    a1 = sum([1.0/i for i in range(1, num_seqs)])
    return float(seg_sites) / a1

def wattersons_theta_site(seg_sites, num_seqs, len_seq):
    """
    Return Watterson's Theta (per site)
    """
    theta_seq = wattersons_theta_seq(seg_sites, num_seqs)
    return theta_seq / len_seq

def read_info(file_in):
    """
    Read a .locs file formatted for ldhat and
    return the number of segregating sites and the total length of the sequence.
    """
    with open(file_in, 'r') as infile:
        linenum = 0
        for line in infile:
            linenum +=1
            line = line.strip()
            splitline = line.split()
            num_variants = int(splitline[0])
            len_seq = 1 + (1000 * float(splitline[1])) #locs file gives the
            #length in kbp and also seems to subtract one from the total length,
            # so we add it back in here
            break
    return num_variants, len_seq

def get_list_scafs(file_in):
    """
    Return a list of scaffolds to walk through. Expects each scaffold on its own
    line in the input file.
    """
    scaf_ls = []
    with open(file_in, 'r') as infile:
        for line in infile:
            scaf_ls.append(line.strip())
    return scaf_ls

def write_output(file_out, theta_ls):
    """
    Write out the Watterson's theta values. Expects a list of tuples with the
    scaffold name as the first element in the tuple and the theta value as the
    second element.
    """
    with open(file_out, 'w') as outfile:
        for thetaVal in theta_ls:
            outStr = str(thetaVal[0]) + '\t' + str(thetaVal[1]) + '\n'
            outfile.write(outStr)

def get_list_thetas(scaf_file, num_haps):
    """
    Requires the file listing the scaffolds, one per line. Retrieves theta for
    each scaffold.
    """
    ls_scafs = get_list_scafs(scaf_file)
    ls_thetas = []
    for scaf in ls_scafs:
        locFile = "/media/walllab/zhanna/owl/ldhat/" + scaf + "_BADOeastGrEq1Mb.ldhat.locs"
        num_variants, length_seq = read_info(locFile)
        thetaVal = wattersons_theta_site(num_variants, num_haps, length_seq)
        ls_thetas.append((scaf, thetaVal))
    return ls_thetas

def mean_theta(ls_thetas):
    thetaSum = 0
    numThetas = 0
    for thetaVal in ls_thetas:
        thetaSum += float(thetaVal[1])
        numThetas += 1
    return thetaSum / numThetas

num_haplotypes = 24 # provide the number of haplotypes.
# It is 2 times the number of individuals.
list_thetas = get_list_thetas(sys.argv[1], num_haplotypes)
write_output(sys.argv[2], list_thetas)
meanTheta = "Mean theta: " + str(mean_theta(list_thetas)) + "\n"
print meanTheta
