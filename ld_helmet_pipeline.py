
def get_vcf_data(scaf):
    input = "/media/walllab/saurabh/data/1000genomes/vcf/phase3_GRCh38_cleanedup_removeDuplicates/{regions}/{regions}.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.vcf"
    out_prefix = "/media/walllab/zhanna/owl/LDhelmet/"
    output = ""
    keep_file = "/media/walllab/zhanna/owl/LDhelmet/{params.region}_1KG_{params.numSeq}Seq.list"
    cmd = "vcftools --vcf " + input + " --keep " + keep_file + " --chr {params.chromosome} --recode --out {params.outprefix}
    Run(cmd)

def vcf2fasta(scaf):
    input =
    output =
    mkdir -p output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles
    /media/walllab/zhanna/bin/vcflib/bin/vcf2fasta --reference /media/walllab/saurabh/data/referenceGenomes/chr{params.chromosome}.fa --prefix output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles/ {input}
    cat output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles/* > {output}
    rm -rf output/{params.region}/{params.numSeq}sequences/chr{params.chromosome}/individualFastaFiles/
