import sys
import os

# Input:
# str in_filename: name of gzipped vcf file we want to filter
# str out_filename: name of vcf file we want to write biallelic sites to
def bcf_biallelic(in_filename, out_filename):
        # preps the system by loading bcftools
        os.system('module load bcftools')
        # executes the biallelic command based on input
        os.system('bcftools view -m2 -M2 -v snps ' + in_filename + ' -o ' + out_filename)

bcf_biallelic(sys.argv[1], sys.argv[2])