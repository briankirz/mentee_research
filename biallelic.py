import sys
import gzip
import numpy as np
import re


def biallelic(vcf_input, vcf_output):
    """
    Takes in input filepath and output filepath. Unzips the input,
    and returns a NumPy array.

    :param vcf_input: filepath for inputs, ex. 'ALL...vcf.gz'
    :param vcf_output: filepath for outputs, ex. '.vcf'
    """

    # unzips input file, gives write permissions to output file
    with gzip.open(vcf_input, 'rb') as gz, open(vcf_output, 'w') as f:
        # header boolean keeps track of when header is passed
        header = True
        line_num = 0
        for line in gz:
            # Indicates continued progress on processing large chromosomes
            if line_num % 10000 == 0:
                print("Running line " + str(line_num))
            # Makes sure to write the header exactly as is. Once we are no longer in the header this short-circuits
            if header and (line[0] == "#"):
                f.write(line)
                line_num += 1
            else:
                # Indicates we are no longer in the header
                header = False
                # Splits the line by tabs (column) and stores it in a list
                split = line.split('\t')
                # Less than 2 alternative alleles ?
                # Site not monomorphic?
                # gating with these two booleans to save speed and cut down on redundant work
                if len(split[4]) == 1 and split[4] != '.':
                    # If VT=SNP is found in INFO, we can write to the output file
                    if 'VT=SNP' in split[7]:
                        f.write(line)
                line_num += 1

    print("Biallelic filtering complete")
    return None


biallelic(sys.argv[1], sys.argv[2])
