import sys
import gzip
import numpy as np

def biallelic(vcf_input, vcf_output):
        """
        Takes in input filepath and output filepath. Unzips the input,
        and returns a NumPy array.

        Read the data of the file into a buffer and use np.frombuffer
        to turn the data into a NumPy array.

        :param vcf_in: filepath for inputs, ex. 'ALL...vcf.gz'
        :param vcf_out: filepath for outputs, ex. '.vcf'
        :return: NumPy array of vcf file
        """

        with gzip.open(vcf_input, 'rb') as gz, open(vcf_output, 'w') as f:
                header = True
                line_num = 0
                for line in gz:
                        if line_num % 10000 == 0:
                                print("Running line " + str(line_num))
                        if header and (line[0] == "#"):
                                f.write(line)
                                line_num += 1
                        else:
                                header = False
                                alt = line.split('\t')[4]
                                if (len(alt) == 1) and (alt != "."):
                                        f.write(line)
                                line_num += 1
        print("Biallelic filtering complete")
        return None


biallelic(sys.argv[1], sys.argv[2])
