import sys
import gzip

def compare_vcf(my_vcf_gz, bcf_vcf_gz):
	'''
	This function takes in two vcfs and compares them to see if they are the same.
	If they are, it returns True (vv)
	'''

	with open(my_vcf_gz, 'rb') as f1, open(bcf_vcf_gz, 'rb') as f2:
		data_line = 0
		# becomes untrue if mismatch found
		identical = True
		f1_data_index = 0
		f2_data_index = 0
		 # finds index of data start in file1
                for line1 in f1:
                        if line1[0] == "#":
                                f1_data_index += 1
                        else:
                                break
                # finds index of data start in file2
                for line2 in f2:
                        if line2[0] == "#":
                                f2_data_index += 1
                        else:
                                break

		
		# starting on data lines, check for identity on all following lines
		for line1 in f1:
			for line2 in f2:
				data_line += 1
				if line1 != line2:
					identical = False
					# print statements help debug
					print("Line " + str(data_line) + ": DIFFERENT")
					print("\tbcf: "+ f1[line1])
					print("\tpython: "+ f2[line2])
					# no need to continue when files are not identical
					break
		print(identical)
		f1.close()
		f2.close()

compare_vcf(sys.argv[1], sys.argv[2])
