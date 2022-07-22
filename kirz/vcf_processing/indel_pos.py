import sys
import gzip

def extract_indel_pos(chr_vcf_gz, out_txt):
	'''
	###################################################################################
	INPUT:
	a gzipped vcf file representing a chromosome
	--------------------------------------------
	OUPUT:
	a text file that contains the indel positions from the vcf file
	the chromosome reference number (ie 1, 2, 3, etc not chr1, chr2, chr3) in column 1
	tab-delineated
	the line number (reference number?) of a single indel
	###################################################################################
	'''
	# unzips input file, gives write permissions to output file
	with gzip.open(chr_vcf_gz, 'r') as gz, open(out_txt, 'w') as f:
		header = True
		line_num = 0
		#placeholder
		chrom = 0
		for line in gz:
			# Indicates continued progress on chromosome processing
			if line_num % 10000 == 0:
				print("Running line " + str(line_num))
			line_num += 1
			if header:
				if not (line[0] == '#'):
					header = False
			if not header:
				# No longer in the header
				header = False
				# Splits the line by tabs (columns) and stores each value in a list
				split = line.split()
				# Stores string data named after column it came from
				#if statement saves time on chromosome calculation
				if chrom == 0:
					chrom_index = 0
					chrom = split[chrom_index]
				pos_index = 1
				info_index = 7
				pos = split[pos_index]
				info = split[info_index]
				# If VT=INDEL is found in INFO, we write the chr reference and indel position to the output file
				if 'VT=INDEL' in info:
					f.write(chrom + '\t' + pos + '\n')
	print('Indel positions recorded for Chromosome ' + str(chrom) + ' complete')
	return None

extract_indel_pos(sys.argv[1], sys.argv[2])
