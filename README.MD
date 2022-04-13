biallelic.py: filters a gzipped vcf file into a new vcf file that only contains biallelic SNPs

bcf_biallelic.sh: shell script that does the same thing as biallelic.py but uses bcftools to filter

compare_vcf.py: returns whether two vcf files (such as those produced by filtering) are identical to each other