#!/bin/bash
#
#SBATCH -J bcf_biallelic                                                        #job name
#SBATCH -N 1                                                                    #ensure that all cores are on one node
#SBATCH -n 4                                                                    #number of cores
#SBATCH -t 1-0                                                                  #runtime in D-HH:MM
#SBATCH --mem 2G                                                                #memory in GB
#SBATCH -o bcf_chr20.vcf                                         #file for STDOUT
#SBATCH -e bcf_chr20.err                                          #file for STDERR
#SBATCH --mail-type=ALL                                                 #type of email notification: BEGIN,END,FAIL
#SBATCH --mail-user=yourUsername@brown.edu              #email where notifications will be sent
#SBATCH --input=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Your command should be formatted like this:
# sbatch bcf_biallelic.sh

###Load the bcftools module and then execute the filtering on the desired file
module load bcftools
# print("bcftools module loaded")
bcftools view -m2 -M2 -i VT="SNP" ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o bcf_chr20.vcf.gz
# NEW: Included -I to skip indels
# bcftools view -m2 -M2 -v -I snps ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
