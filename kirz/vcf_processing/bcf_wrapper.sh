#!/bin/bash
#
#SBATCH -J tgp_filter           # Job name
#SBATCH -N 1                    # Number of cores
#SBATCH -n 4                    # Ensure that all cores are on one machine
#SBATCH -t 6:00:00              # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 4g                # Memory in MB
#SBATCH -o logs/tgp_filter-%A.out # File for STDOUT (with jobid = %j)
#SBATCH -e logs/tgp_filter-%A.err # File for STDERR (with jobid = %j)
#SBATCH --mail-type=ALL         # Type of email notification: BEGIN,END,FAIL,A$
#SBATCH --mail-user=brian_kirz@brown.edu  #Email where notifications will be sent


#module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1

#for CHR in {1..22}; do
#bcftools view -m2 -M2 -i 'INFO/VT="SNP"' -Oz -o iinfo_biallelic_snps_chr${CHR}.vcf.gz /users/bkirz/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#done

#for CHR in {1..22}; do
#bcftools view -m2 -M2 -v snps -Oz -o vsnp_biallelic_snps_chr${CHR}.vcf.gz /users/bkirz/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#done


#for CHR in {1..22}; do
#python indel_pos.py vsnp_biallelic_snps_chr${CHR}.vcf.gz /users/bkirz/scratch/kirz_tutorials/vcf_processing/VSNP_INDELPOS/vsnps_chr${CHR}.txt
#done

#for CHR in {1..22}; do
#python indel_pos.py iinfo_biallelic_snps_chr${CHR}.vcf.gz /users/bkirz/scratch/kirz_tutorials/vcf_processing/IINFO_INDELPOS/iinfo_chr${CHR}.txt
#done


#for CHR in {1..22}; do
#python indel_diff.py VSNP_INDELPOS/vsnps_chr${CHR}.txt IINFO_INDELPOS/iinfo_chr${CHR}.txt /users/bkirz/scratch/kirz_tutorials/vcf_processing/INDEL_DIFF_RESULTS/indel_diff_chr${CHR}.txt
#done
