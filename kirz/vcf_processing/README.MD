# Filtering chromosomes to include only biallelic SNPS
___
### `biallelic.py`
Filters a gzipped vcf file into a new gzipped vcf file that only contains biallelic SNPs

```python
python biallelic.py vcf_input vcf_output
```

### `bcf_biallelic.sh`
shell script that does the same thing as `biallelic.py` but loads the `bcftools` module and uses it to filter using one or both of the following commands:

```bash
bcftools view -m2 -M2 -i VT="SNP" -Oz -o iinfo_CHR.vcf.gz source_CHR.vcf.gz
```

```bash
bcftools view -m2 -M2 -v snps -Oz -o vsnp_CHR.vcf.gz source_CHR.vcf.gz
```

### `compare_vcf.py`
returns whether two vcf files (such as one produced by `biallelic.py` and one by `bcftools` filtering) are identical to each other

```python
python compare_vcf my_biallelic.vcf.gz bcftools_biallelic.vcf.gz
```

### `indel_pos.py`
Takes in a gzipped vcf file (such as one produced by filtering out indels with `bcftools`)
and writes all indels present in the file to a new text file in the following tab-delineated column format:

`CHROM   |   POS`

```python
python indel_pos.py filtered_CHR.vcf.gz out.txt
```

### `indel_diff.py`
Takes in two text files that resulted from indel_pos.py and compares them to look for discrepancies.
If there are any discrepancies in the text files, the text file output records them in the following format:

`CHROM | VERSION (-v or -i) | POS`

```python
python indel_diff.py vsnp_CHR.txt iinfo_CHR.txt indel_diff_CHR.txt
```
