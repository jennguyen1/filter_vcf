# Filters VCF Files

Extracts SNPs from VCF filese and combines into an extracted list. Can be useful for extracting a small subset of SNPs for analysis such as computing a genetic risk score.

* extract_snps.py extracts a list of SNPs from VCF files
  * wrapper for <tabix -p vcf -f {in_file} {c}:{p}-{p}>
* combine_vcf.py combines VCF files

For more information on running the program, type the command
```
python extract_snps.py -h
python combine_vcf.py -h 
```
