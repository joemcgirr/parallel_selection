# Parallel Selection
> This repository contains scripts used to identify genomic regions under selection in multiple lineages with SweeD.
>
> see [SweeD documentation](https://cme.h-its.org/exelixis/resource/download/software/sweed3.0_manual.pdf) for details
>
> Analyses begin with a `.vcf` file generated with samtools mpileup

## parallel_selection.ipynb jupyter notebook generates scripts to:
1. Quality filter `.vcf` to include only biallelic SNPs
2. Run selective sweep scan with SweeD
3. Parse SweeD output
4. Identify gene features sweeping in multiple populations and plot results 



# Commands
## 1. Quality filter `.vcf` to include only biallelic SNPs
> Requires samtools and vcftools
```
bcftools call --threads 4 -m -Ov -o species1_snps.vcf species1_raw_variants.vcf
vcftools --vcf species1_snps.vcf --minQ 20 --maf 0.05 --max-missing 0.9 --remove-indels --min-alleles 2 --max-alleles 2 --recode --out species1_filtered_snps
sed 's/,\*//g' species1_filtered_snps.recode.vcf > species1_filtered_snps_sweed.vcf
```

## 2. Run selective sweep scan with SweeD
```
./SweeD-P -input species1_filtered_snps_sweed.vcf -name species1_grid1000 -strictPolymorphic -folded -grid 1000 -osfs species1_sfs_grid1000 -threads 4

```
## 3. Parse SweeD output
```
sed '/^#/ d' species1_filtered_snps_sweed.vcf | awk '{print $1}' | cat -n | sort -k2 | uniq -f1 -d | sort -n | cut -f2- | nl > species1_chrom_key.txt
```
> After creating `species1_chrom_key.txt` run `ParseSweeD.py` to create a dataframe with results
```
python3 ParseSweeD.py -species1 species2 species3
```

## 4. Identify gene features sweeping in multiple populations and plot results
> specify directories and names for gene feature tables in parallel_selection.R


