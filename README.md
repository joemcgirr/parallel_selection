# Parallel Selection
> This repository contains scripts used to identify genomic regions under selection in multiple lineages.
>
> Analyses begin with a `.vcf` file generated with samtools mpileup

## parallel_selection.ipynb jupyter notebook generates scripts to:
1. Quality filter `.vcf` to include only biallelic SNPs
2. Run selective sweep scan with SweeD
3. Parse SweeD output
4. Identify gene features within sweeping regions 
5. Identify gene features sweeping in multiple populations and plot results


# Commands
## 1. Quality filter `.vcf` to include only biallelic SNPs
> Requires samtools and vcftools
```
bcftools call --threads 4 -m -Ov -o species1_snps.vcf species1_raw_variants.vcf
vcftools --vcf species1_snps.vcf --minQ 20 --maf 0.05 --max-missing 0.9 --remove-indels --min-alleles 2 --max-alleles 2 --recode --out species1_filtered_snps
sed 's/,\*//g' species1_filtered_snps.recode.vcf > species1_filtered_snps_sweed.vcf
```

## 2. deduplicate .bam files with picard.jar
```
java -Xmx10g -jar picard.jar MarkDuplicates INPUT=sample.sort.bam OUTPUT=sample.sort.dedup.bam METRICS_FILE=sample.metrics.txt MAX_FILE_HANDLES=1000
samtools index sample.sort.dedup.bam
```
## 3. call snps with gatk 3.8
```
gatk -T HaplotypeCaller -ERC GVCF -drf DuplicateRead -R reference.fasta -I sample.sort.dedup.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -nct 4 -o sample_raw_variants.g.vcf
gatk -T GenotypeGVCFs -R reference.fasta --variant sample1_raw_variants.g.vcf --variant sample2_raw_variants.g.vcf -o merged_raw_variants.vcf
vcftools --vcf merged_raw_variants.vcf --maf 0.05 --max-missing 0.9 --remove-indels --recode --out filtered_snps.vcf
```
## 4. calculate Fst, Tajima's D, and pi with vcftools
```
vcftools --vcf filtered_snps.vcf --keep populations_1_and_2.txt --out population_1_vs_2.weir.fst --weir-fst-pop population_1.txt --weir-fst-pop population_2.txt
vcftools --vcf filtered_snps.vcf --TajimaD 20000 --out taj_d_20kb_windows.txt 
vcftools --vcf filtered_snps.vcf --window-pi 20000 --out pi_20kb_windows.txt
```
> bash commands to remove negative values from Fst output and calculate genome-wide mean Fst
```
sed \'s/-[0-9].*/0/g\' population_1_vs_2.weir.fst | sed \'s/-nan/0/g\' > population_1_vs_2.weir.fst
awk -F\'\\t\' \'{ sum += $3 } END { print sum / NR }\' population_1_vs_2.weir.fst > genome_wide_avg.txt
```
## 5. calculate Dxy with simon martin scripts
> see [simonhmartin/genomics_general](https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing) for `parseVCF.py` and `popgenWindows.py`
```
python parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | bgzip > output.geno.gz
python popgenWindows.py -w 10000 -m 10 -g output.geno.gz -o popgen_stats.csv -f phased -T 4 -p population_1 sample_1 sample_2 -p population_2.txt sample_3 sample_4
```
## 6. find hard sweeps with SweeD
> see [SweeD documentation](https://cme.h-its.org/exelixis/resource/download/software/sweed3.0_manual.pdf) for details
>
> -eN flag used to specify demography specific to San Salvador pupfish system (population decrease 10kya)
```
SweeD -input filtered_snps.vcf -s 100 -name sweeps.txt -eN .005 .01 -strictPolymorphic -folded -grid 1000 -osfs sweeps_sfs.txt
```
## 7. create RaxML plylogeny
> see [RaxML documentation](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf) for details
>
> first create biallelic `.vcf` and convert to `.phy` with `vcf2phylip.py` from [edgardomortiz](https://github.com/edgardomortiz/vcf2phylip)
```
vcftools --vcf filtered_snps.vcf --min-alleles 2 --man-alleles 2 --recode --out biallelic_snps.vcf

python vcf2phylip.py --input biallelic_snps.vcf

standard-RAxML-master/raxmlHPC -x 76345 -f a -m GTRGAMMA -o outgroup_sample_name -p 28393 -s biallelic_snps.phy -# 10 -n output.name
```
## 8. GWAS with GEMMA
> see [GEMMA documentation](https://www.xzlab.org/software/GEMMAmanual.pdf) for details
>
> first convert biallelic `.vcf` to [plink](http://zzz.bwh.harvard.edu/plink/) format
```
vcftools --vcf biallelic_snps.vcf --plink --out biallelic_snps.plink
plink.exe --file biallelic_snps.plink --make-bed --out biallelic_snps.plink
gwas/gemma-0.98.1-linux-static -bfile biallelic_snps.plink -w 50000000 -s 100000000 -n 5 -rpace 10000 -wpace 100000 -bslmm 1 -o run_1.output
```




