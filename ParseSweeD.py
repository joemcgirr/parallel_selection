# This script generates a data frame contatining results from SweeD,
# a software that identifies signitures of selective sweeps.
# It takes the output file from SweeD and organizes it for easy visualization 
# in R (see https://github.com/joemcgirr/parallel_selection/blob/master/parallel_selection.R).
# It also requires a file that orders chromosomes or scaffolds as they appear in 
# the vcf analyzed by SweeD 
# (see https://github.com/joemcgirr/parallel_selection/blob/master/example_data/Neoformans_chrom_key.txt)
#
#
# Example usage: 
# python3 ParseSweeD.py -species1 species2 species3
#
# Script will look for files named:
# SweeD_Report.<species1-3>_grid1000 and <species1-3>_chrom_key.txt




import argparse

parser = argparse.ArgumentParser(description='Parse SweeD Report to create results table')
parser.add_argument('species_names', metavar='names', nargs='+',type = str,
                    help='one or more species names matching prefix(s) for .vcf file(s). ex: species1 species2 species3')

args = parser.parse_args()

for species_list in vars(args).values():
	for species in species_list:
     
		report = 'SweeD_Report.'+species+'_grid1000'
		key = species+'_chrom_key.txt'
		
		def scaf_dict(key_file):
			infile = open(key_file,'r')
			iterable = iter(infile)
			keys = []
			values = []
			for line in iterable:
				entry = line.strip().split("\t")
				key = entry[0].rstrip()
				keys.append(key)
				values.append(entry[1])
			scaffold = dict(zip(keys,values))
			return scaffold
	
		lookup = scaf_dict(key)
		
		infile = open(report, 'r')
		iter_infile = iter(infile) 
		file = open(report+'_sweeps.txt', 'w')
		# count_windows = 0
		for line in iter_infile:
			if line.startswith("//"):
				sweed_key = lookup[line[2:].strip()]
				continue
			elif line[0].isdigit():
				line = sweed_key + '\t' + line
				file.write(line)
		file.close()
 
 # Note: will throw sweed key error if last line in chrom_key is empty