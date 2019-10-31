import argparse

parser = argparse.ArgumentParser(description='Parse SweeD Report to create results table')
parser.add_argument('species_names', metavar='names', nargs='+',type = str,
                    help='one or more species names matching prefix(s) for .vcf file(s). ex: species1 species2 species3')

# parser = argparse.ArgumentParser(description='Parse SweeD Report to create results table')
# parser.add_argument('-l', '--list', help='delimited list input', type=str)
# args = parser.parse_args()
# my_list = [str(item) for item in args.list.split(',')]

args = parser.parse_args()

for species_list in vars(args).values():
	for species in species_list:
		#print(species)
		#print(type(species))
     
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
 
 # might throw sweed key error if last line in chrom_key is empty