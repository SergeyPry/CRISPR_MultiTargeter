######################
# Exons
######################
import re

# first it is necessary to read the data from the exon sequence file
exon_fh = open('exons_new_BioMart.txt', 'r')

# a giant string to hold the whole file and process
exon_data = ""
exon_data = exon_fh.read()
exon_records = []
exon_records = exon_data.split('>')
exon_data = ""

species = ""
species = "Danio rerio"

# close the filehandle for the file containing the exons 
exon_fh.close()

# create a dictionary to store exons
exons = {}

for rec in exon_records:
	lines = []
	lines = rec.split('\n')
	
	if len(lines[0]) > 10 and len(lines) > 1:
		info = []
		info = lines[0].split('|')

		ExonID = info[2]
		EnsGeneID = info[0]
		EnsGeneID = EnsGeneID.strip('>') # in case the Ensembl ID still has the ">" character
		
		ChrStart = int(info[3])
		ChrEnd = int(info[4])
		strand = int(info[5])
		
		exons[ExonID] = {}
		exons[ExonID]['EnsGeneID'] = EnsGeneID
		exons[ExonID]['ChrStart'] = ChrStart
		exons[ExonID]['ChrEnd'] = ChrEnd
		exons[ExonID]['strand'] = strand
		
		sequence = ""
		
		for line in lines[1:]:
			sequence += line.strip('\n')
		
		exons[ExonID]['sequence'] = sequence
		# an example of the header, for reference purposes
		# >ENSDARG00000092340|ENSDART00000140911|ENSDARE00001018997|42055|42228|-1|32894|44457

# create a dictionary for genes and store the relevant data in it
# open a file containing the basic information about zebrafish genes
gene_info_fh = open('genes-symbols_BioMart.txt')

	
# add the exon information to the dictionary above (preferably sorted according to the coordinates and strand)
exons_genes = {}

for exon in exons.keys():
	geneID = exons[exon]['EnsGeneID']
	
	if geneID in exons_genes:
		# now we need to store both the start and end coordinates of the exons to facilitate subsequent sorting
		if exons[exon]['strand'] == 1:
			exons_genes[geneID]['forw_exons'][exon] = exons[exon]['ChrStart']
		elif exons[exon]['strand'] == -1:
			exons_genes[geneID]['rev_exons'][exon] = exons[exon]['ChrEnd']

	else:
		exons_genes[geneID] = {}
		exons_genes[geneID]['forw_exons'] = {}
		exons_genes[geneID]['rev_exons'] = {}
		# Add exons to the data structure
		if exons[exon]['strand'] == 1:
			exons_genes[geneID]['forw_exons'][exon] = exons[exon]['ChrStart']
		elif exons[exon]['strand'] == -1:
			exons_genes[geneID]['rev_exons'][exon] = exons[exon]['ChrEnd']	
			


# create a file to append the data for genes and their exons
genes_table_fh = open('zebrafish_genes_table.txt', 'w')
gene_exon_wr_fh = open('zebrafish_gene-exon_table.txt', 'w')
failed_exons_wr = open('failed_exons.txt', 'w')


# define counts for successful and failed matches of exons in genes
match_count = 0
fail_count = 0
gene_count = 0
failed_exons = {}

# combine the exons constituting each gene into a single sequence, also dealing with cases of overlapping exons
for line in gene_info_fh:
	# obtain all the necessary data from the file line
	
	(geneID,symbol, strand, start, end)= line.split("\t")
	strand = int(strand)
	start = int(start)
	end = int(end)
	
	# collect exon IDs in the order of their coordinates to later use when merging all exons
	sorted_exonIDs = []
	
	# do the sorting of exons, the basic algorithm being:
	# sorted(mydict, key=lambda key: mydict[key])
	
	# adjust the order of exons according to the strand
	if strand == 1:
		
		# the first coordinate is the ChrStart
		sorted_exonIDs = sorted(exons_genes[geneID]['forw_exons'], key = lambda key: exons_genes[geneID]['forw_exons'][key])
	elif strand == -1:
	
		# the second coordinate is the ChrEnd, which in this case the starting coordinate
		# in addition the order of reversed elements has to be reversed
		sorted_exonIDs = sorted(exons_genes[geneID]['rev_exons'], key = lambda key: exons_genes[geneID]['rev_exons'][key])
		sorted_exonIDs.reverse()
		
		
	
	# Now we need to merge all exons using an algorithm which compares coordinates of neighboring exons
	
	# 1. Initialize the merged sequence
	
	merged_seq = ""
	
	# 2. Go over all of the sorted exons using an index
	# for the position inside a list of exons
	
	# to keep track of where we are in the gene we will have the last_coord variable
	last_coord = 0
	
	for j, val in enumerate(sorted_exonIDs):
		# The algorithm will be a little different depending on the strand
		# of the gene. Therefore, the first 'if' branch will take care of that
		
		if strand == 1:
			# exon merging code for the +1 strand
			
			# this is the first exon
			if j == 0:
				merged_seq += exons[val]['sequence']
				last_coord = exons[val]['ChrEnd']
			
			# the next exon does not overlap with the previous
			elif last_coord < exons[val]['ChrStart']:
				merged_seq += exons[val]['sequence']
				last_coord = exons[val]['ChrEnd']
			
			# the next exon overlaps with the previous one
			elif last_coord > exons[val]['ChrStart']:
				
				# this exon may overlap with the previous one
				# or be within it
				# case 1: it overlaps
				
				if last_coord < exons[val]['ChrEnd']:
					substr_coord = 0
					substr_coord = int(last_coord) - int(exons[val]['ChrStart']) + 1
					merged_seq += exons[val]['sequence'][substr_coord:]
					last_coord = exons[val]['ChrEnd']
				# case 2: it is within the previous exon - no need to do anything, which means it will be skipped
				
			
		elif strand == -1:
			# exon merging code for the -1 strand
			
			# this is the first exon
			if j == 0:
				merged_seq += exons[val]['sequence']
				last_coord = exons[val]['ChrStart']
			
			# the next exon does not overlap with the previous
			
			elif last_coord > exons[val]['ChrEnd']:
				merged_seq += exons[val]['sequence']
				last_coord = exons[val]['ChrStart']
			
			# the next exon overlaps with the previous one
			elif last_coord < exons[val]['ChrEnd']:
				
				# this exon may overlap with the previous one
				# or be within it
				# case 1: it overlaps
				
				if last_coord > exons[val]['ChrStart']:
					substr_coord = 0
					substr_coord = int(exons[val]['ChrEnd']) - int(last_coord) + 1
					merged_seq += exons[val]['sequence'][substr_coord:]
					last_coord = exons[val]['ChrStart']
				# case 2: it is within the previous exon - no need to do anything, which means it will be skipped
	
	
	
	# Now we need to write the results to the Genes table
	genes_table_fh.write(geneID)
	genes_table_fh.write('\t')
	
	if symbol == "":
		genes_table_fh.write('NULL')
	else:
		genes_table_fh.write(symbol)
	genes_table_fh.write('\t')
	genes_table_fh.write(species)
	genes_table_fh.write('\t')
	
	genes_table_fh.write(merged_seq)
	genes_table_fh.write('\n')
	
	# The next step is to check each exon against the merged sequence and find out its position
	# as it has been done in case of genes

	
	genes_exons = {}

	if strand == 1:
		
		for exon in exons_genes[geneID]['forw_exons']:
			match = re.search(exons[exon]['sequence'], merged_seq)
			
			if match != None:
				coord = match.span()
				# store the results of exon sequence match
				genes_exons[exon] = []
				genes_exons[exon].append(coord[0])
				genes_exons[exon].append(coord[1])
				match_count += 1
			else:
				fail_count += 1
				failed_exons_wr.write(exon)
				failed_exons_wr.write('\t')
				failed_exons_wr.write(str(exons[exon]['ChrStart']))
				failed_exons_wr.write('\t')
				failed_exons_wr.write(str(exons[exon]['ChrEnd']))
				failed_exons_wr.write('\t')
				failed_exons_wr.write(geneID)
				failed_exons_wr.write('\n')
	
	elif strand == -1:
		
		for exon in exons_genes[geneID]['rev_exons']:
			match = re.search(exons[exon]['sequence'], merged_seq)
			
			if match != None:
				coord = match.span()
				# store the results of exon sequence match
				genes_exons[exon] = []
				genes_exons[exon].append(coord[0])
				genes_exons[exon].append(coord[1])
			
			if match != None:
				match_count += 1
			else:
				fail_count += 1
	

	# write the results of the match
	for exon in genes_exons.keys():
		gene_exon_wr_fh.write(geneID)
		gene_exon_wr_fh.write('\t')
		gene_exon_wr_fh.write(exon)
		gene_exon_wr_fh.write('\t')
		gene_exon_wr_fh.write(str(genes_exons[exon][0]))
		gene_exon_wr_fh.write('\t')
		gene_exon_wr_fh.write(str(genes_exons[exon][1]))
		gene_exon_wr_fh.write('\n')
			
print("The number of exons is the following: ")
print(len(exons))
print "The number of exons which match their genes is %d: " %match_count
print "The number of exons which do NOT match their genes is %d: " %fail_count

genes_table_fh.close()
gene_exon_wr_fh.close()