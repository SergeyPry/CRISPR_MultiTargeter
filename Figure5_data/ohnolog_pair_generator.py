import os

# first read all the lines from the text file

# create a filehandle
read_fh = open("ohnologues_pairs_zebrafish_E63.txt", 'r')

# create a dictionary
gene_pairs = {}

# process all lines in the file
for line in read_fh:
	# remove the newline
	line = line.strip('\n')
	
	# divide the line into two substrings
	(group1, group2) = line.split('\t')
	
	# divide each of the substrings into separate lists
	list1 = group1.split()
	list2 = group2.split()
	
	# make all possible pairs between the gene(s) in the list1 and the gene(s) in the list2
	# these pairs should consist of gene identifiers separated by a tab
	# store each pair as a key of a dictionary to make sure that each is unique
	
	for gene1 in list1:
		for gene2 in list2:
			# generate a pair from two genes in different groups
			pair = gene1 + '\t' + gene2
			
			# create a new entry in a dictionary
			if pair not in gene_pairs:
				gene_pairs[pair] = 1

read_fh.close()	

# once the dictionary is populated, write its contents to an output file

# open the file to write to  
write_fh = open("ohnologues_pairs_output.txt", 'w')

# code to output all the pairs to the file
for pair in gene_pairs:
	write_fh.write(pair)
	write_fh.write("\n")
