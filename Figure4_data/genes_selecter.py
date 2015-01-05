
read_fh = open("zebrafish_transcript_table.txt", 'r')

alt_genes  = {}

for line in read_fh: 
	# remove the newline
	line = line.strip('\n')
	
	# get the information on the gene and transcript ID
	(trID, geneID, seq) = line.split("\t")
	
	# record the number of transcripts are associated with any specific gene
	if geneID in alt_genes:
		alt_genes[geneID] += 1
	else:
		alt_genes[geneID] = 1

		
# make a new file handle for writing
outputFH = open('genes_with_alternative_transcripts.txt', 'w')
		
# output the genes with multiple transcript isoforms		
for gene in alt_genes:
	if alt_genes[gene] >= 2:
		outputFH.write(gene)
		outputFH.write('\n')
		

outputFH.close()