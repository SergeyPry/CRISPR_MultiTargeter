#!/usr/bin/env python

##########################
# Preliminaries
##########################

# loading all the necessary modules
import cgi, cgitb
import os
import sys
import re
import string
import os.path
import random
import string
import sqlite3 as db

# Biopython
import Bio
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align.Generic import Alignment
from Bio.Seq import Seq

# to enable exception handling
cgitb.enable()

################
# functions
################

dna = set("ATGC")
iupac_dna = set("ATGCWSMKRYBDHVN")

iupac_dict = {}
iupac_dict["A"] = "A"
iupac_dict["T"] = "T"
iupac_dict["C"] = "C"
iupac_dict["G"] = "G"
iupac_dict["W"] = "[AT]"
iupac_dict["S"] = "[GC]"
iupac_dict["M"] = "[AC]"
iupac_dict["K"] = "[GT]"
iupac_dict["R"] = "[AG]"
iupac_dict["Y"] = "[CT]"
iupac_dict["B"] = "[CGT]"
iupac_dict["D"] = "[AGT]"
iupac_dict["H"] = "[ACT]"
iupac_dict["V"] = "[ACG]"
iupac_dict["N"] = "[ATGC]"

def validate(seq, alphabet=dna):
    # "Checks that a sequence only contains values from an alphabet"
    leftover = set(seq.upper()) - alphabet
    return not leftover
	
def validate_pam(pam_seq, alphabet=iupac_dna):
    leftover = set(pam_seq.upper()) - alphabet
    return not leftover

def simple_regex_generate(dinuc, length, pam_seq, oriented, mismatch):
	# initialize the regular expression string 
	regex = ""
	
	# mismatch is allowed
	if mismatch == "Yes":
		# regex string creation based on the regular expression language
		# add the 5'-terminal dinucleotide and the rest of the target
		# according to the target length
		if dinuc == "GG":
			regex += "("
			regex += dinuc
			regex += "[GATCX]{%d}" % 6
			length_right = int(length) -8
			regex += "[GATC]{%d})" % length_right
		elif dinuc == "GN":
			regex += "("
			regex += dinuc
			regex += "[GATCX]{%d}" % 6
			length_right = int(length) -8
			regex += "[GATC]{%d})" % length_right
		else:
			total_length = int(length)
			regex += "([GATCX]{%d}" % 8
			length_right = total_length -8
			regex += "[GATC]{%d})" % length_right
			
	elif mismatch == "No":
		if dinuc == "GG":
			regex += "("
			regex += dinuc
			length_right = int(length) -2
			regex += "[GATC]{%d})" % length_right
		elif dinuc == "GN":
			regex += "("
			regex += dinuc
			length_right = int(length) -2
			regex += "[GATC]{%d})" % length_right
		else:
			total_length = int(length)
			regex += "[GATC]{%d})" % total_length
	
	if pam_seq == "NGG":
	
		if oriented == 3:
			regex += "[GATCX]GG"
		elif oriented == 5:
			badrequest_oriented(pam_seq)
		
	else:
		# assume that this alternative PAM sequence has been previously validated
		# and includes only valid IUPAC DNA characters
		# iupac_dna = set("ATGCWSMKRYBDHVN")
		
		if oriented == 3:
			for c in pam_seq:
				# use the previously-created dictionary to generate a full regular expression
				# for an alternative PAM sequence
				code_letter = c.upper()
				regex += iupac_dict[code_letter]
		elif oriented == 5:
			pam_iupac = ""
			
			for c in pam_seq:
				# use the previously-created dictionary to generate a full regular expression
				# for an alternative PAM sequence
				code_letter = c.upper()
				pam_iupac += iupac_dict[code_letter]
			
			regex = pam_iupac + regex
		
	return(regex)

	
def rc(dna):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = dna.translate(complements)[::-1]
	return rcseq

def compute_best_sgRNA_target(targ_seq, targ_id, sequences_targets):
	
	# initialize the sequence to return
	best_seq = ''
	
	# if there is no mismatch, return the actual target sequence
	if 'X' not in targ_seq:
		best_seq = targ_seq
	else:
	# compute which nucleotide should replace the 'X'
	
		# get the location of 'X' in the target sequence	
		loc = targ_seq.index('X')
		
		# extract which nucleotides are present in all of the sequences in the alignment
		nucs = []
		
		for id in sequences_targets.keys():
			if id != "consensus":
				seq = str(sequences_targets[id]["aligned_matches"][targ_id]["target_seq"])
				nucs.append(seq[loc])
		
		# use some simple decision rule to optimize the choice of the nucleotide to insert into
		# the final targeting sequence, which will be used for the design of an sgRNA
		
		# For example:
		# perfect match = 4
		# R with Y or Y with R = 2
		# R with R or Y with Y = 0
		# think of the biological rationale for this or an empirical reasone to support such a scoring scheme
		
		
		base_scores = {}
		
		for base in 'ATCG':
			# initialize the dictionary
			base_scores[base] = 0
			
			# iterate over all of the bases in different sequences and compute the total score for each nucleotide
			for nuc in nucs:
				if base == nuc:
					base_scores[base] += 4
				elif (base in 'AG') and (nuc in 'AG'):
					base_scores[base] += 2
				elif (base in 'TC') and (nuc in 'TC'):
					base_scores[base] += 2
				elif (base in 'AG') and (nuc in 'TC'):
					base_scores[base] += 0
				elif (base in 'TC') and (nuc in 'AG'):
					base_scores[base] += 0
		
		# compute which base will get a higher score and then sort them by these scores	
		best = ''
		max_score = 0
		
		# iterate over the best_scores dictionary 
		for key in base_scores.keys():
			if base_scores[key] >= max_score:
				best = str(key)
				max_score = base_scores[key]
		
		best_seq = targ_seq.replace('X', best)
	
	return best_seq
	

######################################
# gRNA target search - site definition
######################################

# variables assigned manually for the purposes of systematic common target analysis
# in zebrafish ohnologs

# first two nucleotides - of some importance in current sgRNA molecules
dinuc = "NN"

# overall length
length = 20

# orientation of the PAM sequence
oriented = 3
# mismatches between target sites
mismatch = "Yes"

# define the PAM sequence  
pam_seq = "NGG"

# define species
species = "Danio rerio"

# retrieve the molecule type
mol_type = "ensgenes"

############################################
# Read the input file and store its contents
# for subsequent analysis
############################################

read_fh = open("ohnologues_pairs_output.txt", 'r')

# create a dictionary
gene_pairs = {}

# process all lines in the file
for line in read_fh:
	# remove the newline
	line = line.strip('\n')
	
	# divide the line into two substrings
	(gene1, gene2) = line.split('\t')

	# make all possible pairs between the gene(s) in the list1 and the gene(s) in the list2
	# these pairs should consist of gene identifiers separated by a tab
	# store each pair as a key of a dictionary to make sure that each is unique

	# generate a pair from two genes in different groups
	pair = gene1 + '\t' + gene2
			
	# create a new entry in a dictionary
	if pair not in gene_pairs:
		gene_pairs[pair] = 1

read_fh.close()	

##########################################################
# RETRIEVAL OF DATA FROM THE DATABASE
##########################################################


# 2. Prepare data structures for the subsequent steps in the process, such as a basic sequences dictionary with transcript sequences,
# but also exon structures of each transcript.

# a dictionary containing information about the exons of 
seqs_exons = {}

# initialize a dictionary for storing the results of the run
results = {}

# 1. establish a database connection
conn = db.connect('genomes.db')
cur = conn.cursor()

# iterate the whole algorithm over all available gene pairs
for pair in gene_pairs:

	# generate a list with gene identifiers 
	geneids = pair.split('\t')
	
	# another data structure will only contain the sequence information for input genes or transcripts
	genes_transcripts = {}
	
	# use the resulting list in the analysis 
	for gene_id in geneids:

		# copy the initial value of an identifier
		orig_id = gene_id
		
		# make the identifier uppercase, but then use case-insensitive comparison
		gene_id.upper()
		
		# first check if the gene identifier has a format of an Ensembl ID
		if mol_type == "ensgenes":
			# try to select the data from Genes table on the assumption that the gene identifier is of Ensembl ID format
			cur.execute("SELECT * from Genes WHERE geneid = ? COLLATE NOCASE", (gene_id,))
			row = cur.fetchone()
			
			if row is None:
				break
			else:		
				# retrieve the data from the Genes table
				geneid = row[0]
				symbol = row[1]
				sequence = row[2]
				species = row[3]
				
				# remove the newline from the species variable
				species = species.rstrip()
				sequence = sequence.rstrip()
				
				# store the sequence for this identifier
				genes_transcripts[orig_id] = str(sequence)
				
				# now it is necessary to obtain the exon information for this gene
				seqs_exons[orig_id] = {}
				
				# extract exon sequences from the database
				for row in cur.execute("SELECT exonID, sequence FROM Exons WHERE geneid = ?", (geneid,)):
					# data from the current row
					exonID = row[0]
					exon_seq = row[1]
					
					# store them in the dictionary
					seqs_exons[orig_id][exonID] = str(exon_seq)		
					
	# create a dummy file which will be used as an input for the ClustalW2 program
	tmp_path = "D:\\Python_bioinformatics_projects\\CRISPR_MultiTargeter\\cgi-bin\\tmp\\"

	###########################################################
	# Random folder generation for storing each user's results
	###########################################################

	rand_folder = str(''.join(random.sample(string.letters*8,8)))
	new_folder = tmp_path + rand_folder + '\\'
	os.mkdir(new_folder)

	input_file = new_folder  + 'alignment_input.fas'

	align_input_file = open(input_file, 'w')

	# writing sequences to the file  for alignment the proper Biopython way
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio import SeqIO


	# Use these commands later on in the program to delete the directory:
	# os.remove(input_file)
	# os.rmdir(new_folder)

	# to keep the subsequent code consistent with other scripts,
	# the "transcripts" dictionary to the "sequences" dictionary
	sequences = genes_transcripts

	for key in sequences.keys():
		
		id_line = ">" + key
		
		align_input_file.write(id_line)
		align_input_file.write("\n")
		
		i = 0
		incr = 80
		sequence = sequences[key]
		
		while i <= len(sequence):
			next = i + 80
			
			if next > len(sequence):
				align_input_file.write(sequence[i:len(sequence)])
				align_input_file.write("\n")
			else:
				align_input_file.write(sequence[i:next])
				align_input_file.write("\n")
				
			i += 80
		

	align_input_file.close()


	##########################################################################
	# 2. Align using ClustalW2 and read the alignment into an object
	##########################################################################

	# Now that the sequences are read into the file we need to have some code
	# which will run the alignment program and read the resulting data into a suitable
	# object for downstream processing

	# setup of the alignment software and running ClustalW2

	# Comment:
	# 
	# There is some discussion about whether multiple sequence alignment tool is the best one for aligning
	# two sequences. To establish if I should use one of the pairwise alignment tools, I tried running
	# Smith-Waterman and Needleman-Wunsch tools and found that they produce very gap-rich alignments
	# on somewhat dissimilar sequences. ClustalW2 produced much fewer gaps and was chosen as the primary alignment tool.

	from Bio.Align.Applications import ClustalwCommandline
	clustalw_exe = r"C:\\Program Files\\ClustalW2\\clustalw2.exe"

	# generate an output file in the same folder
	output_file = new_folder  + 'alignment_output.aln'
	dnd_file = new_folder  + 'alignment_input.dnd'

	clustalw_cline = ClustalwCommandline(clustalw_exe, infile= input_file, outfile= output_file)
	assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
	stdout, stderr = clustalw_cline()

	# read the resulting alignment from a file generated by the preceding code
	from Bio import AlignIO

	align = AlignIO.read(output_file, "clustal")

	##########################################################################
	# 3. Calculate the consensus sequence from the alignment and create a data 
	# structure to hold the aligned sequences as well as the consensus 
	# sequence expressed in the DNA and IUPAC ambiguous alphabet for DNA and 
	# gap '-' characters. It must be possible to match a regular expression
	# against the consensus and derive the corresponding matches from the 
	# input sequences.
	##########################################################################

	# using Biopython to calculate the consensus of the
	# alignment created in the previous step
	from Bio.Align import AlignInfo
	summary_align = AlignInfo.SummaryInfo(align)

	# consensus is calculated using the convention that it produces
	# the character 'X' if less than 95 % of residues at that position have the same character

	consensus = summary_align.gap_consensus(threshold=0.95)
	rev_consensus = consensus.reverse_complement()
	
	
	
	# consensus = str(consensus).replace('X', ' ')
	# rev_consensus = str(rev_consensus).replace('X', ' ')
	# consensus = re.sub('[ACTG]', '*', consensus)
	# rev_consensus = re.sub('[ACTG]', '*', rev_consensus)

	# Now it is necessary to set up a data structure which will store both the sequences in
	# the alignment as well as their consensus
	# create separate data structures for both the sense strand and its reverse_complement

	# The following data structures are dictionaries to store identifiers of sequences,
	# their actual sequences and the overall consensus sequence for the alignment
	# this is done for both the main sequence and its reverse complement

	align_cons = {}
	rev_align_cons = {}

	for record in align:
		
		# store sequences in dictionaries for subsequent processing
		align_cons[record.id] = str(record.seq)	
		rev_align_cons[record.id] = str(record.seq.reverse_complement())

	# add consensus sequences
	align_cons["consensus"] = str(consensus)
	rev_align_cons["consensus"] = str(rev_consensus)
	
	##########################################################################
	# 4. Create a regular expression consistent with all the recent studies 
	# and frequent enough to make its identification in similar regions 
	# worthwhile.
	##########################################################################

	# The latest studies of CRISPR targeting specificity indicate that 
	# it is not essential for the 5'-terminal GG sequence to be present in the
	# targeted genomic DNA (Hwang et al., 2013). This result expands the sequence targeting range
	# to 1 CRISPR site every 8 nucleotides. The authors tested two versions of their
	# CRISPRs without GG sequence present in the genome: 18 and 20 nt. In this tool,
	# we will use the average of these values, i.e. 19. This length of the CRISPR target
	# was also used in another study (Hsu et al., 2013), where the targeting specificity of
	# CRISPRs was examined in detail. Their results can be summarized as follows:
	# the 7 5'-most sequences are not very sensitive to mismatches except for these cases:
	#  rA:dC - position 5
	#  rU:dC - position 5
	#  rC:dC - position 5
	#  rU:dT - position 6

	# The regular expression to search for CRISPR target sites can thus be defined as follows
	# Additional details of the matches and mismatch rules between similar but not identical 
	# sequences will be checked later when filtering potential matches

	#  step 1: Construct a regular expression out of the available input
	#  components and compile it.

	sgRNA_regex = simple_regex_generate(dinuc, length, pam_seq, oriented, mismatch)
	# lookahead functionality to enable overlapping matches
	final_regex = "(?=%s)" % sgRNA_regex

	regex = re.compile(final_regex, re.IGNORECASE)

	#############################################################################
	# 5. Run the regular expression on the consensus sequence, record the matches and locations of matches.
	# Using locations, obtain the corresponding sequences from the aligned sequences
	#############################################################################	

	# after getting a match it will be necessary to count the number of Xs and reject if there are more
	# than one. Also, one needs to check for gaps in the aligned sequences

	# find all of the matches for the regular expression shown above

	#######################################
	# Finding targets in the forward strand
	#######################################

	match_count = 1 			# a variable to give numeric identifiers
	sequences_matches_forward = {} # a dictionary with the following structure:

	# At the top, there will be sequence identifiers: sequence IDs and "consensus", which will have several pieces 
	# of information in a data structure as their values.
	# More specifically, "align_seq" = sequence of from the alignment object
	# 					 "original_seq" = the sequence from the original data structure containing sequences
	#					 "aligned_matches" = matches of the sgRNA target regex extracted directly from the aligned sequence
	#					 "orig_matches" = matches of the sgRNA target regex extracted indirectly from the original sequence.
	#					 This last part of this data structure is required to be able to identify unique sequences by removing
	#					 them from the list of all possible matches and then checking all the other matches using some alignment
	#					 or fuzznuc-based approaches

	# initialize a dictionary for the consensus sequence 
	sequences_matches_forward["consensus"] = {}
	sequences_matches_forward["consensus"]["align_seq"] = str(align_cons["consensus"])
	sequences_matches_forward["consensus"]["aligned_matches"] = {}

	# initialize a dictionary for all individual aligned sequence
	for id in align_cons.keys():
		if id != "consensus":
			sequences_matches_forward[id] = {}
			sequences_matches_forward[id]["align_seq"] = align_cons[id]
		
			sequences_matches_forward[id]["original_seq"] = sequences[id]
			sequences_matches_forward[id]["aligned_matches"] = {}
			sequences_matches_forward[id]["orig_matches"] = {}

	# run regular expression matching of the consensus sequence of the alignment
	# The aim of the following code is to locate sgRNA targets common among all of the aligned sequences

	cons = str(consensus)

	for match in regex.finditer(align_cons["consensus"]):

		# The code to do a more detailed processing of the matches
		if match.group(1).count('X') <= 1:	# to exclude the sequences which contain more than 1 non-consensus site
						
	# The targets in the original sequences will have to be verified not to span exon boundaries. Here is the algorithm
	# which will be implemented in the following.
			
	# 1. Do the alignment with the transcript sequences retrieved by the user.
	#
	# 2. Provisionally store them in the dictionary for guide RNA targets.
	# 
	# 3. Go to the original sequences and retrieve both the target sequences themselves and their real coordinates.
	# 
	# 4. Run these target coordinates against exon coordinates inside the transcript.
	#    Target coordinates for the target matches in the reverse DNA strand need to be converted to real coordinates by
	#    subtracting them from the total length of sequence and reversing start and end coordinates.
	# 
	#    The rules for testing if a target is completely inside an exon are the following:
	#
	#    - (start of target is > than the start of exon) AND (end of target is > than end of exon)
	#    - (start of target is < than start of exon) AND (end of target is > than start of exon)
	#    - An alternative is to use the test of string presense in an exon sequence in a particular exon, but this may be
	#    computationally expensive
			
			# populate the nested dictionaries with the data
			start = match.start(1)
			end = match.end(1)

			for record in align:
			#initialization for all of the aligned sequences
				sequences_matches_forward[record.id]["aligned_matches"][match_count] = {}
				sequences_matches_forward[record.id]["orig_matches"][match_count] = {}
			
			# define a boolean variable, which is True if the target does not span two exons in any of the sequences
			cons_inside = True
			
			# Next, iterate over all of the aligned sequences
			for record in align:
			
				# The coordinates for matches in the original sequences are only different by the number of gaps in a part of that sequence
				# just before the match. Therefore, it should be possible to convert the data between the two parts of the dictionary by a simple 
				# subtraction
			
				# collect the essential information for the current match in the current sequence
				seq = str(record.seq)
				gaps = int(seq[:match.start()].count('-'))
				id = record.id
				target_start = start - gaps
				target_end = end - gaps       # PAM sequence has to be inside an exon as well
				
				# define a boolean variable to check if the target is inside a particular exon
				# it is defined as True initially because it is a priori more likely not to overlap two exons
				# once one of the conditions is satisfied, it becomes False and subsequent code execution 
				# occurs differently
				inside = False
				
				# obtain the sequence of the target and check that it is inside a particular exon
				target_seq = seq[start: end]
				
				# check for the molecule type selected by the user
				# the choice of the "refseq" option would prevent execution of the following code
				# and the cons_inside variable will remain True
				
				if mol_type != "refseq":
					# check if this target does not overlap two exons			
					for exonID in seqs_exons[id].keys():
						
						if target_seq in seqs_exons[id][exonID]:
							inside = True
					
					# Update the variable for the target sequence to be completely inside a particular exon to False if there
					# is a sequence which does not satisfy this
					if inside: 
						pass
					else:
						cons_inside = False
			
			
			# Test if there were any sequences, where the target overlapped two exons
			if cons_inside:
				# Consensus sequence initialization and storing the data
				# initialize the dictionary for the consensus and each of the matches
				sequences_matches_forward["consensus"]["aligned_matches"][match_count] = {}

				sequences_matches_forward["consensus"]["aligned_matches"][match_count]["target_seq"] = match.group(1)
				sequences_matches_forward["consensus"]["aligned_matches"][match_count]["coord_target"] = [start, end]
			
				if oriented == 5:
					sequences_matches_forward["consensus"]["aligned_matches"][match_count]["coord_match"] = [start-len(pam_seq), end]
				elif oriented == 3:
					sequences_matches_forward["consensus"]["aligned_matches"][match_count]["coord_match"] = [start, end + len(pam_seq)]


				# Next, iterate over all of the aligned sequences
				for record in align:		
					# matches in the original sequences:
					orig_seq = str(sequences[record.id])
					seq = str(record.seq)

					gaps = int(seq[:match.start()].count('-'))

					sequences_matches_forward[record.id]["orig_matches"][match_count]["target_seq"] = orig_seq[start - gaps: end - gaps]
					sequences_matches_forward[record.id]["orig_matches"][match_count]["coord_target"] = [start - gaps, end - gaps]

					if oriented == 5:
						sequences_matches_forward[record.id]["orig_matches"][match_count]["coord_match"] = [start - gaps - len(pam_seq), end - gaps]
					elif oriented == 3:
						sequences_matches_forward[record.id]["orig_matches"][match_count]["coord_match"] = [start - gaps, end - gaps + len(pam_seq)]


					# matches in the aligned sequences:
					sequences_matches_forward[record.id]["aligned_matches"][match_count]["target_seq"] = seq[start: end]
					sequences_matches_forward[record.id]["aligned_matches"][match_count]["coord_target"] = [start, end]
				
					if oriented == 5:
						sequences_matches_forward[record.id]["aligned_matches"][match_count]["coord_match"] = [start - len(pam_seq), end]
					elif oriented == 3:
						sequences_matches_forward[record.id]["aligned_matches"][match_count]["coord_match"] = [start, end + len(pam_seq)]	
					
				# increment the count of matches
				match_count += 1

	#######################################
	# Finding targets in the reverse strand
	#######################################

	match_count = 1 			# a variable to give numeric identifiers
	sequences_matches_reverse = {} # a dictionary with the following structure:

	# At the top, there will be sequence identifiers: sequence IDs and "consensus", which will have several pieces 
	# of information in a data structure as their values.
	# More specifically, "align_seq" = sequence of from the alignment object
	# 					 "original_seq" = the sequence from the original data structure containing sequences
	#					 "aligned_matches" = matches of the sgRNA target regex extracted directly from the aligned sequence
	#					 "orig_matches" = matches of the sgRNA target regex extracted indirectly from the original sequence.
	#					 This last part of this data structure is required to be able to identify unique sequences by removing
	#					 them from the list of all possible matches and then checking all the other matches using some alignment
	#					 or fuzznuc-based approaches

	# initialize a dictionary for the consensus sequence 
	sequences_matches_reverse["consensus"] = {}
	sequences_matches_reverse["consensus"]["align_seq"] = str(rev_align_cons["consensus"])
	sequences_matches_reverse["consensus"]["aligned_matches"] = {}

	# initialize a dictionary for all individual aligned sequence
	for id in rev_align_cons.keys():
		if id != "consensus":
			sequences_matches_reverse[id] = {}
			sequences_matches_reverse[id]["align_seq"] = rev_align_cons[id]
		
			sequences_matches_reverse[id]["original_seq"] = rc(sequences[id])
			sequences_matches_reverse[id]["aligned_matches"] = {}
			sequences_matches_reverse[id]["orig_matches"] = {}

	# run regular expression matching of the consensus sequence of the alignment
	# The aim of the following code is to locate sgRNA targets common among all of the aligned sequences

	cons = str(rev_consensus)

	for match in regex.finditer(rev_align_cons["consensus"]):

		# The code to do a more detailed processing of the matches
		if match.group(1).count('X') <= 1:	# to exclude the sequences which contain more than 1 non-consensus site
						
			# populate the nested dictionaries with the data
			start = match.start(1)
			end = match.end(1)
			
			# define a boolean variable, which is True if the target does not span two exons in any of the sequences
			cons_inside = True
			
			# Next, iterate over all of the aligned sequences
			for record in align:
			
				# The coordinates for matches in the original sequences are only different by the number of gaps in a part of that sequence
				# just before the match. Therefore, it should be possible to convert the data between the two parts of the dictionary by a simple 
				# subtraction. Since this is the reverse strand, one needs to subtract the resulting coordinates from the length of the sequence
				# to get the target coordinate relative to the forward strand

				# collect the essential information for the current match in the current sequence
				seq = rc(str(record.seq))
				gaps = int(seq[:match.start()].count('-'))			
				id = record.id
				target_start = len(seq) - (end - gaps)
				target_end = len(seq) - (start - gaps)       # PAM sequence has to be inside an exon as well

				# define a boolean variable to check if the target is inside a particular exon
				# it is defined as True initially because it is a priori more likely not to overlap two exons
				# once one of the conditions is satisfied, it becomes False and subsequent code execution 
				# occurs differently
				inside = False
				
				# obtain the sequence of the target and check that it is inside a particular exon
				target_seq = rc(seq[start: end])
				
				
				# check for the molecule type selected by the user
				# the choice of the "refseq" option would prevent execution of the following code
				# and the cons_inside variable will remain True
				
				if mol_type != "refseq":
					# check if this target does not overlap two exons			
					for exonID in seqs_exons[id].keys():
						
						if target_seq in seqs_exons[id][exonID]:
							inside = True
					
					# Update the variable for the target sequence to be completely inside a particular exon to False if there
					# is a sequence which does not satisfy this
					if inside: 
						pass
					else:
						cons_inside = False
					
			# Test if there were any sequences, where the target overlapped two exons
			if cons_inside:
				# Consensus sequence initialization and storing the data
				# initialize the dictionary for the consensus and each of the matches
				sequences_matches_reverse["consensus"]["aligned_matches"][match_count] = {}

				sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["target_seq"] = match.group(1)
				sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["coord_target"] = [start, end]
			
				if oriented == 5:
					sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["coord_match"] = [start-len(pam_seq), end]
				elif oriented == 3:
					sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["coord_match"] = [start, end + len(pam_seq)]

				# Next, iterate over all of the aligned sequences			
				for record in align:
					#initialization for all of the aligned sequences
					sequences_matches_reverse[record.id]["aligned_matches"][match_count] = {}
					sequences_matches_reverse[record.id]["orig_matches"][match_count] = {}
										
					# matches in the original sequences:
					orig_seq = rc(str(sequences[record.id]))
					seq = rc(str(record.seq))

					gaps = int(seq[:match.start()].count('-'))

					sequences_matches_reverse[record.id]["orig_matches"][match_count]["target_seq"] = orig_seq[start - gaps: end - gaps]
					sequences_matches_reverse[record.id]["orig_matches"][match_count]["coord_target"] = [start - gaps, end - gaps]

					if oriented == 5:
						sequences_matches_reverse[record.id]["orig_matches"][match_count]["coord_match"] = [start - gaps - len(pam_seq), end - gaps]
					elif oriented == 3:
						sequences_matches_reverse[record.id]["orig_matches"][match_count]["coord_match"] = [start - gaps, end - gaps + len(pam_seq)]			
						
					# matches in the aligned sequences:

					sequences_matches_reverse[record.id]["aligned_matches"][match_count]["target_seq"] = seq[start: end]
					sequences_matches_reverse[record.id]["aligned_matches"][match_count]["coord_target"] = [start, end]
				
					if oriented == 5:
						sequences_matches_reverse[record.id]["aligned_matches"][match_count]["coord_match"] = [start - len(pam_seq), end]
					elif oriented == 3:
						sequences_matches_reverse[record.id]["aligned_matches"][match_count]["coord_match"] = [start, end + len(pam_seq)]
				
				
				# increment the count of matches
				match_count += 1			
						
	#############################################################
	#  CLEAN THE FILESYSTEM
	#############################################################

	os.remove(input_file)
	os.remove(output_file)
	os.remove(dnd_file)
	os.rmdir(new_folder)
	
	#############################################################
	# RESULTS
	#############################################################
	
	# store the results of each run
	# 1) number of target sites in the reverse and forward orientations
	# 2) update a count variable for the number of pairs with common target sites
	
	# record the data on each site, categorize them into normal and mismatch-carrying sites
	# think how to do this in the most appropriate way for the subsequent experimental analysis
	
	
	
	# obtain the counts of common target sites
	forw_count = len(sequences_matches_forward["consensus"]["aligned_matches"])
	rev_count = len(sequences_matches_reverse["consensus"]["aligned_matches"])
		
	if (forw_count > 0) or (rev_count > 0):
		results[pair] = {}
		results[pair]["forw_count"] = str(forw_count)
		results[pair]["rev_count"] = str(rev_count)
		results[pair]["gene1"] = gene1
		results[pair]["gene2"] = gene2
		results[pair]["direct_targets"] = {}
		results[pair]["reverse_targets"] = {}
		results[pair]["consensus_direct"] = {}
		results[pair]["consensus_reverse"] = {}
		
		print pair
		
		if forw_count > 0:
			
			# consensus sequence of common target sites in the forward orientation in each pair and the corresponding computed best sequence
			# initialize
			results[pair]["consensus_direct"]["consensus_seq"] = {}
			results[pair]["consensus_direct"]["best_seq"] = {}
			
			print "Consensus - forward"
			for targ_id in sequences_matches_forward["consensus"]["aligned_matches"]:
				# record consensus sequence for each target site
				targ_seq = sequences_matches_forward["consensus"]["aligned_matches"][targ_id]["target_seq"]
				results[pair]["consensus_direct"]["consensus_seq"][targ_id] = targ_seq
				
				# record best computed sequence for each target site
				best_target_seq = compute_best_sgRNA_target(targ_seq, targ_id, sequences_matches_forward)		
				results[pair]["consensus_direct"]["best_seq"][targ_id] = best_target_seq

				print targ_id
				
			
			for gene_id in geneids:
				results[pair]["direct_targets"][gene_id]={}
				
				if gene_id in sequences_matches_forward:
					ident = gene_id + "-forward"
					print ident
					for matchID in sequences_matches_forward[gene_id]["aligned_matches"]:
						if sequences_matches_forward[gene_id]["aligned_matches"][matchID]:
							results[pair]["direct_targets"][gene_id][matchID] =  sequences_matches_forward[gene_id]["aligned_matches"][matchID]["target_seq"]	
							print matchID
							
		if rev_count > 0:		
	
			# consensus sequence of common target sites in the reverse orientation in each pair and the corresponding computed best sequence
			# initialize
			results[pair]["consensus_reverse"]["consensus_seq"] = {}
			results[pair]["consensus_reverse"]["best_seq"] = {}
			
			print "Consensus - reverse"
			for targ_id in sequences_matches_reverse["consensus"]["aligned_matches"]:
				# record consensus sequence for each target site
				targ_seq = sequences_matches_reverse["consensus"]["aligned_matches"][targ_id]["target_seq"]
				results[pair]["consensus_reverse"]["consensus_seq"][targ_id] = targ_seq
				
				# record best computed sequence for each target site
				best_target_seq = compute_best_sgRNA_target(targ_seq, targ_id, sequences_matches_reverse)		
				results[pair]["consensus_reverse"]["best_seq"][targ_id] = best_target_seq
				print targ_id
				
			for gene_id in geneids:
				results[pair]["reverse_targets"][gene_id]={}
				
				if gene_id in sequences_matches_reverse:
					ident = gene_id + "-reverse"
					print ident
					
					for matchID in sequences_matches_reverse[gene_id]["aligned_matches"]:
						if sequences_matches_reverse[gene_id]["aligned_matches"][matchID]:
							results[pair]["reverse_targets"][gene_id][matchID] =  sequences_matches_reverse[gene_id]["aligned_matches"][matchID]["target_seq"]	
							print matchID	
			
#############################################################################################
# Output Phase: write key pieces of information to text files for visualization and reporting
#############################################################################################

results_writeFH = open("ohnolog_target_sites_counts.txt", 'w')
details_writeFH = open("ohnolog_target_sites_details.txt", 'w')

###############################
# Output target site count data 
###############################

header = "Gene pair" + "\t" + "Forward target site count"  + "\t" + "Reverse target site count"  + "\n"
results_writeFH.write(header)

for pair in results:
	results_writeFH.write(pair)
	results_writeFH.write("\t")
	results_writeFH.write(results[pair]["forw_count"])
	results_writeFH.write("\t")
	results_writeFH.write(results[pair]["rev_count"])
	results_writeFH.write("\n")
	
results_writeFH.close()

################################
# Output detailed information 
# on identified target sites
################################
details_header = "Gene1" + "\t" + "Gene2"  + "\t" + "Orientation"  + "\t" + "Consensus_site"  + "\t" + "Best_computed_sequence" + "\t" + "Site in Gene1" + "\t" + "Site in Gene2" + "\n"
details_writeFH.write(details_header)
  
for pair in results:
	# get the genes making up the pair
	(gene1, gene2) = pair.split('\t')
	
	# first go over target sites in the forward orientation
	if results[pair]["forw_count"] > 0:
		# define orientation	
		orientation = "forward"
		
		if results[pair]["consensus_direct"]:
		
			for matchID in results[pair]["consensus_direct"]["consensus_seq"]:
				
				# write each piece of data to the text file
				details_writeFH.write(gene1)
				details_writeFH.write("\t")
				details_writeFH.write(gene2)
				details_writeFH.write("\t")
				details_writeFH.write(orientation)
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["consensus_direct"]["consensus_seq"][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["consensus_direct"]["best_seq"][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["direct_targets"][gene1][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["direct_targets"][gene2][matchID])
				details_writeFH.write("\n")
	
	# next process target sites in the reverse orientation
	if results[pair]["rev_count"] > 0:
		
		# define orientation	
		orientation = "reverse"
		
		if results[pair]["consensus_reverse"]:
			for matchID in results[pair]["consensus_reverse"]["consensus_seq"]:
				# write each piece of data to the text file
				details_writeFH.write(gene1)
				details_writeFH.write("\t")
				details_writeFH.write(gene2)
				details_writeFH.write("\t")
				details_writeFH.write(orientation)
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["consensus_reverse"]["consensus_seq"][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["consensus_reverse"]["best_seq"][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["reverse_targets"][gene1][matchID])
				details_writeFH.write("\t")
				details_writeFH.write(results[pair]["reverse_targets"][gene2][matchID])
				details_writeFH.write("\n")
	