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
import tempfile
import os.path
import random
import string
import sqlite3 as db
import math

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
			regex += "("
			total_length = int(length)
			regex += "[GATC]{%d})" % total_length
	
	if pam_seq == "NGG":
	
		if oriented == 3:
			regex += "[GATCX]GG"
			
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

	
def regex_simple_match(sequences, sgRNA_regex, pam_seq, length, oriented):
	# initialize an identical to the sequences_targets dictionary to store the results
	# of regular expression matches
	sequences_matches = {}
	
	# lookahead functionality to enable overlapping matches
	final_regex = "(?=%s)" % sgRNA_regex
	
	regex = re.compile(final_regex, re.IGNORECASE)
	
	# iterate over all of the sequences in the sequences dictionary
	for id in sequences.keys():
		sequences_matches[id] = {}
		
		# store the sequence information for this particular sequence identifier
		sequences_matches[id]["sequence"] = sequences[id]
		sequences_matches[id]["rc_sequence"] = rc(sequences[id])
		
		# initialize the dictionaries needed to store the CRISPR sgRNA matches
		sequences_matches[id]["direct_targets"] = {}
		sequences_matches[id]["reverse_targets"] = {}
		
		# it is now time to run the input regular expression and store information about
		# all the matches

		# for the forward DNA strand:
		match_count = 0

		for match in regex.finditer(sequences_matches[id]["sequence"]):

			# The code to do a more detailed processing of the matches	
			# create a new matchID 
			match_count += 1
				
			# initialize the nested dictionaries
			sequences_matches[id]["direct_targets"][match_count] = {}
			
			# populate the nested dictionaries with the data
			sequences_matches[id]["direct_targets"][match_count]["target_seq"] = match.group(1)
			start = match.start(1)
			end = match.end(1)

			sequences_matches[id]["direct_targets"][match_count]["coord_target"] = [start, end]

			if oriented == 5:
				sequences_matches[id]["direct_targets"][match_count]["full_match"] = sequences[id][start-len(pam_seq): end]
				sequences_matches[id]["direct_targets"][match_count]["coord_match"] = [start-len(pam_seq), end]
			elif oriented == 3:
				sequences_matches[id]["direct_targets"][match_count]["full_match"] = sequences[id][start: end + len(pam_seq)]
				sequences_matches[id]["direct_targets"][match_count]["coord_match"] = [start, end + len(pam_seq)]

		# for the reverse DNA strand:
		match_count = 0
		
		for match in regex.finditer(sequences_matches[id]["rc_sequence"]):

			# The code to do a more detailed processing of the matches	
			# create a new matchID 
			match_count += 1
			
			# initialize the nested dictionaries
			sequences_matches[id]["reverse_targets"][match_count] = {}

			# populate the nested dictionaries with the data
			sequences_matches[id]["reverse_targets"][match_count]["target_seq"] = match.group(1)
			start = match.start(1)
			end = match.end(1)

			sequences_matches[id]["reverse_targets"][match_count]["coord_target"] = [start, end]
			
			# populate the nested dictionaries with the data
			if oriented == 5:
				sequences_matches[id]["reverse_targets"][match_count]["full_match"] = sequences_matches[id]["rc_sequence"][start-len(pam_seq): end]
				sequences_matches[id]["reverse_targets"][match_count]["coord_match"] = [start-len(pam_seq), end]
			elif oriented == 3:
				sequences_matches[id]["reverse_targets"][match_count]["full_match"] = sequences_matches[id]["rc_sequence"][start: end + len(pam_seq)]
				sequences_matches[id]["reverse_targets"][match_count]["coord_match"] = [start, end + len(pam_seq)]

	return(sequences_matches)
	
def match_test(seq1, seq2):
# returns TRUE if two sequences have <= 2 mismatches
	
	# count mismatches
	mismatch = 0
	
	# iterate over both sequences
	for i in range(len(seq1)):
		if mismatch > 2:
			return(False)
			break
		
		if seq1[i] == seq2[i]:
			pass
		else:
			mismatch += 1
	
	# after iterating over both sequences simultaneously
	# check if there are not more than 2 mismatches
	# and return True to mean that these sequences match well
	if mismatch <= 2:
		return(True)
	

	
def compute_output_unique_targets(sequences, sgRNA_regex, gene_id):
	


	# get all of the matches for all of the sequences
	# the structure of this dictionary can be found in the function regex_simple_match
	
	seq_matches = regex_simple_match(sequences, sgRNA_regex, pam_seq, length, oriented)
	
	# define which gene we are working with
	gene_id = gene_id
	
	# then initialize the dictionary to store unique targets for each sequence
	# It will contain unique targets for each sequence with a label whether the target is in the forward or reverse orientation
	#
	# seq_ID => {targ_id => {sequence => "seq"; orientation => "forward"; coordinates => [x,y]}}
	unique_targets = {}
	
	# match each of the targets against all possible targets in both orientation
	# A target is considered unique if there are 2 or more mismatches between itself and ALL of
	# the other possible targets in other sequences in the alignment
	#
	# To speed up the process of matching, the program will start from 5'-end of the target sequence,
	# iterate over this sequence and count the number of mismatches between the two sequence. As soon as the number of
	# mismatches reaches 2, there will be a break command in the loop which will abort further comparison 
	
	for seq_id in seq_matches.keys():
		
		# unique target ID variable
		uniqID = 0
		
		# initialize a dictionary for a particular sequence identifier
		unique_targets[seq_id] = {}
		
		# if no matching targets are found, add the current target to the list of unique targets
		# use the above function match_test to check the quality of the match
			
		# loop over all of the guide RNA targets in the forward strand
		for targ_id in seq_matches[seq_id]["direct_targets"].keys():
			# initialize a boolean variable which will be used to decide if a target is unique or not
			unique = True
			# get the sequence of the target, which is currently being tested
			curr_seq = seq_matches[seq_id]["direct_targets"][targ_id]["target_seq"]
			
			# loop inside over all of the targets in other sequences in both orientations
			for id in seq_matches.keys():
				if id != seq_id:
					# targets in the forward strand of the sequence
					for targ in seq_matches[id]["direct_targets"].keys():
						# get the sequence against which to match the current sequence
						to_match_seq = seq_matches[id]["direct_targets"][targ]["target_seq"]
						
						# do a test for the current pair of sequences.
						# If the test function output is True, then the target is not unique
						if match_test(curr_seq, to_match_seq):
							unique = False
										
					# targets in the forward strand of the sequence
					for targ in seq_matches[id]["reverse_targets"].keys():
						# get the sequence against which to match the current sequence
						to_match_seq = seq_matches[id]["reverse_targets"][targ]["target_seq"]
						
						# do a test for the current pair of sequences.
						# If the test function output is True, then the target is not unique
						if(match_test(curr_seq, to_match_seq)):
							unique = False					
			

			# test this current target for uniqueness and if it is unique, insert it into the 
			if unique:
			
				inside = False
				exon_N = ""

				target_seq = curr_seq
				
				if seq_id in transcr_exons:
					# check if this target does not overlap two exons			
					for exonID in sorted(gene_exons, key = lambda key: gene_exons[key]["start"]):
						if target_seq in gene_exons[exonID]["exon_seq"]:
							inside = True
							
							# determine exon number where this match occurred and exit the loop
							exon_N = exonID
							break
			
				if inside:
					uniqID += 1
					# initialize a dictionary for a specific target
					unique_targets[seq_id][uniqID] = {}
					# fill this dictionary with the required information
					unique_targets[seq_id][uniqID]["sequence"] = curr_seq
					unique_targets[seq_id][uniqID]["orientation"] = "forward"
					unique_targets[seq_id][uniqID]["coordinates"] = seq_matches[seq_id]["direct_targets"][targ_id]["coord_target"]

					if exon_N:
						unique_targets[seq_id][uniqID]["exon_number"] = exon_N
					exon_N = ""
				

		# loop over all of the guide RNA targets in the reverse strand
		for targ_id in seq_matches[seq_id]["reverse_targets"].keys():
			# initialize a boolean variable which will be used to decide if a target is unique or not
			unique = True
			# get the sequence of the target, which is currently being tested
			curr_seq = seq_matches[seq_id]["reverse_targets"][targ_id]["target_seq"]
			
			# loop inside over all of the targets in other sequences in both orientations
			for id in seq_matches.keys():
				if id != seq_id:
					# targets in the forward strand of the sequence
					for targ in seq_matches[id]["direct_targets"].keys():
						# get the sequence against which to match the current sequence
						to_match_seq = seq_matches[id]["direct_targets"][targ]["target_seq"]
						
						# do a test for the current pair of sequences.
						# If the test function output is True, then the target is not unique
						if match_test(curr_seq, to_match_seq):
							unique = False
										
					# targets in the reverse strand of the sequence
					for targ in seq_matches[id]["reverse_targets"].keys():
						# get the sequence against which to match the current sequence
						to_match_seq = seq_matches[id]["reverse_targets"][targ]["target_seq"]
						
						# do a test for the current pair of sequences.
						# If the test function output is True, then the target is not unique
						if(match_test(curr_seq, to_match_seq)):
							unique = False					
			

			# test this current target for uniqueness and if it is unique, insert it into the dictionary
			if unique:
			
				inside = False
				exon_N = ""

				target_seq = rc(curr_seq)
				
				# check if this target does not overlap two exons	
				if seq_id in transcr_exons:				
					for exonID in sorted(gene_exons, key = lambda key: gene_exons[key]["start"]):
						if target_seq in gene_exons[exonID]["exon_seq"]:
							inside = True

							# determine exon number where this match occurred and exit the loop
							exon_N = exonID
							break
			
				if inside:
					uniqID += 1
					# initialize a dictionary for a specific target
					unique_targets[seq_id][uniqID] = {}
					# fill this dictionary with the required information
					unique_targets[seq_id][uniqID]["sequence"] = curr_seq
					unique_targets[seq_id][uniqID]["orientation"] = "reverse"
					unique_targets[seq_id][uniqID]["coordinates"] = seq_matches[seq_id]["reverse_targets"][targ_id]["coord_target"]

					if exon_N:
						unique_targets[seq_id][uniqID]["exon_number"] = exon_N
					exon_N = ""
					
	####################################
	# Process and output current targets
	####################################
	
	# output a sequence with highlighted targets
	unique_targets_sequences = {}
	
	for seq_id in unique_targets.keys():
		
		if len(unique_targets[seq_id]) == 0:
			pass
		else:
					
			# move over the unique_targets dictionary and other relevant data and 
			# produce the dictionary compatible with the original highlight_targets_output function
			
			unique_targets_sequences[seq_id] = {}
			unique_targets_sequences[seq_id]["sequence"] = seq_matches[seq_id]["sequence"]
			unique_targets_sequences[seq_id]["rc_sequence"] = seq_matches[seq_id]["rc_sequence"]
			unique_targets_sequences[seq_id]["direct_targets"] = {}
			unique_targets_sequences[seq_id]["reverse_targets"] = {}
			
			
			for targ_id in unique_targets[seq_id].keys():
				# check the orientation of each target and populate the corresponding dictionary of unique_targets_sequences
				
				if unique_targets[seq_id][targ_id]["orientation"] == 'forward':
					unique_targets_sequences[seq_id]["direct_targets"][targ_id] = {}
					unique_targets_sequences[seq_id]["direct_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					unique_targets_sequences[seq_id]["direct_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]
						
				elif unique_targets[seq_id][targ_id]["orientation"] == "reverse":
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id] = {}
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]

	# open a file to append the overall data for that particular gene
	overallData_fh = open('total_unique_sites.txt', 'a')

	# detailed data for target sites counts
	detailData_fh = open('detailed_unique_sites.txt', 'a')
					
	
		
	########################
	# OVERALL RESULTS
	########################
	# Part 1 of output: write the GENE, TOTAL # of transcripts and # of transcripts with unique sites
	
	# gene identifier
	overallData_fh.write(gene_id)
	overallData_fh.write('\t')
	
	# total number of transcripts
	overallData_fh.write(str(len(sequences)))
	overallData_fh.write('\t')	
	
	# number of transcripts with unique sites
	transcripts_w_sites = 0
	
	for seq_id in unique_targets:
		if unique_targets[seq_id]:
			transcripts_w_sites = transcripts_w_sites + 1

	overallData_fh.write(str(transcripts_w_sites))
	overallData_fh.write('\t')
	overallData_fh.write('\n')
	
	#########################
	# DETAILED SITES COUNTS
	#########################
	# Part 2: write geneID, transcriptID, count of sense unique target sites and count of anti-sense unique target sites
	# for those transcripts where unique target sites were identified
	for seq_id in unique_targets:
		if unique_targets[seq_id]:
		
			# gene identifier
			detailData_fh.write(gene_id)
			detailData_fh.write('\t')
			
			# transcript identifier
			detailData_fh.write(seq_id)
			detailData_fh.write('\t')
			
			# number of unique target sites in the sense orientation
			if (seq_id in unique_targets_sequences) and ("direct_targets" in unique_targets_sequences[seq_id]):
				detailData_fh.write(str(len(unique_targets_sequences[seq_id]["direct_targets"])))

			detailData_fh.write('\t')			

			# number of unique target sites in the antisense orientation
			if (seq_id in unique_targets_sequences) and ("reverse_targets" in unique_targets_sequences[seq_id]):
				detailData_fh.write(str(len(unique_targets_sequences[seq_id]["reverse_targets"])))

				
			detailData_fh.write('\t')
			detailData_fh.write('\n')
				
################################
# form data variable assignment 
################################

# form fields related to the properties of the CRISPR guide RNA

# first two nucleotides - of some importance in current sgRNA molecules
dinuc = "NN"

# overall length
length = 20

# orientation of the PAM sequence
oriented = 3
oriented = int(oriented)

# There is no point in allowing mismatches in this script, but to keep functions the same
# a mismatch variable has been created and assigned "No" value
mismatch = "No"

# an option where an user can choose whether to use the fixed most common PAM sequence
pam_seq = "NGG"

# obtain the species value from the form
species = "Danio rerio"

# retrieve the molecule type
mol_type = "ensgenes"

########################################
# READ ALL IDENTIFIERS FOR THIS ANALYSIS
########################################
genes_all = []

# read all of the genes from the file
read_fh = open('genes_with_alternative_transcripts.txt', 'r')

for ensg in read_fh:
	ensg = ensg.rstrip('\n')
	genes_all.append(ensg)
 
########################################
# END OF READ
########################################

GENE_COUNT = 0

for gene_id in genes_all:

	GENE_COUNT  += 1
	
	print GENE_COUNT

	########################################
	# VALIDATION FOR AN IDENTIFIER	
	########################################
	# - pattern matching of an identifier to find out what kind of identifier was entered
	# - database query validation of an identifier and the species name

	# clean up the gene_id variable by removing all possible spaces from it
	gene_id = ''.join(gene_id.split())

	# make the identifier uppercase, but then use case-insensitive comparison
	gene_id.upper()

	# also, establish a database connection
	# initialize and populate the dictionary mapping species names and separate species databases
	DBs_species = {}
	DBs_species["Homo sapiens"] = 'human.db'
	DBs_species["Mus musculus"] = 'mouse.db'
	DBs_species["Rattus norvegicus"] = 'rat.db'
	DBs_species["Gallus gallus"] = 'chicken.db'
	DBs_species["Xenopus tropicalis"] = 'xenopus.db'
	DBs_species["Danio rerio"] = 'zebrafish.db'
	DBs_species["Oryzias latipes"] = 'medaka.db'
	DBs_species["Drosophila melanogaster"] = 'drosophila.db'
	DBs_species["Caenorhabditis elegans"] = 'Celegans.db'
	DBs_species["Arabidopsis thaliana"] = 'arabidopsis.db'
	DBs_species["Zea mays"] = 'zea_mays.db'
	DBs_species["Oriza sativa japonica"] = 'rice.db'

	# initialize the connection to the relevant database
	database = DBs_species[species]
	conn = db.connect(database)
	cur = conn.cursor()


	# first check if the gene identifier has a format of an Ensembl ID
	if mol_type == "ensgenes":

		# try to select the data from Genes table on the assumption that the gene identifier is of Ensembl ID format
		cur.execute("SELECT * from Genes WHERE geneid = ? AND species = ? COLLATE NOCASE", (gene_id, species))
		row = cur.fetchone()
		
		if row is None:
			pass
		else: # retrieve all of the data from the Genes database for this particular gene
			gene_id = row[0]
			symbol = row[1]
			# gene_sequence = row[2]   -  in this case, the task is to retrieve the sequences of individual transcripts
			# and not the merged gene sequence
			species = row[3]
			
	elif mol_type == "genes":
	# assume that the gene identifier is a gene symbol

		cur.execute("SELECT * from Genes WHERE symbol = ? COLLATE NOCASE AND species = ? COLLATE NOCASE", (gene_id, species))
		row = cur.fetchone()

		if row is None:
			pass
		else: # retrieve all of the data from the Genes database for this particular gene
			gene_id = row[0]
			symbol = row[1]
			# gene_sequence = row[2]   -  in this case, the task is to retrieve the sequences of individual transcripts
			# and not the merged gene sequence
			species = row[3]

	##########################################################
	# RETRIEVAL OF DATA FROM THE DATABASE
	##########################################################

	# 1. Prepare data structures for the subsequent steps in the process, such as a basic sequences dictionary with transcript sequences,
	# but also exon structures of each transcript.

	# another data structure will only contain the sequence information for transcripts
	# under the transcriptID keys
	transcripts = {}


	# 2. Retrieve all the possible data items from the database and verify that this has worked well.
	#    Consider raising some error if the process does not go as planned.
	
	transcr_exons = {}
	
	# retrieve transcript information
	for row in cur.execute("SELECT * FROM Transcripts WHERE geneid = ?", (gene_id,)):
		transcriptID = row[0]
		gene_id = row[1]
		sequence = row[2]
		sequence = sequence.rstrip()
		
		# update the dictionary with transcripts and exons
		transcr_exons[transcriptID] = {}
		transcr_exons[transcriptID]["gene"] = gene_id
		transcr_exons[transcriptID]["sequence"] = str(sequence)
		
		# update the dictionary for the alignment
		transcripts[transcriptID] = str(sequence)
		
	# retrieve exon information
	# since a gene contains multiple exons, to insert their sequences into the transcripts dictionary,
	# one needs to extract all exons for one gene, store exon sequences in a separate dictionary and 
	# then retrieve the sequences from there

	# store exon sequences and coordinates in this dictionary

	# the program retrieves all the exons of the input
	# and stores its sequence and coordinates

	gene_exons = {}
		
	# extract exon sequences from the database and their coordinates in the gene
	for row in cur.execute("SELECT exonID, sequence, genestart, geneend FROM Exons WHERE geneid = ?", (gene_id,)):
		# data from the current row
		exonID = row[0]
		exon_seq = row[1]
		start = row[2]
		end = row[3]
		
		#print exonID
		#print start
		#print end
		
		if end:
			end = end.rstrip()
		

		
		# check for validity of the data
		if (start is not None) and (end is not None):

			# store these data in the dictionary				
			gene_exons[exonID] = {}
			gene_exons[exonID]["exon_seq"] = exon_seq
		
		# store the data pieces in a dictionary
			if start == "0":
				gene_exons[exonID]["start"] = 0
			else:
				gene_exons[exonID]["start"] = int(start)

			gene_exons[exonID]["end"] = int(end)
			
	# close the connection to the database
	conn.close()		
	
	# to keep the subsequent code consistent with other scripts,
	# the "transcripts" dictionary to the "sequences" dictionary
	sequences = transcripts
	
	
	########################################################################################
	# Output Phase: write relevant data items to text files
	########################################################################################
	
	#  Construct a regular expression out of the available input
	#  components and compile it.
	sgRNA_regex = simple_regex_generate(dinuc, length, pam_seq, oriented, mismatch)
	
	# compute unique sgRNA targets in each sequence
	compute_output_unique_targets(sequences, sgRNA_regex, gene_id)

