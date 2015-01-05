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
import sqlite3 as db
import math

# Biopython
import Bio
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align.Generic import Alignment
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO

# to enable exception handling
cgitb.enable()

################
# functions
################

def calcDoenchScore(seq):
	params = [
	# pasted/typed table from PDF and converted to zero-based positions
	(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
	(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
	(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
	(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
	(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
	(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
	(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
	(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
	(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
	(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
	(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
	(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
	(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
	(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
	(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
	(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
	(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
	(26,'GT',0.11787758),(28,'GG',-0.69774)]
	 
	intercept =  0.59763615
	gcHigh    = -0.1665878
	gcLow     = -0.2026259

	score = intercept
 
	guideSeq = seq[4:24]
	gcCount = guideSeq.count("G") + guideSeq.count("C")
	if gcCount <= 10:
		gcWeight = gcLow
	if gcCount > 10:
		gcWeight = gcHigh
	score += abs(10-gcCount)*gcWeight

	for pos, modelSeq, weight in params:
		subSeq = seq[pos:pos+len(modelSeq)]
		if subSeq==modelSeq:
			score += weight
	return 1.0/(1.0+math.exp(-score))


def print_markers(starts, curr_slice):

	# define the lower boundary for starts
	prev_slice = curr_slice - 100
	
	# store start sites in the current segment
	segment = []
	

	# short line
	short_line = False

	# filter the start sites
	for s in starts:

		if s >= prev_slice and s < curr_slice:
			segment.append(s)

	# start printing spaces and markers  if there are any start sites in the current segment
	if segment:
		
		for s in segment:		
			if s == prev_slice:
				short_line = True

		# initiate the line for markers
		if short_line:
			line = "       "
		else:
			line = "        "

		sys.stdout.write("%s" %line)

		# initialize the most recent printed position
		last = prev_slice	

		for site in segment:
			# define a coordinate of site in the segment
			coord = site - last
		
			num_spaces = coord -1
			last = site

			line = " "*num_spaces
			sys.stdout.write("%s" %line)
			sys.stdout.write("<strong><font style='background-color: green; color:white'>&gt;</font></strong>")
		
		

		num_spaces = curr_slice - last
		line = " "*num_spaces
		sys.stdout.write("%s" %line)
		sys.stdout.write("\n")


def print_markers_double(starts_forward, starts_reverse, curr_slice):

	# define the lower boundary for starts
	prev_slice = curr_slice - 100
	
	# store start sites in the current segment for both forward and reverse strand sites
	segment_forward = []
	segment_reverse = []
	
	# filter the start sites in both orientations
	for s in starts_forward:

		if s >= prev_slice and s < curr_slice:
			segment_forward.append(s)

	for s in starts_reverse:

		if s >= prev_slice and s < curr_slice:
			segment_reverse.append(s)


	# start printing spaces and markers  if there are any start sites in the current segment
	if segment_forward or segment_reverse:
		# the main idea of this function is to take the positions where relevant markers
		# need to be inserted, insert them and then convert back to the string
		
		# initialize a string variable of 100 spaces
		line = " "*100

		line_list = list(line)
		
		# iterate over all forward strand sites in this segment
		for loc in segment_forward:

			pos = loc - prev_slice
			line_list[pos] = "<strong><font style='background-color: green; color:white'>&gt;</font></strong>"

		# iterate over all forward strand sites in this segment
		for loc in segment_reverse:

			pos = loc - prev_slice
			line_list[pos] = "<strong><font style='background-color: magenta; color:white'>&lt;</font></strong>"

		line = "".join(line_list)
		line = "       " + line
	
		sys.stdout.write("%s" %line)
		sys.stdout.write("\n")



def GCpercent(seq):
	# make sure that the sequence is in the string format
	seq = str(seq)
	seq = seq.upper()

	# calculate the percentage
	sum = float(seq.count("G") + seq.count("C"))
	size = float(len(seq))
	GCperc = sum*100.00/size

	# output the result
	return GCperc

def DNA_RNA_Tm(s, dnac=50, saltc=50): 
## Returns DNA/RNA tm using nearest neighbor thermodynamics. 
 
# dnac is DNA concentration [nM] 
# saltc is salt concentration [mM]. 

# Copyright 2004-2008 by Sebastian Bassi. 
# All rights reserved. 
# This code is part of the Biopython distribution and governed by its 
# license.  Please see the LICENSE file that should have been included 
# as part of this package. 
# Code derived from Bio.SeqUtils.MeltingTemp

	def overcount(st, p): 
	# """Returns how many p are on st, works even for overlapping""" 
		ocu = 0 
		x = 0 
		while True: 
			try: 
				i = st.index(p, x) 
			except ValueError: 
				break 
			ocu += 1 
			x = i + 1 
		return ocu

		
	dh = 0  # DeltaH. Enthalpy 
	ds = 0  # deltaS Entropy 
				
	def tercorr(stri): 
		deltah = 0
		deltas = 0
		
		# # initiation parameters
		# dh = -1.9
		# ds = 3.9
		
		dhL = dh + 2*1.9
		dsL = ds - 2*3.9
 
		return dsL, dhL 

	# initialize other parameters	
	R = 1.987  # universal gas constant in Cal/degrees C*Mol 
	k = (dnac/4.0)*1e-9 
	
	sup = str(s).upper()  # turn any Seq object into a string (need index method)
	vs, vh = tercorr(sup)
	
	# Nearest-neighbour calculations
	
	# Enthalpy
	vh = vh+(overcount(sup, "AA"))*11.5+(overcount(sup, "TT"))*7.8 + (overcount(sup, "AT"))*8.3 + (overcount(sup, "TA"))*7.8 + (overcount(sup, "CA"))*10.4 + (overcount(sup, "TG"))*9+ (overcount(sup, "GT"))*5.9 + (overcount(sup, "AC"))*7.8

	vh = vh + (overcount(sup, "CT"))*9.1 + (overcount(sup, "AG"))*7  + (overcount(sup, "GA"))*8.6 + (overcount(sup, "TC"))*5.5 
	vh = vh + (overcount(sup, "CG"))*16.3 + (overcount(sup, "GC"))*8 + (overcount(sup, "GG"))*9.3 + (overcount(sup, "CC"))*12.8

	# Entropy
	vs = vs + (overcount(sup, "AA"))*36.4 + (overcount(sup, "TT"))*21.9 +(overcount(sup, "AT"))*23.9+(overcount(sup, "TA"))*23.2
	vs = vs + (overcount(sup, "CA"))*28.4 + (overcount(sup, "TG"))*26.1 + (overcount(sup, "GT"))*12.3 + (overcount(sup, "AC"))*21.6
	vs = vs + (overcount(sup, "CT"))*23.5 + (overcount(sup, "AG"))*19.7 + (overcount(sup, "GA"))*22.9 + (overcount(sup, "TC"))*13.5
	vs = vs + (overcount(sup, "CG"))*47.1 + (overcount(sup, "GC"))*17.1 + (overcount(sup, "GG"))*23.2 + (overcount(sup, "CC"))*31.9

	ds = vs 
	dh = vh 

	ds = ds-0.368*(len(s)-1)*math.log(saltc/1e3) 
	tm = ((1000* (-dh))/(-ds+(R * (math.log(k)))))-273.15 
	# print("ds=%f" % ds) 
	# print("dh=%f" % dh) 
	return tm


def crispr_targeter_header(sequences):
	
	seq_count = 3*len(sequences)  # maximum number of necessary sections which will have to be hidden and expanded
	
	print """
<!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Simple CRISPR guide RNA target search in an input sequence</title>

<link href="http://www.multicrispr.net/style.css" media="screen" rel="Stylesheet" type="text/css">
<script type="text/javascript" src="http://www.multicrispr.net/view.js"></script>
<style type="text/css">
"""

	for i in range(seq_count):
		print "#ToggleTarget%d {display: none;}" %(i+1)

	print """
	</style>

	<script type="text/javascript">
	"""

	i = 0

	for i in range(seq_count):
		print """
	function Toggle%d() {
		var el = document.getElementById("ToggleTarget%d");
		if (el.style.display == "block") {
			el.style.display = "none";
		}
		else {
			el.style.display = "block";
		}
	}
	""" %(i+1, i + 1)



	print """
	</script>
	</head>
	<body id="main_body">

	<div id="everything">
		<div id="masthead">
			<div id="masthead-inner">          
				<div class="title">
					<a href="http://www.multicrispr.net/index.html"> <img src="http://www.multicrispr.net/CRISPR Targeter.jpg" style="width: 120px; height: 70px; vertical-align: middle;" border="0"> CRISPR MultiTargeter</a>
				</div>
			</div>
		</div>
		<p>                       </p>
		<div id="form_container">
	"""

def crispr_targeter_footer():
	print """

<p style="font-size: 10pt"><font style='color: red'>Off-target analysis guide:</font> <br />
1. Copy an output table from the text area into a spreadsheet program.<br />
2. Select a column from a spreadsheet program containing the target sequences. The parameters were described on the input page.<br />
3. Construct your target rule such N20NGG (for type II sgRNAs).<br />
4. Decide on your mismatch parameter (number of mismatches above which off-targets are not considered).<br />
5. Remember your target genome.<br />
6. Go to either <strong><a href="http://gt-scan.braembl.org.au/gt-scan/submit">GT-Scan<a>: </strong> or <strong><a href="http://www.rgenome.net/cas-offinder/">Cas-OFFinder<a>: </strong> and perform your off-target analysis.<br />
</p>
<p></p>

</div>

<div id="detailed_text">
  <div id="acks">
     <p> CRISPR MultiTargeter was developed by Sergey Prykhozhij at the IWK Health Centre and Dalhousie University. The design of the logo and the buttons was by Vinothkumar Rajan. The latest update was on the 5th of January 2015.
 </p></div>
</div>
</div>
</body>
</html>
"""

def error_header():

	print "Content-type: text/html\n\n"

	print """
<!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Error Page!</title>


<link href="http://www.multicrispr.net/style.css" media="screen" rel="Stylesheet" type="text/css">
<script type="text/javascript" src="http://www.multicrispr.net/view.js"></script>

</head>
<body id="main_body">

<div id="everything">
    <div id="masthead">
        <div id="masthead-inner">          
			<div class="title">
                <a href="http://www.multicrispr.net/index.html"> <img src="http://www.multicrispr.net/CRISPR Targeter.jpg" style="width: 120px; height: 70px; vertical-align: middle;" border="0"> CRISPR MultiTargeter</a>
            </div>
        </div>
    </div>
	<div id="form_container">
"""

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

def badrequest_seq(bad):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The sequence <span style='color: #ff0000;'>'%s'</span> you entered is not compatible with this program. Please try again later.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()

def badrequest_ids(bad):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>It is not possible to process sequence identifiers <span style='color: #ff0000;'>'%s'</span>  when 'Sequence Input' option is chosen. Please use the correct molecule type for your identifiers.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()	

def badrequest_species(bad):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>It is not possible to process sequence identifiers <span style='color: #ff0000;'>'%s'</span>  when 'Species' option is not chosen. Please choose the correct species.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()
	
def badrequest_fasta():	
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>This sequence is not a correct multi-FASTA sequence. Please try again later.</p>")
	print("</body></html>")
	# Get out of here:
	return sys.exit()
	
def badrequest_input(bad):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>Your input <span style='color: #ff0000;'>'%s'</span> has to contain at least 2 sequences. Please try again later.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()
	
	
def badrequest_pam(pam_seq):
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The PAM sequence <span style='color: #ff0000;'>'%s'</span> you entered is not correct. Please use the NGG option or enter your own PAM sequence.</p>" %pam_seq)
	print("</body></html>")
	return sys.exit()
	
def badrequest_oriented(pam_seq):
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The PAM sequence <span style='color: #ff0000;'>'%s'</span> you entered does not allow for this orientation. Please use the NGG option and 3' orientation.</p>" %pam_seq)
	print("</body></html>")
	return sys.exit()	
	
def badrequest_length(bad):
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The Target Length <span style='color: #ff0000;'>'%s'</span> you entered is not correct. Please enter an integer number.</p>" %bad)
	print("</body></html>")
	return sys.exit()

def badrequest_spacer(bad):
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The Spacer Length for nickase targeting <span style='color: #ff0000;'>'%s'</span> you entered is not correct. Please enter an integer number.</p>" %bad)
	print("</body></html>")
	return sys.exit()

def simple_regex_generate(dinuc, length, pam_seq, oriented):
	# initialize the regular expression string 
	regex = ""
	

	# regex string creation based on the regular expression language
	# add the 5'-terminal dinucleotide and the rest of the target
	# according to the target length

	if dinuc == "GG":
		regex += "("
		regex += dinuc
		length_right = int(length) -2
		regex += "[GATC]{%d})" % length_right
	elif dinuc == "GN":
		regex += "("
		regex += "G"
		regex += "[GATC]"
		length_right = int(length) -2
		regex += "[GATC]{%d})" % length_right
	else:
		total_length = int(length)
		regex += "([GATC]{%d})" % total_length
	
	# add the PAM sequence
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

def double_regex_generate(dinuc, length, pam_seq, min_spacer, max_spacer):
	# initialize the regular expression string 
	regex = ""
	
	##############################################
	# Add the right half of the regular expression 
	##############################################
	
	# regex string creation based on the regular expression language
	# add the 5'-terminal dinucleotide and the rest of the target
	# according to the target length

	if dinuc == "GG":
		regex += "("
		regex += dinuc
		length_right = int(length) -2
		regex += "[GATC]{%d})" % length_right
	elif dinuc == "GN":
		regex += "("
		regex += "G"
		regex += "[GATC]"
		length_right = int(length) -2
		regex += "[GATC]{%d})" % length_right
	else:
		total_length = int(length)
		regex += "([GATC]{%d})" % total_length
	
	# add the PAM sequence
	if pam_seq == "NGG":
	
		regex += "[GATC]GG"
	
	else:
		# assume that this alternative PAM sequence has been previously validated
		# and includes only valid IUPAC DNA characters
		# iupac_dna = set("ATGCWSMKRYBDHVN")
		
		for c in pam_seq:
			# use the previously-created dictionary to generate a full regular expression
			# for an alternative PAM sequence
			code_letter = c.upper()
			regex += iupac_dict[code_letter]
	
	###################################################
	# Add the spacer sequence between the half-targets
	###################################################	
	spacer = "[GATC]{%d,%d}" %(min_spacer, max_spacer)
	regex = spacer + regex
	
	##############################################
	# Add the left half of the regular expression 
	##############################################
	
	left_regex = ""
	
	# add the PAM sequence
	if pam_seq == "NGG":
	
		left_regex += "CC[GATC]"
	
	else:
		# assume that this alternative PAM sequence has been previously validated
		# and includes only valid IUPAC DNA characters
		# iupac_dna = set("ATGCWSMKRYBDHVN")
		
		for c in rc(pam_seq):
			# use the previously-created dictionary to generate a full regular expression
			# for an alternative PAM sequence
			code_letter = c.upper()
			left_regex += iupac_dict[code_letter]
	
	# add the left target
	if dinuc == "GG":
		left_regex += "("
		length_right = int(length) -2
		left_regex += "[GATC]{%d}" % length_right
		left_regex += rc(dinuc)
		left_regex += ")"
		
	elif dinuc == "GN":
		left_regex += "("
		length_right = int(length) -2
		left_regex += "[GATC]{%d}" % length_right
		left_regex += "[GATC]"
		left_regex += "C"
		left_regex += ")"		
		
	else:
		total_length = int(length)
		left_regex += "([GATC]{%d})" % total_length
	
	###############################################
	# Combine all the parts and return
	###############################################
	
	regex = left_regex + regex
	
	return(regex)	

	
def rc(dna):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = dna.translate(complements)[::-1]
	return rcseq

	
def regex_simple_match(sequences, sgRNA_regex, pam_seq, length):
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

			# the idea here is to record the coordinates of the target site only
                     # other coordinates will be derived
			start = match.start(1)
			end = match.end(1)

			################################################################
			# Block of code to check if the current match overlaps two exons
			################################################################
			
			inside = False
			exon_N = ""
			
			# check if the current target overlaps two exons
			# only if the molecule type is not "refseq"
			if (mol_type != "refseq") and (text_input == ""): 
				target_seq = sequences[id][start: end]
				
				if id in seqs_exons:
					# check if this target does not overlap two exons
					
					# sort the exons based on their starting coordinate			
					for exonID in sorted(coords_exons[id], key = lambda key: coords_exons[id][key]["start"]):
						if target_seq in seqs_exons[id][exonID]:
							inside = True
							
							# determine exon number where this match occurred
							exon_N = exonID
							break
			else:
			# if the molecule type is "refseq", the program gives the inside variable a True value
				inside = True
			
			###################################
			# End of Block
			###################################
			
			
			if inside:
			
				# initialize the nested dictionaries
				sequences_matches[id]["direct_targets"][match_count] = {}
				
				# populate the nested dictionaries with the data
				sequences_matches[id]["direct_targets"][match_count]["target_seq"] = match.group(1)
				
				if oriented == 5:
					sequences_matches[id]["direct_targets"][match_count]["full_match"] = sequences_matches[id]["sequence"][start-len(pam_seq): end]
					sequences_matches[id]["direct_targets"][match_count]["coord_match"] = [start-len(pam_seq), end]
				elif oriented == 3:
					sequences_matches[id]["direct_targets"][match_count]["full_match"] = sequences_matches[id]["sequence"][start: end + len(pam_seq)]
					sequences_matches[id]["direct_targets"][match_count]["coord_match"] = [start, end + len(pam_seq)]
					
					# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
					# together with its sequence
 
					if (pam_seq == "NGG") and len(match.group(1)) == 20:
						# check that the target site is sufficiently far from the end of the sequence
						if end < len(sequences_matches[id]["sequence"]) - 6:
							sequences_matches[id]["direct_targets"][match_count]["seq_forscore"] = sequences_matches[id]["sequence"][start-4: end + 6]


				sequences_matches[id]["direct_targets"][match_count]["coord_target"] = [start, end]

				if exon_N:
					sequences_matches[id]["direct_targets"][match_count]["exon_number"] = exon_N
				exon_N = ""			

		# for the reverse DNA strand:
		match_count = 0
		
		for match in regex.finditer(sequences_matches[id]["rc_sequence"]):

			# The code to do a more detailed processing of the matches	
			# create a new matchID 
			match_count += 1
			
			# the idea here is to record the coordinates of the target site only
                    # other coordinates will be derived
			start = match.start(1)
			end = match.end(1) 
			################################################################
			# Block of code to check if the current match overlaps two exons
			################################################################
			
			inside = False
			exon_N = ""
					
			# check if the current target overlaps two exons
			# only if the molecule type is not "refseq"
			if (mol_type != "refseq") and (text_input == ""):
				target_seq =  rc(sequences_matches[id]["rc_sequence"][start: end])
				
				if id in seqs_exons:
					# check if this target does not overlap two exons			
					for exonID in seqs_exons[id].keys():
						if target_seq in seqs_exons[id][exonID]:
							inside = True

							# determine exon number where this match occurred
							exon_N = exonID
							break


			else:
			# if the molecule type is "refseq", the program gives the inside variable a True value
				inside = True
			
			###################################
			# End of Block
			###################################
			
			if inside:			
				# initialize the nested dictionaries
				sequences_matches[id]["reverse_targets"][match_count] = {}

				# populate the nested dictionaries with the data
				sequences_matches[id]["reverse_targets"][match_count]["target_seq"] = match.group(1)
				
				if oriented == 5:
					sequences_matches[id]["reverse_targets"][match_count]["full_match"] = sequences_matches[id]["rc_sequence"][start-len(pam_seq): end]
					sequences_matches[id]["reverse_targets"][match_count]["coord_match"] = [start-len(pam_seq), end]
				elif oriented == 3:
					sequences_matches[id]["reverse_targets"][match_count]["full_match"] = sequences_matches[id]["rc_sequence"][start: end + len(pam_seq)]
					sequences_matches[id]["reverse_targets"][match_count]["coord_match"] = [start, end + len(pam_seq)]

					# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
					# together with its sequence
 
					if (pam_seq == "NGG") and len(match.group(1)) == 20:
						# check that the target site is sufficiently far from the end of the sequence
						if end < len(sequences_matches[id]["rc_sequence"]) - 6:
							sequences_matches[id]["reverse_targets"][match_count]["seq_forscore"] = sequences_matches[id]["rc_sequence"][start-4: end + 6]

					
				sequences_matches[id]["reverse_targets"][match_count]["coord_target"] = [start, end]

				if exon_N:
					sequences_matches[id]["reverse_targets"][match_count]["exon_number"] = exon_N
				exon_N = ""


	return(sequences_matches)
	
def regex_double_match(sequences, sgRNA_regex, pam_seq, length):
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
		pair_count = 1
		
		
		for match in regex.finditer(sequences_matches[id]["sequence"]):

			# The code to do a more detailed processing of the matches	
			# create a new matchID 
			match_count += 1
			start1 = match.start(1) - 3
			end2 = match.end(2) + 3
			
			
			################################################################
			# Block of code to check if the current match overlaps two exons
			################################################################
			
			inside = False
			exon_N = ""

			# check if the current target overlaps two exons
			# only if the molecule type is not "refseq"
			if (mol_type != "refseq") and (text_input == ""):
				target_seq = sequences_matches[id]["sequence"][start1:end2]
				
				if id in seqs_exons:
					# check if this target does not overlap two exons			
					for exonID in seqs_exons[id].keys():
						if target_seq in seqs_exons[id][exonID]:
							inside = True

							# determine exon number where this match occurred
							exon_N = exonID
							break

			else:
			# if the molecule type is "refseq", the program gives the inside variable a True value
				inside = True
			
			###################################
			# End of Block
			###################################			
			
			if inside:
			
				# initialize the nested dictionaries
				sequences_matches[id]["direct_targets"][match_count] = {}
				
				# populate the nested dictionaries with the data
				
				#############
				# FIRST MATCH
				#############

				# the first match is in the reverse strand of DNA so the target sequence will be from the reverse strand 
				sequences_matches[id]["direct_targets"][match_count]["target_seq"] = rc(match.group(1))
				start = match.start(1)
				end = match.end(1)
				sequences_matches[id]["direct_targets"][match_count]["coord_target"] = [start, end]
				sequences_matches[id]["direct_targets"][match_count]["pair_id"] = pair_count

				if (pam_seq == "NGG") and len(match.group(1)) == 20:
					sequences_matches[id]["direct_targets"][match_count]["seq_forscore"] = rc(sequences_matches[id]["sequence"][start-6: end + 4])

				if exon_N:
					sequences_matches[id]["direct_targets"][match_count]["exon_number"] = exon_N

				##############
				# SECOND MATCH
				##############
				
				# iterate the match identifier and move on to the next half target
				match_count += 1
					
				# initialize the nested dictionaries
				sequences_matches[id]["direct_targets"][match_count] = {}
							
				sequences_matches[id]["direct_targets"][match_count]["target_seq"] = match.group(2)
				start = match.start(2)
				end = match.end(2)
				sequences_matches[id]["direct_targets"][match_count]["coord_target"] = [start, end]
				sequences_matches[id]["direct_targets"][match_count]["pair_id"] = pair_count

				if (pam_seq == "NGG") and len(match.group(2)) == 20:

					# check that the target site is sufficiently far from the end of the sequence
					if end < len(sequences_matches[id]["sequence"]) - 6:
						sequences_matches[id]["direct_targets"][match_count]["seq_forscore"] = sequences_matches[id]["sequence"][start-4: end + 6]

				if exon_N:
					sequences_matches[id]["direct_targets"][match_count]["exon_number"] = exon_N				

				pair_count += 1
				exon_N = ""

		# for the reverse DNA strand:
		match_count = 0
		pair_count = 1
		
		for match in regex.finditer(sequences_matches[id]["rc_sequence"]):

			# The code to do a more detailed processing of the matches	
			# create a new matchID 
			match_count += 1
			
			start1 = match.start(1) - 3
			end2 = match.end(2) + 3
			################################################################
			# Block of code to check if the current match overlaps two exons
			################################################################
			
			inside = False
			exon_N = ""
			
			# check if the current target overlaps two exons
			# only if the molecule type is not "refseq"
			if (mol_type != "refseq") and (text_input == ""):
				target_seq = rc(sequences_matches[id]["rc_sequence"][start1:end2])
				
				if id in seqs_exons:
					# check if this target does not overlap two exons			
					for exonID in seqs_exons[id].keys():
						if target_seq in seqs_exons[id][exonID]:
							inside = True

							# determine exon number where this match occurred
							exon_N = exonID
							break

			else:
			# if the molecule type is "refseq", the program gives the inside variable a True value
				inside = True
			
			###################################
			# End of Block
			###################################		
			
			if inside:
				##############
				# FIRST MATCH
				##############
				
				# initialize the nested dictionaries
				sequences_matches[id]["reverse_targets"][match_count] = {}
				
				# populate the nested dictionaries with the data
				
				# the first match is in the reverse strand of DNA so the target sequence will be from the reverse strand 
				sequences_matches[id]["reverse_targets"][match_count]["target_seq"] = rc(match.group(1))
				start = match.start(1)
				end = match.end(1)
				sequences_matches[id]["reverse_targets"][match_count]["coord_target"] = [start, end]
				sequences_matches[id]["reverse_targets"][match_count]["pair_id"] = pair_count

				if (pam_seq == "NGG") and len(match.group(1)) == 20:
					sequences_matches[id]["reverse_targets"][match_count]["seq_forscore"] = rc(sequences_matches[id]["rc_sequence"][start-6: end + 4])

				if exon_N:
					sequences_matches[id]["reverse_targets"][match_count]["exon_number"] = exon_N

				##############
				# SECOND MATCH
				##############
			
				# iterate the match identifier and move on to the next half target
				match_count += 1
				
				# initialize the nested dictionaries
				sequences_matches[id]["reverse_targets"][match_count] = {}
				
				sequences_matches[id]["reverse_targets"][match_count]["target_seq"] = match.group(2)
				start = match.start(2)
				end = match.end(2)
				sequences_matches[id]["reverse_targets"][match_count]["coord_target"] = [start, end]

				if (pam_seq == "NGG") and len(match.group(2)) == 20:

					# check that the target site is sufficiently far from the end of the sequence
					if end < len(sequences_matches[id]["rc_sequence"]) - 6:
						sequences_matches[id]["reverse_targets"][match_count]["seq_forscore"] = sequences_matches[id]["rc_sequence"][start-4: end + 6]
				
				if exon_N:
					sequences_matches[id]["reverse_targets"][match_count]["exon_number"] = exon_N


				sequences_matches[id]["reverse_targets"][match_count]["pair_id"] = pair_count
				
				pair_count += 1
				exon_N = ""

	return(sequences_matches)

	
def highlight_targets_output(sequences_targets, design_type): 
 	# make a variable for the design type
	type = design_type # single or double
	
	id_count = 1
	
	for seq in sequences_targets.keys():
		# print basic headers 
		
		print """
		<h4><a href="javascript:Toggle%d();">The CRISPR sgRNA targets in %s sequence - expand or hide</a></h4>
		<div id="ToggleTarget%d">
		""" %(id_count, seq[0:40], id_count)
	
		id_count += 1
		
		#####################
		# Ensembl link output
		#####################
		

		if (species != "Oriza sativa japonica") and (species != "different") and (species != "default"):
			(genus, sp) = species.split()


		animals = ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Gallus gallus", "Xenopus tropicalis", "Danio rerio", "Oryzias latipes", "Drosophila melanogaster", "Caenorhabditis elegans"]
		plants = ["Arabidopsis thaliana", "Oriza sativa japonica", "Zea mays"]

		if text_input == "":

			# deal with different species
			if species in animals:
				print("<h4><a style='color: #0000FF;' href='http://useast.ensembl.org/%s_%s/Gene/Summary?db=core;g=%s'>%s Ensembl Link</a></h4>" %(genus, sp, seq, seq))
			elif species in plants:
				print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(seq, seq))
			elif species == "different":
				
				if "AT" in seq:
					print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(seq, seq))
				elif "OS" in seq:
					print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(seq, seq))
				elif "GRMZM" in seq:
					print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(seq, seq))
				else:
					print("<h4><a style='color: #0000FF;' href='http://useast.ensembl.org/Gene/Summary?db=core;g=%s'>%s Ensembl Link</a></h4>" %(seq, seq))
				

		#####################
		# END of OUTPUT
		#####################
		
		print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Visual View of sgRNA Targets</h3></div>")
		print("<p style='font-family: Times New Roman, Times, serif; font-size: 18px;'><strong>Guide RNA targets have been highlighted. Long stretches of highlighted sequence are due to target overlap.<br />Please see tables for individual sequences.</strong></p>")
		################################
		# Forward strand
		################################
		
		print("<div style='background-color:#0000DD; color: white; width: 500px; height: 25px; margin:10px;'><p class='normalp'><strong>Targets in the forward strand:</strong></p></div>")
		
		# print the <pre> element which will contain all of the sequence
		print("<pre style='margin:10px; padding:5px; line-height: 1.6em;'>")
		
		# store the sequence from the input dictionary in a variable
		# to make it easier to work with these data
		sequence = sequences_targets[seq]["sequence"]
		direct_targets = {}
		direct_targets = sequences_targets[seq]["direct_targets"]
	
		# do a test for the possibility that no targets of sgRNA have been identified
		if len(direct_targets) == 0:
			print("This sequence does not contain any sgRNA targets in the sense orientation.")
	
		# iterate over all the targets of CRISPR sgRNAs in the current sequence
		# 
		# Use the following styles for the sgRNA targets and adjacent PAM sequences:
		# <pre>
		# <strong><font style="BACKGROUND-COLOR: #0000CC; color: white">sgRNA target</font></strong>
		# </pre>
		

		######################################################
		#  Start site marker labeling
		######################################################
		
		# Labeling 
		# <strong><font style="background-color: green; color:white">&rArr;</font></strong>

		# generate a list of target start sites
		start_coords = []

		# looping over sgRNA targets:
		for targ_id in sorted(direct_targets, key = lambda key: direct_targets[key]["coord_target"][0]):
			begin = direct_targets[targ_id]["coord_target"][0]
			start_coords.append(begin)

		# prepare the data for the case that we have to label nickase sites
		if type == "double":
			start_forward = []
			start_reverse = []
			
			for targ_id in sorted(direct_targets, key = lambda key: direct_targets[key]["coord_target"][0]):

				if is_odd(targ_id):
					start_reverse.append(direct_targets[targ_id]["coord_target"][1] - 1)	
				else:
					start_forward.append(direct_targets[targ_id]["coord_target"][0])

		######################################################
		# End of start coordinates list of coordinates
		######################################################

		# initialize the positions to use while producing output
		curr_pos = 0
		curr_slice = 100

		# print the first line of markers if they are available
		if type == "single":
			print_markers(start_coords, curr_slice)
		elif type == "double":
			print_markers_double(start_forward, start_reverse, curr_slice)				


		# print the initial coordinate
		sys.stdout.write('{0:5d}'.format(curr_slice - 99))
		sys.stdout.write("  ")
	
		# looping over sgRNA targets:
		for targ_id in sorted(direct_targets, key = lambda key: direct_targets[key]["coord_target"][0]):
		
			# Case 1: the start position of the target is smaller than the current slice
			
			if direct_targets[targ_id]["coord_target"][0] < curr_slice:
				# print the sequence and the target and prepare the ground for the remaining sequence
				
				# first print the sequence leading up to the sgRNA match
				first = direct_targets[targ_id]["coord_target"][0]
				
				# check for an overlap between the previous target and the new one
				if first < curr_pos:
					first = curr_pos
				else:
					sys.stdout.write("%s" %sequence[curr_pos: first]) # highlighting is not necessary
					curr_pos = first
				
				# check where the target ends
				if direct_targets[targ_id]["coord_target"][1] < curr_slice:
					last = direct_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])      # mark with highlighted font
					curr_pos = last
					
				elif direct_targets[targ_id]["coord_target"][1] == curr_slice:
					last = direct_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])	# mark with highlighted font
					sys.stdout.write("\n")
					
					# increment the relevant positions
					curr_pos = curr_slice
					curr_slice += 100
					
					# print the markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)	
					

					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
				elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = direct_targets[targ_id]["coord_target"][1]


					# print the markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)

					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					curr_pos = last
				
			# Case 2: the start position of the target is equal or larger than the current slice  
			if direct_targets[targ_id]["coord_target"][0] >= curr_slice:
			
				first = direct_targets[targ_id]["coord_target"][0]
				sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
				sys.stdout.write("\n")
				
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers
				if type == "single":
					print_markers(start_coords, curr_slice)
				elif type == "double":
					print_markers_double(start_forward, start_reverse, curr_slice)	


				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")				
				
				# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
				while first >= curr_slice:
					
					sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100

					# print the markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)

					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")				

				sys.stdout.write("%s" %sequence[curr_pos: first]) # highlighting is not necessary
				curr_pos = first
				
				# check where the target ends
				if direct_targets[targ_id]["coord_target"][1] < curr_slice:
					last = direct_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])      # mark with highlighted font
					curr_pos = last
					
				elif direct_targets[targ_id]["coord_target"][1] == curr_slice:
					last = direct_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])	# mark with highlighted font
					sys.stdout.write("\n")
					
					# increment the relevant positions
					curr_pos = curr_slice
					curr_slice += 100

					# print the markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)

					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
				elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = direct_targets[targ_id]["coord_target"][1]

					# print the markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)	

					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")	
					
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					curr_pos = last
		
		# finishing the end of the sequence
		 
		if curr_slice >= len(sequence):
			# output the rest of the sequence
			sys.stdout.write("%s" %sequence[curr_pos: ])
			sys.stdout.write("\n")
		else:
			# finish the current line
			sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
			sys.stdout.write("\n")
			# update the variables
			curr_pos = curr_slice
			curr_slice += 100
			
			sys.stdout.write('{0:5d}'.format(curr_slice - 99))
			sys.stdout.write("  ")				
			
			
			# go over any other remaining lines and print them
			while curr_slice < len(sequence):
				# finish the current line
				sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
				sys.stdout.write("\n")
				# update the variables
				curr_pos = curr_slice
				curr_slice += 100
				
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")
				
			if curr_slice >= len(sequence):
				# output the rest of the sequence
				sys.stdout.write("%s" %sequence[curr_pos: ])
				sys.stdout.write("\n")

		print("</pre>")

		##################
		# Reverse strand
		##################

		print("<div style='background-color: #CF5300; color: white; width: 500px; height: 25px; margin:10px;'><p class='normalp'><strong>Targets in the reverse strand:</strong></p></div>")
		
		# print the <pre> element which will contain all of the sequence
		print("<pre style='margin:10px; padding:5px;line-height: 1.6em;'>")
		
		# store the sequence from the input dictionary in a variable
		# to make it easier to work with these data
		sequence = sequences_targets[seq]["rc_sequence"]
		reverse_targets = {}
		reverse_targets = sequences_targets[seq]["reverse_targets"]
	
		# do a test for the possibility that no targets of sgRNA have been identified
		if len(reverse_targets) == 0:
			print("This sequence does not contain any sgRNA targets in the anti-sense orientation.")
	
		# iterate over all the targets of CRISPR sgRNAs in the current sequence
		# 
		# Use the following styles for the sgRNA targets and adjacent PAM sequences:
		# <pre>
		# <strong><font style="BACKGROUND-COLOR: #0000CC; color: white">sgRNA target</font></strong>
		# </pre>

		######################################################
		#  Start site marker labeling
		######################################################
		
		# Labeling 
		# <strong><font style="background-color: green; color:white">&rArr;</font></strong>

		# generate a list of target start sites
		start_coords = []

		# looping over sgRNA targets:
		for targ_id in sorted(reverse_targets, key = lambda key: reverse_targets[key]["coord_target"][0]):
			begin = reverse_targets[targ_id]["coord_target"][0]
			start_coords.append(begin)

		# prepare the data for the case that we have to label nickase sites
		if type == "double":
			start_forward = []
			start_reverse = []
			
			for targ_id in sorted(reverse_targets, key = lambda key: reverse_targets[key]["coord_target"][0]):

				if is_odd(targ_id):
					start_reverse.append(reverse_targets[targ_id]["coord_target"][1] - 1)	
				else:
					start_forward.append(reverse_targets[targ_id]["coord_target"][0])


		######################################################
		# End of start coordinates list of coordinates
		######################################################

		# initialize the positions to use while producing output
		curr_pos = 0
		curr_slice = 100

		# print the first line of markers if they are available
		if type == "single":
			print_markers(start_coords, curr_slice)
		elif type == "double":
			print_markers_double(start_forward, start_reverse, curr_slice)

		# initialize the positions to use while producing output
		curr_pos = 0
		curr_slice = 100
		
		# print the initial coordinate
		sys.stdout.write('{0:5d}'.format(curr_slice - 99))
		sys.stdout.write("  ")		
		
		# looping over sgRNA targets:
		for targ_id in sorted(reverse_targets, key = lambda key: reverse_targets[key]["coord_target"][0]):
			# Case 1: the start position of the target is smaller than the current slice
			
			if reverse_targets[targ_id]["coord_target"][0] < curr_slice:
				
				# first print the sequence leading up to the sgRNA match
				first = reverse_targets[targ_id]["coord_target"][0]
				
				# print the sequence and the target and prepare the ground for the remaining sequence
				
				# check for an overlap between the previous target and the new one
				if first < curr_pos:
					first = curr_pos
				else:
					sys.stdout.write("%s" %sequence[curr_pos: first]) # highlighting is not necessary
					curr_pos = first
				
				# check where the target ends
				if reverse_targets[targ_id]["coord_target"][1] < curr_slice:
					last = reverse_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])      # mark with highlighted font
					curr_pos = last
					
				elif reverse_targets[targ_id]["coord_target"][1] == curr_slice:
					last = reverse_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])	# mark with highlighted font
					sys.stdout.write("\n")
					
					# increment the relevant positions
					curr_pos = curr_slice
					curr_slice += 100
					
					# print markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)
					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")	
					
				elif reverse_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = reverse_targets[targ_id]["coord_target"][1]

					# print markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)

					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					curr_pos = last
				
			# Case 2: the start position of the target is equal or larger than the current slice  
			if reverse_targets[targ_id]["coord_target"][0] >= curr_slice:
			
				first = reverse_targets[targ_id]["coord_target"][0]
				sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
				sys.stdout.write("\n")
				
				curr_pos = curr_slice
				curr_slice += 100

				# print markers
				if type == "single":
					print_markers(start_coords, curr_slice)
				elif type == "double":
					print_markers_double(start_forward, start_reverse, curr_slice)
				
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")
				
				# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
				while first >= curr_slice:
					
					sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100

					# print markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)
					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
				sys.stdout.write("%s" %sequence[curr_pos: first]) # highlighting is not necessary
				curr_pos = first
				
				# check where the target ends
				if reverse_targets[targ_id]["coord_target"][1] < curr_slice:
					last = reverse_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])      # mark with highlighted font
					curr_pos = last
					
				elif reverse_targets[targ_id]["coord_target"][1] == curr_slice:
					last = reverse_targets[targ_id]["coord_target"][1]
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])	# mark with highlighted font
					sys.stdout.write("\n")
					
					# increment the relevant positions
					curr_pos = curr_slice
					curr_slice += 100

					# print markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)
					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
				elif reverse_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = reverse_targets[targ_id]["coord_target"][1]

					# print markers
					if type == "single":
						print_markers(start_coords, curr_slice)
					elif type == "double":
						print_markers_double(start_forward, start_reverse, curr_slice)

					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					curr_pos = last
		
		# finishing the end of the sequence
		 
		if curr_slice >= len(sequence):
			# output the rest of the sequence
			sys.stdout.write("%s" %sequence[curr_pos: ])
			sys.stdout.write("\n")
		else:
			# finish the current line
			sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
			sys.stdout.write("\n")
			# update the variables
			curr_pos = curr_slice
			curr_slice += 100

			# print markers
			if type == "single":
				print_markers(start_coords, curr_slice)
			elif type == "double":
				print_markers_double(start_forward, start_reverse, curr_slice)
			
			sys.stdout.write('{0:5d}'.format(curr_slice - 99))
			sys.stdout.write("  ")			
			
			# go over any other remaining lines and print them
			while curr_slice < len(sequence):
				# finish the current line
				sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
				sys.stdout.write("\n")
				# update the variables
				curr_pos = curr_slice
				curr_slice += 100

				# print markers
				if type == "single":
					print_markers(start_coords, curr_slice)
				elif type == "double":
					print_markers_double(start_forward, start_reverse, curr_slice)
				
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")
				
			if curr_slice >= len(sequence):
				# output the rest of the sequence
				sys.stdout.write("%s" %sequence[curr_pos: ])
				sys.stdout.write("\n")

		print("</pre>")
	
		################################################################
		# Print all of the targets in both orientations of the sequence
		################################################################
		if type == "single":
			if len(direct_targets) > 0 or len(reverse_targets) > 0:
				targets_table_output(direct_targets, reverse_targets, id_count)
				id_count = id_count + 2
		elif type == "double":
			if len(direct_targets) > 0 or len(reverse_targets) > 0:
				double_targets_table_output(direct_targets, reverse_targets, id_count)
				id_count = id_count + 2
		
		
		print("</div>")
		print("<br />")
		
		# the color labeling of sgRNA target sequences to use when the complete function will be written
		# <pre>The due date is <strong><font
        # style="BACKGROUND-COLOR: #0000CC; color: white"> next </font></strong><strong><font
        # style="BACKGROUND-COLOR: #CF5300; color: white">week</strong>.</pre>


def targets_table_output(direct_targets, reverse_targets, id_count):
	id_count = id_count
	
	# first print the colored heading
	print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Table View of sgRNA Targets</h3></div>")
	
	
	if direct_targets: # checks that the direct_targets dictionary not empty
		print """
			<h4><a href="javascript:Toggle%d();">Table of the targets in the forward strand - expand or hide</a></h4>
			<div id="ToggleTarget%d">
			""" %(id_count, id_count)
			
		id_count += 1
			
		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# all of the target data goes in here
		print """
				<table >
					<tr>
						<td>
							Target ID
						</td>
						<td>
							Target sequence
						</td>
		"""

		if (pam_seq == "NGG") and (length == 20):		
			print """
						<td >
							NNNN-Target-PAM-NNN
						</td>
			"""

		print """
						<td>
							Target start
						</td>
						<td>
							Target end
						</td>
						<td>
							GC% of Target site
						</td>
						<td>
							Tm of sgRNA:DNA
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in direct_targets.keys():
			if "exon_number" in direct_targets[targ_id].keys():
				exon_info_available = True

		if exon_info_available:
			print """
						<td>
							Exon number
						</td>

			"""
		
		if (pam_seq == "NGG") and (length == 20):
			print """
						<td>
							Score
						</td>
			"""

		print("</tr>")
							
		for targ_id in sorted(direct_targets, key= lambda key: direct_targets[key]["coord_target"][0]):
			
			# generate a row
			print("<tr>")
			print('<td>%d</td>' %targ_id)
			print('<td>%s</td>' %direct_targets[targ_id]["target_seq"])
			
			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					print('<td>%s</td>' %direct_targets[targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")
			
			print('<td>%d</td>' %(direct_targets[targ_id]["coord_target"][0] + 1))
			print('<td>%d</td>' %(direct_targets[targ_id]["coord_target"][1] + 1))

			print('<td>')
			print "{0:3.2f}".format(GCpercent(direct_targets[targ_id]["target_seq"]))
			print('</td>')

			print('<td>')
			print "{0:3.2f}".format(DNA_RNA_Tm(direct_targets[targ_id]["target_seq"]))
			print('</td>')
			
			if "exon_number" in direct_targets[targ_id].keys():
				print('<td>%s</td>' %direct_targets[targ_id]["exon_number"])
			
			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					print('<td>')
					print "{0:3.2f}".format(calcDoenchScore(direct_targets[targ_id]["seq_forscore"]))
					print('</td>')
				else:
					print('<td>%s</td>' % "")
				
			print("</tr>")

		print("</table>")	
		print("</div>")
		print("<p class ='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>")


		#output the textarea with the same information as above
		
		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=160 class=mono>
		"""

		sys.stdout.write('%3s' % "ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
			sys.stdout.write('\t')

		sys.stdout.write('%5s' % "Start")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "End")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		# initialize the boolean variable exon_info_available
		exon_info_available = False

		# print the header for Exon number 
		for targ_id in direct_targets.keys():
			if "exon_number" in direct_targets[targ_id].keys():
				exon_info_available = True
				
		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('\t')
			sys.stdout.write('%6s' % "Score")
				
		sys.stdout.write('\n')

		for targ_id in sorted(direct_targets, key= lambda key: direct_targets[key]["coord_target"][0]):
			
			# generate a row
			sys.stdout.write('%3d' %targ_id)
			sys.stdout.write('\t')
			sys.stdout.write('%25s' % direct_targets[targ_id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					sys.stdout.write('%35s' % direct_targets[targ_id]["seq_forscore"])
					sys.stdout.write('\t')
				else:
					sys.stdout.write('%35s' % "")
					sys.stdout.write('\t')

			sys.stdout.write('%5d' %(direct_targets[targ_id]["coord_target"][0] + 1))
			sys.stdout.write('\t')
			sys.stdout.write('%4d' %(direct_targets[targ_id]["coord_target"][1] + 1))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(GCpercent(direct_targets[targ_id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(direct_targets[targ_id]["target_seq"])))


			if "exon_number" in direct_targets[targ_id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %direct_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					sys.stdout.write('\t')
					sys.stdout.write("{0:3.2f}".format(calcDoenchScore(direct_targets[targ_id]["seq_forscore"])))
				else:
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "")
			
			sys.stdout.write('\n')

		print("</textarea>")
		print("</div>")
	
	if reverse_targets: # checks that the reverse_targets dictionary not empty
		print """
			<h4><a href="javascript:Toggle%d();">Table of the targets in the reverse strand - expand or hide</a></h4>
			<div id="ToggleTarget%d">
			""" %(id_count, id_count)
			
		id_count += 1	
		
		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# all of the target data goes in here
		print """
				<table >
					<tr>
						<td>
							Target ID
						</td>
						<td>
							Target sequence
						</td>
		"""

		if (pam_seq == "NGG") and (length == 20):		
			print """
						<td >
							NNNN-Target-PAM-NNN
						</td>
			"""

		print """
						<td>
							Target start
						</td>
						<td>
							Target end
						</td>
						<td>
							GC% of Target site
						</td>
						<td>
							Tm of sgRNA:DNA
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in reverse_targets.keys():
			if "exon_number" in reverse_targets[targ_id].keys():
				exon_info_available = True

		if exon_info_available:
			print """
						<td>
							Exon number
						</td>
			"""

		if (pam_seq == "NGG") and (length == 20):
			print """
						<td>
							Score
						</td>
			"""

		print("</tr>")
		
		for targ_id in sorted(reverse_targets, key= lambda key: reverse_targets[key]["coord_target"][0]):
			
			# generate a row
			print("<tr>")
			# generate all of the table cells
			print('<td>%d</td>' %targ_id)
			print('<td>%s</td>' %reverse_targets[targ_id]["target_seq"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					print('<td>%s</td>' %reverse_targets[targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")

			print('<td>%d</td>' %(reverse_targets[targ_id]["coord_target"][0] + 1))
			print('<td>%d</td>' %(reverse_targets[targ_id]["coord_target"][1] + 1))

			print('<td>')
			print "{0:3.2f}".format(GCpercent(reverse_targets[targ_id]["target_seq"]))
			print('</td>')

			print('<td>')
			print "{0:3.2f}".format(DNA_RNA_Tm(reverse_targets[targ_id]["target_seq"]))
			print('</td>')

			if "exon_number" in reverse_targets[targ_id].keys():
				print('<td>%s</td>' %reverse_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					print('<td>')
					print "{0:3.2f}".format(calcDoenchScore(reverse_targets[targ_id]["seq_forscore"]))
					print('</td>')
				else:
					print('<td>%s</td>' % "")
			print("</tr>")
			
		print("</table>")	
		print('</div>')
		print("<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>")

		#output the textarea with the same information as above	
		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=160 class=mono>
		"""
		sys.stdout.write('%3s' % "ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')

		if (pam_seq == "NGG") and (length == 20):
			
			sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
			sys.stdout.write('\t')


		sys.stdout.write('%5s' % "Start")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "End")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		# initialize the boolean variable exon_info_available
		exon_info_available = False

		# print the header for Exon number 
		for targ_id in reverse_targets.keys():
			if "exon_number" in reverse_targets[targ_id].keys():
				exon_info_available = True
				
		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):

			sys.stdout.write('\t')
			sys.stdout.write('%6s' % "Score")
	
		sys.stdout.write('\n')

		for targ_id in sorted(reverse_targets, key= lambda key: reverse_targets[key]["coord_target"][0]):
			
			# generate a row
			sys.stdout.write('%3d' %targ_id)
			sys.stdout.write('\t')
			sys.stdout.write('%25s' % reverse_targets[targ_id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					sys.stdout.write('%35s' % reverse_targets[targ_id]["seq_forscore"])
					sys.stdout.write('\t')
				else:
					sys.stdout.write('%35s' % "")
					sys.stdout.write('\t')

			sys.stdout.write('%5d' %(reverse_targets[targ_id]["coord_target"][0] + 1))
			sys.stdout.write('\t')
			sys.stdout.write('%4d' %(reverse_targets[targ_id]["coord_target"][1] + 1))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(GCpercent(reverse_targets[targ_id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(reverse_targets[targ_id]["target_seq"])))


			if "exon_number" in reverse_targets[targ_id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %reverse_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					sys.stdout.write('\t')
					sys.stdout.write("{0:3.2f}".format(calcDoenchScore(reverse_targets[targ_id]["seq_forscore"])))
				else:	
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "")

			sys.stdout.write('\n')

		print("</textarea>")
		print('</div>')

def is_odd(num):
   return num % 2 != 0	
	
	
def double_targets_table_output(direct_targets, reverse_targets, id_count):
	id_count = id_count
	
	# first print the colored heading
	print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Table View of sgRNA Targets</h3></div>")

	if direct_targets:	
		print """
		<h4><a href="javascript:Toggle%d();">Table of sgRNA pairs for nickase targeting in the forward strand - expand or hide</a></h4>
		<div id="ToggleTarget%d">
		""" %(id_count, id_count)
		
		id_count += 1	

		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# all of the target data goes in here
		print """
				<table >
					<tr>
						<td>
							Pair ID
						</td>					

						<td>
							Target sequence
						</td>
		"""

		if (pam_seq == "NGG") and (length == 20):		
			print """
						<td >
							NNNN-Target-PAM-NNN
						</td>
			"""

		print """
						<td>
							Target start
						</td>
						<td>
							Target end
						</td>
						<td>
							Strand
						</td>
						<td>
							GC% of Target site
						</td>
						<td>
							Tm of sgRNA:DNA
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in direct_targets.keys():
			if "exon_number" in direct_targets[targ_id].keys():
				exon_info_available = True

		if exon_info_available:
			print """
						<td>
							Exon number
						</td>

			"""

		if (pam_seq == "NGG") and (length == 20):
			print """
						<td>
							Score
						</td>
			"""


		print("</tr>")

		for targ_id in sorted(direct_targets, key= lambda key: direct_targets[key]["pair_id"]):
		
			# generate a row
			print("<tr>")
		
			print('<td>%s</td>' %direct_targets[targ_id]["pair_id"])
			print('<td>%s</td>' %direct_targets[targ_id]["target_seq"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					print('<td>%s</td>' %direct_targets[targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")

			print('<td>%d</td>' %(direct_targets[targ_id]["coord_target"][0] + 1))
			print('<td>%d</td>' %(direct_targets[targ_id]["coord_target"][1] + 1))
		
			if is_odd(targ_id):
				print('<td>reverse</td>')
			else:
				print('<td>forward</td>')	

			print('<td>')
			print "{0:2.2f}".format(GCpercent(direct_targets[targ_id]["target_seq"]))
			print('</td>')

			print('<td>')
			print "{0:2.2f}".format(DNA_RNA_Tm(direct_targets[targ_id]["target_seq"]))
			print('</td>')

			if "exon_number" in direct_targets[targ_id].keys():
				print('<td>%s</td>' %direct_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					print('<td>')
					print "{0:3.2f}".format(calcDoenchScore(direct_targets[targ_id]["seq_forscore"]))
					print('</td>')
				else:
					print('<td>%s</td>' % "")

			print("</tr>")

		print("</table>")
		print("</div>")

		print("<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>")

		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=160 class=mono>
		"""
		sys.stdout.write('%7s' % "Pair ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
			sys.stdout.write('\t')

		sys.stdout.write('%5s' % "Start")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "End")
		sys.stdout.write('\t')
		sys.stdout.write('%7s' % "Strand")
		sys.stdout.write('\t')						
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		# initialize the boolean variable exon_info_available
		exon_info_available = False

		# print the header for Exon number 
		for targ_id in direct_targets.keys():
			if "exon_number" in direct_targets[targ_id].keys():
				exon_info_available = True
				
		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('\t')
			sys.stdout.write('%6s' % "Score")

		sys.stdout.write('\n')


		for targ_id in sorted(direct_targets, key= lambda key: direct_targets[key]["pair_id"]):
			
			# generate a row
			sys.stdout.write('%7d' %direct_targets[targ_id]["pair_id"])
			sys.stdout.write('\t')
			sys.stdout.write('%25s' %direct_targets[targ_id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					sys.stdout.write('%35s' % direct_targets[targ_id]["seq_forscore"])
					sys.stdout.write('\t')
				else:	
					sys.stdout.write('%35s' % "")
					sys.stdout.write('\t')

			sys.stdout.write('%5d' %(direct_targets[targ_id]["coord_target"][0] + 1))
			sys.stdout.write('\t')
			sys.stdout.write('%5d' %(direct_targets[targ_id]["coord_target"][1] + 1))
			sys.stdout.write('\t')

			if is_odd(targ_id):
				sys.stdout.write('%7s' % "reverse")
				sys.stdout.write('\t')	
			else:
				sys.stdout.write('%7s' % "forward")
				sys.stdout.write('\t')

			sys.stdout.write("{0:3.2f}".format(GCpercent(direct_targets[targ_id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:3.2f}".format(DNA_RNA_Tm(direct_targets[targ_id]["target_seq"])))

			if "exon_number" in direct_targets[targ_id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %direct_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in direct_targets[targ_id]:
					sys.stdout.write('\t')
					sys.stdout.write("{0:3.2f}".format(calcDoenchScore(direct_targets[targ_id]["seq_forscore"])))
				else:	
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "")

			sys.stdout.write('\n')

		print("</textarea>")
		print("</div>")

######################
# Reverse strand
######################

	if reverse_targets:	
		print """
		<h4><a href="javascript:Toggle%d();">Table of sgRNA pairs for nickase targeting in the reverse strand - expand or hide</a></h4>
		<div id="ToggleTarget%d">
		""" %(id_count, id_count)
		
		id_count += 1	

		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# all of the target data goes in here
		# all of the target data goes in here
		print """
				<table >
					<tr>
						<td>
							Pair ID
						</td>					

						<td>
							Target sequence
						</td>
		"""

		if (pam_seq == "NGG") and (length == 20):		
			print """
						<td >
							NNNN-Target-PAM-NNN
						</td>
			"""

		print """
						<td>
							Target start
						</td>
						<td>
							Target end
						</td>
		
						<td>
							Strand
						</td>
		
						<td>
							GC% of Target site
						</td>
						<td>
							Tm of sgRNA:DNA
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in reverse_targets.keys():
			if "exon_number" in reverse_targets[targ_id].keys():
				exon_info_available = True

		if exon_info_available:
			print """
						<td>
							Exon number
						</td>

			"""

		if (pam_seq == "NGG") and (length == 20):
			print """
						<td>
							Score
						</td>
			"""


		print("</tr>")

		for targ_id in sorted(reverse_targets, key= lambda key: reverse_targets[key]["pair_id"]):
		
			# generate a row
			print("<tr>")
		
			print('<td>%s</td>' %reverse_targets[targ_id]["pair_id"])
			print('<td>%s</td>' %reverse_targets[targ_id]["target_seq"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					print('<td>%s</td>' %reverse_targets[targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")

			print('<td>%d</td>' %(reverse_targets[targ_id]["coord_target"][0] + 1))
			print('<td>%d</td>' %(reverse_targets[targ_id]["coord_target"][1] + 1))
		
			if is_odd(targ_id):
				print('<td>forward</td>')
			else:
				print('<td>reverse</td>')	

			print('<td>')
			print "{0:2.2f}".format(GCpercent(reverse_targets[targ_id]["target_seq"]))
			print('</td>')

			print('<td>')
			print "{0:2.2f}".format(DNA_RNA_Tm(reverse_targets[targ_id]["target_seq"]))
			print('</td>')

			if "exon_number" in reverse_targets[targ_id].keys():
				print('<td>%s</td>' %reverse_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					print('<td>')
					print "{0:3.2f}".format(calcDoenchScore(reverse_targets[targ_id]["seq_forscore"]))
					print('</td>')
				else:
					print('<td>%s</td>' % "")

			print("</tr>")

		print("</table>")
		print("</div>")

		print("<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>")

		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=160 class=mono>
		"""
		sys.stdout.write('%7s' % "Pair ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')

		if (pam_seq == "NGG") and (length == 20):
			if "seq_forscore" in reverse_targets[targ_id]:
				sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
				sys.stdout.write('\t')

		sys.stdout.write('%5s' % "Start")
		sys.stdout.write('\t')							
		sys.stdout.write('%4s' % "End")
		sys.stdout.write('\t')
		sys.stdout.write('%7s' % "Strand")
		sys.stdout.write('\t')						
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		# initialize the boolean variable exon_info_available
		exon_info_available = False

		# print the header for Exon number 
		for targ_id in reverse_targets.keys():
			if "exon_number" in reverse_targets[targ_id].keys():
				exon_info_available = True
				
		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):
			if "seq_forscore" in reverse_targets[targ_id]:
				sys.stdout.write('\t')
				sys.stdout.write('%6s' % "Score")

		sys.stdout.write('\n')


		for targ_id in sorted(reverse_targets, key= lambda key: reverse_targets[key]["pair_id"]):
			
			# generate a row
			sys.stdout.write('%7d' %reverse_targets[targ_id]["pair_id"])
			sys.stdout.write('\t')
			sys.stdout.write('%25s' %reverse_targets[targ_id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					sys.stdout.write('%35s' % reverse_targets[targ_id]["seq_forscore"])
					sys.stdout.write('\t')
				else:	
					sys.stdout.write('%35s' % "")
					sys.stdout.write('\t')

			sys.stdout.write('%5d' %(reverse_targets[targ_id]["coord_target"][0] + 1))
			sys.stdout.write('\t')
			sys.stdout.write('%5d' %(reverse_targets[targ_id]["coord_target"][1] + 1))
			sys.stdout.write('\t')

			if is_odd(targ_id):
				sys.stdout.write('%7s' % "forward")
				sys.stdout.write('\t')	
			else:
				sys.stdout.write('%7s' % "reverse")
				sys.stdout.write('\t')

			sys.stdout.write("{0:3.2f}".format(GCpercent(reverse_targets[targ_id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:3.2f}".format(DNA_RNA_Tm(reverse_targets[targ_id]["target_seq"])))

			if "exon_number" in reverse_targets[targ_id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %reverse_targets[targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				if "seq_forscore" in reverse_targets[targ_id]:
					sys.stdout.write('\t')
					sys.stdout.write("{0:3.2f}".format(calcDoenchScore(reverse_targets[targ_id]["seq_forscore"])))
				else:	
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "")

			sys.stdout.write('\n')

		print("</textarea>")
		print("</div>")

################################
# form data variable assignment 
################################

# obtain the information from the form
form = cgi.FieldStorage(keep_blank_values=1)

dinuc = form["5dinuc"].value
length = form["target_length"].value

# error checking of the target length input
try:
	length = int(length)
except ValueError:
	badrequest_length(length)

# orientation of the PAM sequence
oriented = form["oriented"].value
oriented = int(oriented)

# an option where an user can choose whether to use the fixed most common PAM sequence
pam_value = form["PAM"].value
pam_value = int(pam_value)

if pam_value == 1:
	pam_seq = "NGG"
elif pam_value == 2:
	pam_seq = form["PAM_seq"].value
	
# type of gRNA design	
design_type = form["design_type"].value

# In case the user choses nickase design it is also necessary to collect spacer parameters
min_spacer = form["min_spacer_length"].value

# check the min_spacer variable
try:
	min_spacer = int(min_spacer)
except ValueError:
	badrequest_spacer(min_spacer)

max_spacer = form["max_spacer_length"].value

# check the max_spacer variable
try:
	max_spacer = int(max_spacer)
except ValueError:
	badrequest_spacer(max_spacer)

# test if the max_spacer is > than min_spacer
if max_spacer >= min_spacer:
	pass
else: 
	badrequest_spacer(max_spacer)

# textarea sequence input
text_input = form["input_seq"].value

# input of the file to upload
file_input = form["seqdatafile"].value

# identifier input
geneids_input = form["input_ids"].value
	
# selected species retrieval
species = form["species"].value

# retrieve the molecule type
mol_type = form["moltype"].value


#####################################
# error checking and data validation
#####################################

# first validate the sequence data entered by the user

# create a dictionary to store sequences in case of multi FASTA file
sequences = {}

# case 1: the textarea was submitted empty
if text_input == "":
	
	# check that the provided input file has any information
	# and assign its contents to the text_input
	
	if file_input != "":
		text_input = file_input
	else:
	# if there were no identifiers entered, output an error message
		if geneids_input == "":
			badrequest_seq(text_input)
			
# proceed if there is some sequence input	
if text_input != "":

	# case 2: the user entered FASTA file containing 1 or more sequences	
	if ">" in text_input:
	
		# temporary variable to store data from a split operation
		seqs = []
		seqs = text_input.split('>')
		
		if len(seqs) >= 1:
			
			if ">" != text_input[0]:
				badrequest_fasta()
		
		
		for s in seqs:
		
			if s != "":
				lines = s.split('\n')
				id = lines[0]
				sequence = ""
				
				# processing each line and removing newlines and whitespace characters from the ENDS of each line
				for line in lines[1:]:
					line.strip('\n')
					sequence += line.strip()
				
				# removing whitespace from the middle of the sequence
				sequence = ''.join(sequence.split()) 
				sequence = sequence.upper()
				
				# store the sequence in the dictionary
				sequences[id] = sequence
				
				# checking if the sequence is a DNA sequence
				if not validate(sequence):
					pass
					badrequest_seq(sequence)

	# Input contains some text which needs to be verified as DNA or something else
	if ">" not in text_input:
	# write the code to deal with the case of a single sequence
		# initializing variables
		sequence = ""
		
		# getting all of the data
		lines = text_input.split('\n')

		# processing each line and removing newlines and whitespace characters from the ENDS of each line
		for line in lines:
			line.strip('\n')
			sequence += line.strip()
			
		# removing whitespace from the middle of the sequence
		sequence = ''.join(sequence.split())
		
		# checking if the sequence is a DNA sequence
		if not validate(sequence):
			pass
			badrequest_seq(sequence)
		
		sequences["User input"] = sequence

# test the size of input sequences
for key in sequences.keys():
	if len(sequences[key])>= 50000:
		badrequest_size(key)	
	
# check and produce an error if the user did not specify their own PAM sequence
if not validate_pam(pam_seq):
	badrequest_pam(pam_seq)


if len(pam_seq) == 0:
	if pam_value == 2:
		badrequest_pam(pam_seq)

if (species == "default") and (mol_type != "refseq"):
	if mol_type != "sequence_input":
		badrequest_species(species)		

#######################
# multigene script code
#######################

# Among the molecule types available to the user to choose there is a "Sequence Input"
# option, which is completely equivalent to using one of the options for entering sequences

if mol_type == "sequence_input":
	if text_input != "":
		pass #
	else:
		badrequest_ids(geneids_input) # output an error message saying that it is not correct to use identifiers and sequences together

# first validate the gene identifier entered by the user
# the basic validation should be checking if the field was not empty
# the next step will be to check that it is possible to retrieve data from a database using this identifier, 
# which can also be part of the error checking process
		
# case 1: the sequence identifiers field contains data
if (geneids_input != "") and (text_input == ""):
	# case 2: the number of identifiers is too low	
		
	# generate a list of gene or transcript identifiers
	geneids_list = geneids_input.split("\r\n")

	# define a count variable for the number of potentially valid identifiers
	id_count = 0

	# a list for potentially valid and cleaned-up identifiers
	geneids = []

	# check how many identifiers were entered by the user
	for gene_id in geneids_list:
		
		# clean up the gene_id variable by removing all possible spaces from it
		gene_id = ''.join(gene_id.split())
		
		# check the length of this identifier before proceeding
		if len(gene_id) >= 2:
			# add the new cleaned identifier to the geneids list and increment the count of potentially valid identifiers
			geneids.append(gene_id)
			id_count += 1

	# case 3: any of the identifiers are not correct
		
	########################################
	# ADD MORE VALIDATION FOR IDENTIFIERS
	########################################
	# - pattern matching of an identifier to find out what kind of identifier was entered
	# - database query validation of an identifier and the species name

	# also, establish a database connection
	os.chdir("/home/multicri/cgi-bin/")

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
	if species != 'default':
		database = DBs_species[species]
		conn = db.connect(database)
		cur = conn.cursor()
	elif (species == "default") and (mol_type == "refseq"):
		pass
	elif (species == "default") and (mol_type == "sequence_input"):
		pass
	else:
		badrequest_species(geneids_input)


	##########################################################
	# RETRIEVAL OF DATA FROM THE DATABASE
	##########################################################

	# 1. Prepare data structures for the subsequent steps in the process, such as a basic sequences dictionary with transcript sequences,
	# but also exon structures of each transcript.

	# a dictionary containing exon sequence information
	seqs_exons = {}

	# a dictionary containing gene or transcript coordinates of the exons
	coords_exons = {}
	# gene => exonID =>dict(start => x; end => y; exon_number)
	
	# a brief algorithm for generating this dictionary
	
	# 1. retrieve information from the database
	# 2. store the location information in the dictionary
	# 3. once the dictionary is populated for this gene, sort it according to the start location of the exons
	# 4. iterate over the exons dictionary in the ascending order of "start" and assign the exon_number based on this order


	# 2. Retrieve all the possible data items from the database and verify that this has worked well.
	#    Consider raising some error if the process does not go as planned.	
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
				# generate an error	
				badrequest_seq(gene_id)
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
				sequences[orig_id] = str(sequence)
				
				# now it is necessary to obtain the exon information for this gene
				seqs_exons[orig_id] = {}
				coords_exons[orig_id] = {}				

				# extract exon sequences from the database and their coordinates in the gene
				for row in cur.execute("SELECT exonID, sequence, genestart, geneend FROM Exons WHERE geneid = ?", (geneid,)):
					# data from the current row
					exonID = row[0]
					exon_seq = row[1]
					start = row[2]
					end = row[3]
				
					if end is not None:
						end = end.rstrip()
													
					if (start is not None) and (end is not None):

						# store them in the dictionary
						seqs_exons[orig_id][exonID] = str(exon_seq)

						# store them in the dictionary
						coords_exons[orig_id][exonID] = {}	

						if start == "0":
							coords_exons[orig_id][exonID]["start"] = 0
						else:
							coords_exons[orig_id][exonID]["start"] = int(start)

						coords_exons[orig_id][exonID]["end"] = int(end)			
										

		elif mol_type == "transcripts":
			# try to select the data from Transcripts table on the assumption that the transcript identifier is of Ensembl ID format
			cur.execute("SELECT * from Transcripts WHERE transcriptID = ? COLLATE NOCASE", (gene_id,))
			row = cur.fetchone()
			
			if row is None:
				# generate an error	
				badrequest_seq(gene_id)
			else:			
				# retrieve the data from the Transcripts table
				transcriptID = row[0]
				geneid = row[1]
				sequence = row[2]
				
				# remove the newline from the sequence variable
				sequence = sequence.rstrip()
				
				# store the sequence for this identifier
				sequences[orig_id] = str(sequence)

				# now it is necessary to obtain the exon information for this transcript
				seqs_exons[orig_id] = {}
				coords_exons[orig_id] = {}	
				
				# extract exon sequences from the database and their coordinates in the gene
				for row in cur.execute("SELECT exonID, sequence, genestart, geneend FROM Exons WHERE geneid = ?", (geneid,)):
					# data from the current row
					exonID = row[0]
					exon_seq = row[1]
					start = row[2]
					end = row[3]

					if end is not None:
						end = end.rstrip()
																	
					if (start is not None) and (end is not None):

						# store them in the dictionary
						seqs_exons[orig_id][exonID] = str(exon_seq)

						# store them in the dictionary
						coords_exons[orig_id][exonID] = {}	

						if start == "0":
							coords_exons[orig_id][exonID]["start"] = 0
						else:
							coords_exons[orig_id][exonID]["start"] = int(start)

						coords_exons[orig_id][exonID]["end"] = int(end)	
			
		elif mol_type == "genes":
		# assume that the gene identifier is a gene symbol
			cur.execute("SELECT * from Genes WHERE symbol = ? COLLATE NOCASE AND species = ? COLLATE NOCASE", (gene_id, species))
			row = cur.fetchone()

			if row is None:
				# generate an error	
				badrequest_seq(gene_id)
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
				sequences[orig_id] = str(sequence)
				
				# now it is necessary to obtain the exon information for this gene
				seqs_exons[orig_id] = {}
				coords_exons[orig_id] = {}
				
				# extract exon sequences from the database and their coordinates in the gene
				for row in cur.execute("SELECT exonID, sequence, genestart, geneend FROM Exons WHERE geneid = ?", (geneid,)):
					# data from the current row
					exonID = row[0]
					exon_seq = row[1]
					start = row[2]
					end = row[3]

					if end is not None:
						end = end.rstrip()
				
					if (start is not None) and (end is not None):

						# store them in the dictionary
						seqs_exons[orig_id][exonID] = str(exon_seq)

						# store them in the dictionary
						coords_exons[orig_id][exonID] = {}	

						if start == "0":
							coords_exons[orig_id][exonID]["start"] = 0
						else:
							coords_exons[orig_id][exonID]["start"] = int(start)
						coords_exons[orig_id][exonID]["end"] = int(end)	

		elif mol_type == "refseq":
			
			######################################
			# Initial database retrieval version
			######################################

			# check that the sequence for this accession number is stored in the database
			#cur.execute("SELECT * from Refseq WHERE refseqid = ? COLLATE NOCASE", (gene_id,))
			#row = cur.fetchone()

			#if row is None:
				# generate an error	
			#	badrequest_seq(gene_id)
			#else:			
			#	# retrieve the data from the Refseq table
			#	refseqid = row[0]
			#	description = row[1] # need to decide what to do with this information item
			#	sequence = row[2]
			#	sequence = sequence.rstrip()
				
				# store the sequence for this identifier
			#	sequences[orig_id] = str(sequence)


			######################################
			# Entrez retrieval version
			######################################

			rec = Entrez.read(Entrez.esearch(db="nucleotide", term=gene_id))

			if len(rec) != 0:
				handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=rec["IdList"][0])
			
				for seq_record in SeqIO.parse(handle, "gb"):
					sequences[orig_id] = str(seq_record.seq)

				
	# close the connection to the database
	if mol_type != "refseq":
		conn.close()	

#######################
# END of multigene code
#######################
		

###############################################################
# CRISPR sgRNA target finding and preparing results for output
###############################################################

# make a regex and perform matching according to the type of gRNA design 
if design_type == "single":

	# We have all the data components to search for sgRNA targets inside 
	# the sequences we have

	#  step 1: Construct a regular expression out of the available input
	#  components. This will be an operation to be used in several scripts 
	#  and for multiple sequences so a function will be created

	sgRNA_regex = simple_regex_generate(dinuc, length, pam_seq, oriented)

	# step 2: run the constructed regular expression on all input sequences and 
	# store the results in some data structure for future use when producing the output

	# a dictionary to store the results of CRISPR sgRNA ing against all of the 
	# user input sequences
		
	
	# Plan the structure of the main webpage and how to make it clear, appealing and useful
	
	print "Content-type: text/html\n\n"

	crispr_targeter_header(sequences)

	sequences_targets = {}
	sequences_targets = regex_simple_match(sequences, sgRNA_regex, pam_seq, length)

	############################
	# output webpage generation
	############################

	# produce the output, both sequence with highlighted target sites and a list in table format
	# sgRNA_matches_output(sequences_targets)
	print """
	<p>           </p>
	<div class="form_description">
	<h3>The results of the search for CRISPR sgRNA targets in input sequences</h3>
	<hr style="height:3px;border-width:4; width: 1120px; color:gray;background-color:gray">
	</div>	
	"""

	highlight_targets_output(sequences_targets, design_type)

	crispr_targeter_footer()

elif design_type == "double":
	# We have all the data components to search for sgRNA targets inside 
	# the sequences we have

	#  step 1: Construct a regular expression out of the available input
	#  components. This will be an operation to be used in several scripts 
	#  and for multiple sequences so a function will be created

	sgRNA_regex = double_regex_generate(dinuc, length, pam_seq, min_spacer, max_spacer)

	# step 2: run the constructed regular expression on all input sequences and 
	# store the results in some data structure for future use when producing the output

	# a dictionary to store the results of CRISPR sgRNA matching against all of the 
	# user input sequences

	# Plan the structure of the main webpage and how to make it clear, appealing and useful

	print "Content-type: text/html\n\n"

	crispr_targeter_header(sequences)

	sequences_targets = {}
	sequences_targets = regex_double_match(sequences, sgRNA_regex, pam_seq, length)

	############################
	# output webpage generation
	############################

	# produce the output, both sequence with highlighted target sites and a list in table format
	# sgRNA_matches_output(sequences_targets)
	print """
	<p>           </p>
	<div class="form_description">
	<h3>The results of the search for CRISPR sgRNA targets in input sequences</h3>
	<hr style="height:3px;border-width:4; width: 1050px; color:gray;background-color:gray">
	</div>	
	"""


	highlight_targets_output(sequences_targets, design_type)

	crispr_targeter_footer()