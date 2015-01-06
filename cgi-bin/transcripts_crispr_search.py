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
		line = " "*30
		line = line +  "\t"

		sys.stdout.write("%s" %line)

		# shift the markers by 1 position
		shift = " "
		sys.stdout.write("%s" %shift)

		# initialize the most recent printed position
		last = prev_slice	

		for site in segment:
			# define a coordinate of site in the segment
			coord = site - last
		
			num_spaces = coord - 1 
			last = site

			line = " "*num_spaces
			sys.stdout.write("%s" %line)
			sys.stdout.write("<strong><font style='background-color: green; color:white'>&gt;</font></strong>")
		
		num_spaces = curr_slice - last
		line = " "*num_spaces
		sys.stdout.write("%s" %line)
		sys.stdout.write("\n")

def print_markers_uniq(starts, curr_slice):

	# define the lower boundary for starts
	prev_slice = curr_slice - 100
	
	# store start sites in the current segment
	segment = []

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
		line = " "*8

		sys.stdout.write("%s" %line)

		# initialize the most recent printed position
		last = prev_slice	

		for site in segment:
			# define a coordinate of site in the segment
			coord = site - last
		
			num_spaces = coord - 1 
			last = site

			line = " "*num_spaces
			sys.stdout.write("%s" %line)
			sys.stdout.write("<strong><font style='background-color: green; color:white'>&gt;</font></strong>")
		
		num_spaces = curr_slice - last
		line = " "*num_spaces
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
	seq_count = 3*len(sequences) + 6  # maximum number of necessary sections which will have to be hidden and expanded

	print """
<!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Find common and unique CRISPR guide RNA targets for different transcripts of one gene</title>


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
     <p> CRISPR MultiTargeter was developed by Sergey Prykhozhij at the IWK Health Centre and Dalhousie University. The design of the logo and the buttons was by Vinothkumar Rajan. The latest update was on the 5th of January 2015.</p></div>
</div>
</div>
</body>
</html>
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
	print("<p>The gene identifier <span style='color: #ff0000;'>'%s'</span> you entered is not valid. Please try again later.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()
	
def badrequest_input(gene):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>The gene %s has only one transcript. The analysis in this script is not applicable to this gene.</p>" %gene)
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
	
def badrequest_species(bad):
	#Display an error message
	error_header()
	print("<h2>Error!</h2>")
	print("<p>It is not possible to process sequence identifiers <span style='color: #ff0000;'>'%s'</span>  when 'Species' option is not chosen. Please choose the correct species.</p>" %bad)
	print("</body></html>")
	# Get out of here:
	return sys.exit()


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



def block_output(sequences_targets, aln_block):
	j = 0
	
	for id in sorted(sequences_targets.keys()): # assume that the lines in the alignment block are ordered the same way as sorted keys
												# in a dictionary
		if id != "consensus":
			id_line = str(id)
			
			# adjust the length of all identifiers to the same value to keep alignment 
			if len(id_line) > 30:
				id_line = id_line[:30]
			elif len(id_line) == 30:
				pass
			elif len(id_line) < 30:
				id_line = id_line + (30-len(id_line))*' '
			
			
			sys.stdout.write(id_line)
			sys.stdout.write("\t")
			sys.stdout.write("".join(aln_block[j][1:]))
			sys.stdout.write("\n")
			j += 1				
	cons_line = str(aln_block[-1][0]) + 21*' '
	sys.stdout.write(cons_line)
	sys.stdout.write("\t")
	
	cons_line = "".join(aln_block[-1][1:])
	cons_line = cons_line.replace('X', ' ')
	sys.stdout.write(cons_line)
	sys.stdout.write("\n")
			
	# newline to separate alignment blocks
	sys.stdout.write("\n")				
		

def initialize_block(sequences_targets):
	# initialize a 2-dimensional list to store blocks of highlighted alignment pieces for output
	aln_block = []
	
	length = len(sequences_targets.keys())
	i = 0
	
	while i < length:
		aln_block.append([])
		i += 1
	
	# populate the very first cell of each of the inside list with the sequence identifier
	# then use the same order for alignment block output
	i = 0
	
	for id in sorted(sequences_targets.keys()):
		if id != "consensus":
			aln_block[i].append(id)
			i += 1
		
	aln_block[i].append("consensus")
	
	return(aln_block)

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
	
def highlight_targets_output(sequences_targets_forw, sequences_targets_rev): 
 	
	id_count = 1
	
	# print basic headers
	print( """
	<p>           </p>
	<div class="form_description">
	<h3>Common and unique CRISPR sgRNA targets in different transcripts of %s gene</h3>
	<hr style="height:3px;border-width:4; width: 1120px; color:gray;background-color:gray">
	</div>	
	""" %entered_geneid)

	#####################
	# Ensembl link output
	#####################
		

	if species != "Oriza sativa japonica" :
		(genus, sp) = species.split()


	animals = ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Gallus gallus", "Xenopus tropicalis", "Danio rerio", "Oryzias latipes", "Drosophila melanogaster", "Caenorhabditis elegans"]
	plants = ["Arabidopsis thaliana", "Oriza sativa japonica", "Zea mays"]

	

	# deal with different species
	if species in animals:
		print("<h4><a style='color: #0000FF;' href='http://useast.ensembl.org/%s_%s/Gene/Summary?db=core;g=%s'>%s Ensembl Link</a></h4>" %(genus, sp, entered_geneid, entered_geneid))
	elif species in plants:
		print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(entered_geneid, entered_geneid))
				
	elif species == "different":
				
		if "AT" in seq:
			print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(entered_geneid, entered_geneid))
		elif "OS" in seq:
			print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(entered_geneid, entered_geneid))
		elif "GRMZM" in seq:
			print("<h4><a style='color: #0000FF;' href='http://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit'>%s EnsemblPlants Link</a></h4>" %(entered_geneid, entered_geneid))
		else:
			print("<h4><a style='color: #0000FF;' href='http://useast.ensembl.org/Gene/Summary?db=core;g=%s'>%s Ensembl Link</a></h4>" %(entered_geneid, entered_geneid))
	#####################
	# END of OUTPUT
	#####################   

	print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px; vertical-align: middle;'><h3>Visual View of Common sgRNA Targets</h3></div>")
		
	################################
	# Forward strand
	################################
		
	print("<div style='background-color:#0000DD; color: white; width: 500px; height: 25px; margin:10px;'><p class='normalp'><strong>Targets in the forward strand:</strong></p></div>")
		
	# print the <pre> element which will contain all of the sequences and highlighting
	print """
	<a href="javascript:Toggle%d();"><h4><strong>Please click to expand or hide the alignment</strong></h4></a>
	<div id="ToggleTarget%d">
	""" % (id_count, id_count)
	
	print("<p style='font-family: Times New Roman, Times, serif; font-size: 18px;'><strong>Guide RNA targets have been highlighted. Long stretches of highlighted sequence are due to target overlap.<br />Please see tables for individual sequences.</strong></p>")
	print("<pre style='margin:10px; padding:5px; line-height: 1.6em; font-size: 8pt'>")
		
	# store the consensus sequence from the input dictionary in a variable
	sequence = sequences_targets_forw["consensus"]["align_seq"]
	
	# initialize the dictionary to use as an alias of sgRNA target matches
	direct_targets = {}
	direct_targets = sequences_targets_forw["consensus"]["aligned_matches"]
	
	# do a test for the possibility that no targets of sgRNA have been identified
	if not sequences_targets_forw["consensus"]["aligned_matches"]:
		print("This multiple sequence alignment does not contain any common sgRNA targets in the sense orientation.")
	
	# iterate over all the targets of CRISPR sgRNAs in the consensus sequence
	# 
	# Use the following styles for the sgRNA targets and adjacent PAM sequences:
	# <strong><font style="BACKGROUND-COLOR: #0000CC; color: white">sgRNA target</font></strong>

	######################################################
	#  Start site marker labeling
	######################################################
		
	# generate a list of target start sites
	start_coords = []

	# looping over sgRNA targets:
	for targ_id in sorted(direct_targets, key = lambda key: direct_targets[key]["coord_target"][0]):
		begin = direct_targets[targ_id]["coord_target"][0]
		start_coords.append(begin)

	######################################################
	# End of start coordinates list of coordinates
	######################################################

	
	# initialize the positions to use while producing output
	curr_pos = 0
	curr_slice = 100

	# print the markers if there are any sites within this segment
	print_markers(start_coords, curr_slice)

	# initialize a 2-dimensional list to store blocks of highlighted alignment pieces for output
	aln_block = initialize_block(sequences_targets_forw)
	
	# looping over sgRNA targets:	 
	for targ_id in direct_targets.keys():
						
		# Case 1: the start position of the target is smaller than the current slice
		if direct_targets[targ_id]["coord_target"][0] < curr_slice:
			# print the sequence and the target and prepare the ground for the remaining sequence
			# first print the sequence leading up to the sgRNA match
			
			first = direct_targets[targ_id]["coord_target"][0]
			# check for an overlap between the previous target and the new one
			if first < curr_pos:
				first = curr_pos
			else:				
			
				########################
				# Alignment block update
				########################

				j = 0
		
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: first])
						
						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append("%s" %sequence[curr_pos: first])

				########################
				# Alignment block update
				########################
			
				curr_pos = first
				
			# check where the target ends
			if direct_targets[targ_id]["coord_target"][1] < curr_slice:
				last = direct_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" % sequences_targets_forw[id]["align_seq"][curr_pos: last])
						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################
				
				curr_pos = last
				
			elif direct_targets[targ_id]["coord_target"][1] == curr_slice:
				last = direct_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])
						
						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################

				
				
				########################
				# BLOCK OUTPUT
				########################
				
				# alignment block output
				block_output(sequences_targets_forw, aln_block)
				
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_forw)
				
				########################
				# BLOCK OUTPUT
				########################					

				# increment the relevant positions
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
					
			elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
				last = curr_slice
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])

						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################				
				
				
				
				########################
				# BLOCK OUTPUT
				########################
				
				# alignment block output
				block_output(sequences_targets_forw, aln_block)
				
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_forw)
				
				########################
				# BLOCK OUTPUT
				########################					

					
				curr_pos = curr_slice
				curr_slice += 100
				last = direct_targets[targ_id]["coord_target"][1]
					
 
				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
					
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])

						j += 1
			
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################
				
				curr_pos = last
				
		# Case 2: the start position of the target is equal or larger than the current slice  
		if direct_targets[targ_id]["coord_target"][0] >= curr_slice:
			
			first = direct_targets[targ_id]["coord_target"][0]

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			
			for id in sorted(sequences_targets_forw.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: curr_slice])

					j += 1
			
			aln_block[j].append("%s" %sequence[curr_pos: curr_slice])
				
			########################
			# Alignment block update
			########################			
			
			########################
			# BLOCK OUTPUT
			########################
				
			# alignment block output
			block_output(sequences_targets_forw, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_forw)
				
			########################
			# BLOCK OUTPUT
			########################
			
			curr_pos = curr_slice
			curr_slice += 100

			# print the markers if there are any sites within this segment
			print_markers(start_coords, curr_slice)
				
			# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
			while first >= curr_slice:

				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("%s"  %sequences_targets_forw[id]["align_seq"][curr_pos: curr_slice])

						j += 1
				
				aln_block[j].append("%s" %sequence[curr_pos: curr_slice])
					
				########################
				# Alignment block update
				########################			
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_forw, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_forw)
					
				########################
				# BLOCK OUTPUT
				########################
				
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_forw.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: first])

					j += 1
			
			aln_block[j].append("%s" %sequence[curr_pos: first])
				
			########################
			# Alignment block update
			########################				
				
			curr_pos = first
				
			# check where the target ends
			if direct_targets[targ_id]["coord_target"][1] < curr_slice:
				last = direct_targets[targ_id]["coord_target"][1]
				
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				curr_pos = last
					
			elif direct_targets[targ_id]["coord_target"][1] == curr_slice:
				last = direct_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_forw, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_forw)
					
				########################
				# BLOCK OUTPUT
				########################
					
				# increment the relevant positions
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
					
			elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
				last = curr_slice
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_forw, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_forw)
					
				########################
				# BLOCK OUTPUT
				########################
				
				curr_pos = curr_slice
				curr_slice += 100
				last = direct_targets[targ_id]["coord_target"][1]

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_forw.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequences_targets_forw[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################    # mark with highlighted font
				curr_pos = last
		
	
	# finishing the end of the sequence
		 
	if curr_slice >= len(sequence):
		# output the rest of the sequence

		########################
		# Alignment block update
		########################
		# now update the aln_block list to store this sequence segment for all of the aligned sequences
		# include highlighting of the sequences
		j = 0
		for id in sorted(sequences_targets_forw.keys()):
			if id != "consensus":
				aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: ])
				
				j += 1
				
		aln_block[-1].append( "%s" %sequence[curr_pos: ])
					
		#######################
		# Alignment block update
		#######################				
				
		########################
		# BLOCK OUTPUT
		########################
					
		# alignment block output
		block_output(sequences_targets_forw, aln_block)
			
		# re-initialize the aln_block variable
		aln_block = initialize_block(sequences_targets_forw)
					
		########################
		# BLOCK OUTPUT
		########################

	else:
		# finish the current line

		########################
		# Alignment block update
		########################
		# now update the aln_block list to store this sequence segment for all of the aligned sequences
		# include highlighting of the sequences
		j = 0
		for id in sorted(sequences_targets_forw.keys()):
			if id != "consensus":
				aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: curr_slice])
				
				j += 1
				
		aln_block[-1].append( "%s" %sequence[curr_pos: curr_slice])
					
		#######################
		# Alignment block update
		#######################				
				
		########################
		# BLOCK OUTPUT
		########################
					
		# alignment block output
		block_output(sequences_targets_forw, aln_block)
			
		# re-initialize the aln_block variable
		aln_block = initialize_block(sequences_targets_forw)
					
		########################
		# BLOCK OUTPUT
		########################		
		
		# update the variables
		curr_pos = curr_slice
		curr_slice += 100

		# print the markers if there are any sites within this segment
		print_markers(start_coords, curr_slice)
			
		# go over any other remaining lines and print them
		while curr_slice < len(sequence):
			# finish the current line
			
			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_forw.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: curr_slice])
					
					j += 1
					
			aln_block[-1].append( "%s" %sequence[curr_pos: curr_slice])
						
			#######################
			# Alignment block update
			#######################				
					
			########################
			# BLOCK OUTPUT
			########################
						
			# alignment block output
			block_output(sequences_targets_forw, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_forw)
						
			########################
			# BLOCK OUTPUT
			########################				
			
			# update the variables
			curr_pos = curr_slice
			curr_slice += 100

			# print the markers if there are any sites within this segment
			print_markers(start_coords, curr_slice)
			
		if curr_slice >= len(sequence):
			# output the rest of the sequence

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_forw.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_forw[id]["align_seq"][curr_pos: ])
					
					j += 1
					
			aln_block[-1].append( "%s" %sequence[curr_pos: ])
						
			#######################
			# Alignment block update
			#######################				
					
			########################
			# BLOCK OUTPUT
			########################
						
			# alignment block output
			block_output(sequences_targets_forw, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_forw)
						
			########################
			# BLOCK OUTPUT
			########################			
	
	print("</pre>")
	print("</div>")
	
	
	# ##################
	# # Reverse strand
	# ##################

	print("<div style='background-color: #CF5300; color: white; width: 500px; height: 25px; margin:10px;'><p class='normalp'><strong>Targets in the reverse strand:</strong></p></div>")
	
	id_count += 1
	
	print """
	<a href="javascript:Toggle%d();"><h4><strong>Please click to expand or hide the alignment</strong></h4></a>
	<div id="ToggleTarget%d">
	""" % (id_count, id_count)
	
	print("<p style='font-family: Times New Roman, Times, serif; font-size: 18px;'><strong>Guide RNA targets have been highlighted. Long stretches of highlighted sequence are due to target overlap.<br />Please see tables for individual sequences.</strong></p>")
	
	# print the <pre> element which will contain all of the sequences and highlighting
	print("<pre style='margin:10px; padding:5px; line-height: 1.6em; font-size: 8pt'>")
		
	# store the consensus sequence from the input dictionary in a variable
	sequence = sequences_targets_rev["consensus"]["align_seq"]
	 
	# initialize the dictionary to use as an alias of sgRNA target matches
	reverse_targets = {}
	reverse_targets = sequences_targets_rev["consensus"]["aligned_matches"]
	
	# do a test for the possibility that no targets of sgRNA have been identified
	if not sequences_targets_rev["consensus"]["aligned_matches"]:
		print("This multiple sequence alignment does not contain any common sgRNA targets in the sense orientation.")
	 
	
	# iterate over all the targets of CRISPR sgRNAs in the consensus sequence
	# 
	# Use the following styles for the sgRNA targets and adjacent PAM sequences:
	# <strong><font style="BACKGROUND-COLOR: #CF5300; color: white">sgRNA target</font></strong>

	######################################################
	#  Start site marker labeling
	######################################################

	# generate a list of target start sites
	start_coords = []

	# looping over sgRNA targets:
	for targ_id in sorted(reverse_targets, key = lambda key: reverse_targets[key]["coord_target"][0]):
		begin = reverse_targets[targ_id]["coord_target"][0]
		start_coords.append(begin)

	######################################################
	# End of start coordinates list of coordinates
	######################################################
	
	# initialize the positions to use while producing output
	curr_pos = 0
	curr_slice = 100

	# print the markers if there are any sites within this segment
	print_markers(start_coords, curr_slice)

	# initialize a 2-dimensional list to store blocks of highlighted alignment pieces for output
	aln_block = initialize_block(sequences_targets_rev)
	
	# looping over sgRNA targets:	 
	for targ_id in reverse_targets.keys():
						
		# Case 1: the start position of the target is smaller than the current slice
		if reverse_targets[targ_id]["coord_target"][0] < curr_slice:
			# print the sequence and the target and prepare the ground for the remaining sequence
			# first print the sequence leading up to the sgRNA match
			
			first = reverse_targets[targ_id]["coord_target"][0]

			# check for an overlap between the previous target and the new one
			if first < curr_pos:
				first = curr_pos
			else:	
			
				########################
				# Alignment block update
				########################

				j = 0
		
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: first])
						
						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append("%s" %sequence[curr_pos: first])

				########################
				# Alignment block update
				########################
			
				curr_pos = first
				
			# check where the target ends
			if reverse_targets[targ_id]["coord_target"][1] < curr_slice:
				last = reverse_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])

						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################
				
				curr_pos = last
				
			elif reverse_targets[targ_id]["coord_target"][1] == curr_slice:
				last = reverse_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])

						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################

				
				
				########################
				# BLOCK OUTPUT
				########################
				
				# alignment block output
				block_output(sequences_targets_rev, aln_block)
				
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_rev)
				
				########################
				# BLOCK OUTPUT
				########################					

				# increment the relevant positions
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
				
			elif reverse_targets[targ_id]["coord_target"][1] > curr_slice:
				last = curr_slice
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						# update the sequence item
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])

						# move to the next item in the list
						j += 1
				# update the consensus
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################				
				
				
				
				########################
				# BLOCK OUTPUT
				########################
				
				# alignment block output
				block_output(sequences_targets_rev, aln_block)
				
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_rev)
				
				########################
				# BLOCK OUTPUT
				########################					

					
				curr_pos = curr_slice
				curr_slice += 100
				last = reverse_targets[targ_id]["coord_target"][1]

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)					
 					
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
	
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])

						j += 1
			
				aln_block[j].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
				
				########################
				# Alignment block update
				########################
				
				curr_pos = last
				
		# Case 2: the start position of the target is equal or larger than the current slice  
		if reverse_targets[targ_id]["coord_target"][0] >= curr_slice:
			
			first = reverse_targets[targ_id]["coord_target"][0]

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0	
			
			for id in sorted(sequences_targets_rev.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: curr_slice])
					j += 1
			
			aln_block[j].append("%s" %sequence[curr_pos: curr_slice])
				
			########################
			# Alignment block update
			########################			
			
			########################
			# BLOCK OUTPUT
			########################
				
			# alignment block output
			block_output(sequences_targets_rev, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_rev)
				
			########################
			# BLOCK OUTPUT
			########################
			
			curr_pos = curr_slice
			curr_slice += 100

			# print the markers if there are any sites within this segment
			print_markers(start_coords, curr_slice)
				
			# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
			while first >= curr_slice:

				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("%s"  %sequences_targets_rev[id]["align_seq"][curr_pos: curr_slice])

						j += 1
				
				aln_block[j].append("%s" %sequence[curr_pos: curr_slice])
					
				########################
				# Alignment block update
				########################			
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_rev, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_rev)
					
				########################
				# BLOCK OUTPUT
				########################
				
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_rev.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: first])

					j += 1
			
			aln_block[j].append("%s" %sequence[curr_pos: first])
				
			########################
			# Alignment block update
			########################				
				
			curr_pos = first
				
			# check where the target ends
			if reverse_targets[targ_id]["coord_target"][1] < curr_slice:
				last = reverse_targets[targ_id]["coord_target"][1]
				
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				curr_pos = last
					
			elif reverse_targets[targ_id]["coord_target"][1] == curr_slice:
				last = reverse_targets[targ_id]["coord_target"][1]
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_rev, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_rev)
					
				########################
				# BLOCK OUTPUT
				########################
					
				# increment the relevant positions
				curr_pos = curr_slice
				curr_slice += 100

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
					
			elif reverse_targets[targ_id]["coord_target"][1] > curr_slice:
				last = curr_slice
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################				
				
				########################
				# BLOCK OUTPUT
				########################
					
				# alignment block output
				block_output(sequences_targets_rev, aln_block)
					
				# re-initialize the aln_block variable
				aln_block = initialize_block(sequences_targets_rev)
					
				########################
				# BLOCK OUTPUT
				########################
				
				curr_pos = curr_slice
				curr_slice += 100
				last = reverse_targets[targ_id]["coord_target"][1]

				# print the markers if there are any sites within this segment
				print_markers(start_coords, curr_slice)
				
				########################
				# Alignment block update
				########################
				# now update the aln_block list to store this sequence segment for all of the aligned sequences
				# include highlighting of the sequences
				j = 0
				for id in sorted(sequences_targets_rev.keys()):
					if id != "consensus":
						aln_block[j].append("<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequences_targets_rev[id]["align_seq"][curr_pos: last])
				
						j += 1
				
				aln_block[-1].append( "<strong><font style='BACKGROUND-COLOR: #CF5300; color: white'>%s</font></strong>" %sequence[curr_pos: last])
					
				########################
				# Alignment block update
				########################    # mark with highlighted font
				curr_pos = last
		
	
	# finishing the end of the sequence
		 
	if curr_slice >= len(sequence):
		# output the rest of the sequence

		########################
		# Alignment block update
		########################
		# now update the aln_block list to store this sequence segment for all of the aligned sequences
		# include highlighting of the sequences
		j = 0
		for id in sorted(sequences_targets_rev.keys()):
			if id != "consensus":
				aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: ])
				
				j += 1
				
		aln_block[-1].append( "%s" %sequence[curr_pos: ])
					
		#######################
		# Alignment block update
		#######################				
				
		########################
		# BLOCK OUTPUT
		########################
					
		# alignment block output
		block_output(sequences_targets_rev, aln_block)
			
		# re-initialize the aln_block variable
		aln_block = initialize_block(sequences_targets_rev)
					
		########################
		# BLOCK OUTPUT
		########################

	else:
		# finish the current line

		########################
		# Alignment block update
		########################
		# now update the aln_block list to store this sequence segment for all of the aligned sequences
		# include highlighting of the sequences
		j = 0
		for id in sorted(sequences_targets_rev.keys()):
			if id != "consensus":
				aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: curr_slice])
				
				j += 1
				
		aln_block[-1].append( "%s" %sequence[curr_pos: curr_slice])
					
		#######################
		# Alignment block update
		#######################				
				
		########################
		# BLOCK OUTPUT
		########################
					
		# alignment block output
		block_output(sequences_targets_rev, aln_block)
			
		# re-initialize the aln_block variable
		aln_block = initialize_block(sequences_targets_rev)
					
		########################
		# BLOCK OUTPUT
		########################		
		
		# update the variables
		curr_pos = curr_slice
		curr_slice += 100

		# print the markers if there are any sites within this segment
		print_markers(start_coords, curr_slice)
	
		# go over any other remaining lines and print them
		while curr_slice < len(sequence):
			# finish the current line
			
			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_rev.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: curr_slice])
					
					j += 1
					
			aln_block[-1].append( "%s" %sequence[curr_pos: curr_slice])
						
			#######################
			# Alignment block update
			#######################				
					
			########################
			# BLOCK OUTPUT
			########################
						
			# alignment block output
			block_output(sequences_targets_rev, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_rev)
						
			########################
			# BLOCK OUTPUT
			########################				
			
			# update the variables
			curr_pos = curr_slice
			curr_slice += 100

			# print the markers if there are any sites within this segment
			print_markers(start_coords, curr_slice)
			
		if curr_slice >= len(sequence):
			# output the rest of the sequence

			########################
			# Alignment block update
			########################
			# now update the aln_block list to store this sequence segment for all of the aligned sequences
			# include highlighting of the sequences
			j = 0
			for id in sorted(sequences_targets_rev.keys()):
				if id != "consensus":
					aln_block[j].append("%s" %sequences_targets_rev[id]["align_seq"][curr_pos: ])
					
					j += 1
					
			aln_block[-1].append( "%s" %sequence[curr_pos: ])
						
			#######################
			# Alignment block update
			#######################				
					
			########################
			# BLOCK OUTPUT
			########################
						
			# alignment block output
			block_output(sequences_targets_rev, aln_block)
				
			# re-initialize the aln_block variable
			aln_block = initialize_block(sequences_targets_rev)
						
			########################
			# BLOCK OUTPUT
			########################			
	
	print("</pre>")
	print("</div>")
	
	##########################################################################
	# Print all of the targets in both orientations of the consensus sequence
	# followed by the corresponding output for the individual sequences
	##########################################################################
	
	targets_table_output(sequences_targets_forw, sequences_targets_rev)
	
	# The following is an excerpt from the main part of the script to use as a reference for the data structures
	# of the sequences and potential guide RNA matches

	# # populate the nested dictionaries with the data
		
		# # Consensus sequence first
		# sequences_matches_forward["consensus"]["aligned_matches"][match_count]["target_seq"] = match.group(1)
		# sequences_matches_forward["consensus"]["aligned_matches"][match_count]["coord_target"] = [match.start(), match.end() - len(pam_seq)]
		# sequences_matches_forward["consensus"]["aligned_matches"][match_count]["coord_match"] = [match.start(), match.end()]
		
		# # Next, iterate over all of the aligned sequences
		# for record in align:
			# seq = str(record.seq)
			# sequences_matches_forward[record.id]["aligned_matches"][match_count]["target_seq"] = seq[match.start(): match.end() - len(pam_seq)]
			# sequences_matches_forward[record.id]["aligned_matches"][match_count]["coord_target"] = [match.start(), match.end() - len(pam_seq)]
			# sequences_matches_forward[record.id]["aligned_matches"][match_count]["coord_match"] = [match.start(), match.end()]			

			# # The coordinates for matches in the original sequences are only different by the number of gaps in a part of that sequence
			# # just before the match. Therefore, it should be possible to convert the data between the two parts of the dictionary by a simple 
			# # subtraction
			
			# gaps = int(seq[:match.start()].count('-'))
			
			# sequences_matches_forward[record.id]["orig_matches"][match_count]["target_seq"] = seq[match.start() - gaps: match.end() - len(pam_seq)-gaps]
			# sequences_matches_forward[record.id]["orig_matches"][match_count]["coord_target"] = [match.start() - gaps, match.end() - gaps - len(pam_seq)]
			# sequences_matches_forward[record.id]["orig_matches"][match_count]["coord_match"] = [match.start()- gaps, match.end() - gaps]	

def targets_table_output(sequences_targets_forw, sequences_targets_rev):
	# print colored headers
	print """
	<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Table View of common sgRNA Targets. Click on the links to view.</h3></div>
	"""
	
	id_count = 3 # a variable to assign particular tables to certain style and javascript function id numbers
	
	##################################################################
	# individual sequences tables
	##################################################################
	
	if len(sequences_targets_forw["consensus"]["aligned_matches"]) == 0:
		# no need to print another message saying that there are no targets available for output
		# since this has already been done above
		pass
	else:	
		print """
		<h4><a href="javascript:Toggle%d();">Common sgRNA targets in the forward strand of different transcripts - expand or hide</a></h4>
		<div id="ToggleTarget%d">
		""" %(id_count, id_count)
		
		id_count += 1
		
		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# header for the forward strand consensus sequence target data
		print """
				<table >
					<tr>						
						<td>
							Target ID number
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
							GC% of sgRNA
						</td>

						<td>
							Tm of sgRNA:DNA duplex
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in sorted(sequences_targets_forw["consensus"]["aligned_matches"].keys()):
			if "exon_number" in sequences_targets_forw["consensus"]["aligned_matches"][targ_id].keys():
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

		print("<tr>")


		print("</tr>")
		
		# target data output for this table goes here
		for targ_id in sorted(sequences_targets_forw["consensus"]["aligned_matches"].keys()):
			
			# generate a row
			print("<tr>")
			print('<td>%d</td>' %targ_id)
			print('<td>%s</td>' %sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"])
			print('<td>')
			print "{0:2.2f}".format(GCpercent(sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"]))
			print('</td>')

			if (pam_seq == "NGG") and (length == 20):
				
				for key in sequences_targets_forw.keys():
					seq_key = key

				if "seq_forscore" in sequences_targets_forw[seq_key]["aligned_matches"][targ_id]:
					print('<td>%s</td>' %sequences_targets_forw[seq_key]["aligned_matches"][targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")

			print('<td>')
			print "{0:2.2f}".format(DNA_RNA_Tm(sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"]))
			print('</td>')
			
			if "exon_number" in sequences_targets_forw["consensus"]["aligned_matches"][targ_id].keys():
				print('<td>%s</td>' %sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				# define a data structure to store scores for targets in individual sequences
				scores = []

				# iterate over all sequences in the alignment
				for seq_key in sequences_targets_forw.keys():
					if seq_key != "consensus":

						# calculate the score
						if "seq_forscore" in sequences_targets_forw[seq_key]["aligned_matches"][targ_id]:
							scores.append(calcDoenchScore(sequences_targets_forw[seq_key]["aligned_matches"][targ_id]["seq_forscore"]))

				if len(scores) > 0:
					print('<td>')
					print "{0:2.2f}".format(sum(scores)/len(scores))
					print('</td>')
				else:
					print('<td>%s</td>' % "")

			print("</tr>")
			
		print """	
		</table>
		</div>
		<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>
		"""

		#output the textarea with the same information as above	
		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=175 class=mono>
		"""

		sys.stdout.write('%9s' % "Target ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')				

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
			sys.stdout.write('\t')
			
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('\t')
			sys.stdout.write('%5s' % "Score")


		sys.stdout.write('\n')


		for targ_id in sorted(sequences_targets_forw["consensus"]["aligned_matches"].keys()):
			# generate a row
			sys.stdout.write('%9d' %targ_id)
			sys.stdout.write('\t')
			sys.stdout.write('%25s' %sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):
				
				for key in sequences_targets_forw.keys():
					seq_key = key

				if (pam_seq == "NGG") and (length == 20):
					if "seq_forscore" in sequences_targets_forw[seq_key]["aligned_matches"][targ_id]:
						sys.stdout.write('%35s' % sequences_targets_forw[seq_key]["aligned_matches"][targ_id]["seq_forscore"])
						sys.stdout.write('\t')
					else:
						sys.stdout.write('%35s' % "")
						sys.stdout.write('\t')


			sys.stdout.write("{0:2.2f}".format(GCpercent(sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["target_seq"])))
			
			if "exon_number" in sequences_targets_forw["consensus"]["aligned_matches"][targ_id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %sequences_targets_forw["consensus"]["aligned_matches"][targ_id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				# define a data structure to store scores for targets in individual sequences
				scores = []

				# iterate over all sequences in the alignment
				for seq_key in sequences_targets_forw.keys():
					if seq_key != "consensus":

						# calculate the score
						if "seq_forscore" in sequences_targets_forw[seq_key]["aligned_matches"][targ_id]:
							scores.append(calcDoenchScore(sequences_targets_forw[seq_key]["aligned_matches"][targ_id]["seq_forscore"]))

				sys.stdout.write('\t')

				if len(scores) > 0:
					sys.stdout.write("{0:3.2f}".format(sum(scores)/len(scores)))
				else:
					sys.stdout.write("")

			sys.stdout.write('\n')

		print("</textarea>")		
		print("</div>")
		print("<br />")

	if len(sequences_targets_rev["consensus"]["aligned_matches"]) == 0:
		# no need to print another message saying that there are no targets available for output
		# since this has already been done above
		pass
	else:	
		print """
		<h4><a href="javascript:Toggle%d();">Common sgRNA targets in the reverse strand of different transcripts - expand or hide</a></h4>
		<div id="ToggleTarget%d">
		""" %(id_count, id_count)
		
		id_count += 1
		
		# print the <div> element defining the pretty table style
		print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
		# header for the reverse strand consensus sequence target data
		print """
				<table >
					<tr>
						
						<td>
							Target ID number
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
							GC% of sgRNA
						</td>

						<td>
							Tm of sgRNA:DNA duplex
						</td>
		"""

		# Now it is necessary to make this last column of the table conditional
		# upon the availability of exon location info for this sequence item
		
		# initialize the boolean variable exon_info_available
		exon_info_available = False

		for targ_id in sorted(sequences_targets_rev["consensus"]["aligned_matches"].keys()):
			if "exon_number" in sequences_targets_rev["consensus"]["aligned_matches"][targ_id].keys():
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

		
		# target data output for this table goes here
		for id in sorted(sequences_targets_rev["consensus"]["aligned_matches"].keys()):
			# generate a row
			print("<tr>")
			print('<td>%d</td>' %id)
			print('<td>%s</td>' %sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"])


			if (pam_seq == "NGG") and (length == 20):
				
				for key in sequences_targets_rev.keys():
					seq_key = key

				if "seq_forscore" in sequences_targets_rev[seq_key]["aligned_matches"][targ_id]:
					print('<td>%s</td>' %sequences_targets_rev[seq_key]["aligned_matches"][targ_id]["seq_forscore"])
				else:
					print('<td>%s</td>' % "")

			print('<td>')
			print "{0:2.2f}".format(GCpercent(sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"]))
			print('</td>')

			print('<td>')
			print "{0:2.2f}".format(DNA_RNA_Tm(sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"]))
			print('</td>')

			if "exon_number" in sequences_targets_rev["consensus"]["aligned_matches"][id].keys():
				print('<td>%s</td>' %sequences_targets_rev["consensus"]["aligned_matches"][id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				# define a data structure to store scores for targets in individual sequences
				scores = []

				# iterate over all sequences in the alignment
				for seq_key in sequences_targets_rev.keys():
					if seq_key != "consensus":

						# calculate the score
						if "seq_forscore" in sequences_targets_rev[seq_key]["aligned_matches"][targ_id]:
							scores.append(calcDoenchScore(sequences_targets_rev[seq_key]["aligned_matches"][targ_id]["seq_forscore"]))

				if len(scores) > 0:
					print('<td>')
					print "{0:3.2f}".format(sum(scores)/len(scores))
					print('</td>')
				else:
					print('<td>%s</td>' % "")

			print("</tr>")
		
		print """	
		</table>
		</div>
		<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>
		"""

		#output the textarea with the same information as above	
		print """
		<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
		(<i>e.g.</i> Excel) or text editors.</p>
		<textarea rows=10 cols=175 class=mono>
		"""
		sys.stdout.write('%9s' % "Target ID")
		sys.stdout.write('\t')
		sys.stdout.write('%25s' % "Target sequence")
		sys.stdout.write('\t')

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('%35s' % "NNNN-Target-PAM-NNN")
			sys.stdout.write('\t')
						
		sys.stdout.write('%4s' % "GC %")
		sys.stdout.write('\t')
		sys.stdout.write('%4s' % "Tm")

		if exon_info_available:
			sys.stdout.write('\t')
			sys.stdout.write('%20s' % "Exon")

		if (pam_seq == "NGG") and (length == 20):
			sys.stdout.write('\t')
			sys.stdout.write('%5s' % "Score")

		sys.stdout.write('\n')




		# target data output for this table goes here
		for id in sorted(sequences_targets_rev["consensus"]["aligned_matches"].keys()):
			# generate a row
			sys.stdout.write('%9d' %id)
			sys.stdout.write('\t')
			sys.stdout.write('%25s' %sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"])
			sys.stdout.write('\t')

			if (pam_seq == "NGG") and (length == 20):

				for key in sequences_targets_rev.keys():
					seq_key = key

				if "seq_forscore" in sequences_targets_rev[seq_key]["aligned_matches"][targ_id]:
					sys.stdout.write('%35s' % sequences_targets_rev[seq_key]["aligned_matches"][targ_id]["seq_forscore"])
					sys.stdout.write('\t')
				else:
					sys.stdout.write('%35s' % "")
					sys.stdout.write('\t')


			sys.stdout.write("{0:2.2f}".format(GCpercent(sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"])))
			sys.stdout.write('\t')
			sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(sequences_targets_rev["consensus"]["aligned_matches"][id]["target_seq"])))

			if "exon_number" in sequences_targets_rev["consensus"]["aligned_matches"][id].keys():
				sys.stdout.write('\t')
				sys.stdout.write('%20s' %sequences_targets_rev["consensus"]["aligned_matches"][id]["exon_number"])

			if (pam_seq == "NGG") and (length == 20):
				# define a data structure to store scores for targets in individual sequences
				scores = []

				# iterate over all sequences in the alignment
				for seq_key in sequences_targets_rev.keys():
					if seq_key != "consensus":

						# calculate the score
						if "seq_forscore" in sequences_targets_rev[seq_key]["aligned_matches"][targ_id]:
							scores.append(calcDoenchScore(sequences_targets_rev[seq_key]["aligned_matches"][targ_id]["seq_forscore"]))

				sys.stdout.write('\t')

				if len(scores) > 0:
					sys.stdout.write("{0:3.2f}".format(sum(scores)/len(scores)))
				else:
					sys.stdout.write("")


			sys.stdout.write('\n')

		print("</textarea>")		
		print("</div>")
		print("<br />")
	
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

				# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
				# together with its sequence
 
				if (pam_seq == "NGG") and len(match.group(1)) == 20:
					# check that the target site is sufficiently far from the end of the sequence
					if end < len(sequences_matches[id]["sequence"]) - 6:
						sequences_matches[id]["direct_targets"][match_count]["seq_forscore"] = sequences_matches[id]["sequence"][start-4: end + 6]


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

				# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
				# together with its sequence
 
				if (pam_seq == "NGG") and len(match.group(1)) == 20:
					# check that the target site is sufficiently far from the end of the sequence
					if end < len(sequences_matches[id]["rc_sequence"]) - 6:
						sequences_matches[id]["reverse_targets"][match_count]["seq_forscore"] = sequences_matches[id]["rc_sequence"][start-4: end + 6]

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
	

	
def compute_output_unique_targets(sequences, sgRNA_regex):

	# get all of the matches for all of the sequences
	# the structure of this dictionary can be found in the function regex_simple_match
	
	seq_matches = regex_simple_match(sequences, sgRNA_regex, pam_seq, length, oriented)
	
	id_count = 7	# an id number for tables for unique sequences starting after the first 4 possible tables for common targets
	
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

					# assign the sequence for scoring to the new unique targets
					if (pam_seq == "NGG") and len(curr_seq) == 20:
						if "seq_forscore" in seq_matches[seq_id]["direct_targets"][targ_id]:
							unique_targets[seq_id][uniqID]["seq_forscore"] = seq_matches[seq_id]["direct_targets"][targ_id]["seq_forscore"]
				

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

					# assign the sequence for scoring to the new unique targets
					if (pam_seq == "NGG") and len(curr_seq) == 20:
						if "seq_forscore" in seq_matches[seq_id]["reverse_targets"][targ_id]:
							unique_targets[seq_id][uniqID]["seq_forscore"] = seq_matches[seq_id]["reverse_targets"][targ_id]["seq_forscore"]
				
########################################################################
# Exon filtering algorithm for unique sgRNA targets 
########################################################################
# 1. Take the dictionary of unique sgRNA targets built in the above code
# 2. Iterate over all sequence identifiers and all computed targets
# 
# 3. Retrieve coordinates for the simple match of a target in the sequence. 
#    Re-compute these coordinates in case when guide RNA matches in the reverse strand.
# 
# 4. Run these target coordinates against exon coordinates inside the transcript.
# 
#    The rules for testing if a target is completely inside an exon are the following:
#
#    - (start of target is > than the start of exon) AND (end of target is > than end of exon)
#    - (start of target is < than start of exon) AND (end of target is > than start of exon)
#    - An alternative is to use the test of string presense in an exon sequence in a particular exon, but this may be
#    computationally expensive

# # define the dictionary to store the exon-filtered unique targets
	# unique_targets_filtered = {}
	
# # iterate over all sequences and all computed targets
	# for	seq in unique_targets.keys():
		# # make a new unique identifier for targets which passed the exon filtering algorithm
		# new_uniqID = 1
		
		# if len(unique_targets[seq]) > 0:
		
			# # initialize a dictionary for a particular sequence
			# unique_targets_filtered[seq] = {}
			
			# for target in unique_targets[seq].keys():
				# # retrieve the start and end of the full match
				# # define the coordinates according to the orientation of the match
				# if unique_targets[seq][target]["orientation"] == "forward":
				
					# # obtain the sequence of the target
					# target_seq = seq_matches[seq]["sequence"][unique_targets[seq][target]["coordinates"][0]: unique_targets[seq][target]["coordinates"][1] + len(pam_seq)]

				# elif unique_targets[seq][target]["orientation"] == "reverse":
					
					# # obtain the sequence of the target
					# target_seq = rc(seq_matches[seq]["rc_sequence"][unique_targets[seq][target]["coordinates"][0]: unique_targets[seq][target]["coordinates"][1] + len(pam_seq)])
		
				# # define a boolean variable to check if the target is inside a particular exon
				# # it is defined as True initially because it is a priori more likely not to overlap two exons
				# # once one of the conditions is satisfied, it becomes False and subsequent code execution 
				# # occurs differently
				# inside = False

				# # check if this target does not overlap two exons			
				# for exonID in transcr_exons[seq]["exons"].keys():
					
					# if target_seq in transcr_exons[id]["exons"][exonID]["exon_seq"]:
						# inside = True
					
				# # proceed only if the match lies completely inside an exon of a particular transcript
				# if inside:
					# # initialize a dictionary for a particular target
					# unique_targets_filtered[seq][new_uniqID] = {}
					
					# # store all the relevant information in the new dictionary
					# unique_targets_filtered[seq][new_uniqID]["sequence"] = unique_targets[seq][target]["sequence"]
					# unique_targets_filtered[seq][new_uniqID]["orientation"] = unique_targets[seq][target]["orientation"]
					# unique_targets_filtered[seq][new_uniqID]["coordinates"] = unique_targets[seq][target]["coordinates"]
					
					# # increment the new_uniqID variable
					# new_uniqID += 1
				
				
	########################################################################
	# End of algorithm
	########################################################################
	
	
	####################################################################################
	# Output all of the unique targets in the form of tables, a table for each sequence 
	####################################################################################
	
	print("<div style='background-color:#914D87; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Unique sgRNA Targets</h3></div>")

	# output a sequence with highlighted targets
	unique_targets_sequences = {}
	
	for seq_id in unique_targets.keys():
		
		if len(unique_targets[seq_id]) == 0:
			print("<h4>The sequence %s does not have any unique sgRNA targets</h4>" %seq_id)
			print("<br />")
		else:
			# provide a link to the javascript program handling hide/expand event
			print """
			<h4><a href="javascript:Toggle%d();">Unique sgRNA targets for the sequence %s - expand or hide</a></h4>
			<div id="ToggleTarget%d">
			""" %(id_count, seq_id, id_count)
			id_count += 1
			
			# a dictionary to store information on the sgRNA targets for the current sequence
			current_sequence = {}
			
			# move over the unique_targets dictionary and other relevant data and 
			# produce the dictionary compatible with the original highlight_targets_output function
			
			unique_targets_sequences[seq_id] = {}
			unique_targets_sequences[seq_id]["sequence"] = seq_matches[seq_id]["sequence"]
			unique_targets_sequences[seq_id]["rc_sequence"] = seq_matches[seq_id]["rc_sequence"]
			unique_targets_sequences[seq_id]["direct_targets"] = {}
			unique_targets_sequences[seq_id]["reverse_targets"] = {}

			# initialize and update the current_sequence dictionary
			current_sequence[seq_id] = {}
			current_sequence[seq_id]["sequence"] = seq_matches[seq_id]["sequence"]
			current_sequence[seq_id]["rc_sequence"] = seq_matches[seq_id]["rc_sequence"]
			current_sequence[seq_id]["direct_targets"] = {}
			current_sequence[seq_id]["reverse_targets"] = {}
			 
			
			for targ_id in unique_targets[seq_id].keys():
				# check the orientation of each target and populate the corresponding dictionary of unique_targets_sequences
				
				if unique_targets[seq_id][targ_id]["orientation"] == 'forward':
					unique_targets_sequences[seq_id]["direct_targets"][targ_id] = {}
					unique_targets_sequences[seq_id]["direct_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					unique_targets_sequences[seq_id]["direct_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]

					current_sequence[seq_id]["direct_targets"][targ_id] = {}
					current_sequence[seq_id]["direct_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					current_sequence[seq_id]["direct_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]					

					
	
				elif unique_targets[seq_id][targ_id]["orientation"] == "reverse":
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id] = {}
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					unique_targets_sequences[seq_id]["reverse_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]


					current_sequence[seq_id]["reverse_targets"][targ_id] = {}
					current_sequence[seq_id]["reverse_targets"][targ_id]["target_seq"] = unique_targets[seq_id][targ_id]["sequence"]
					current_sequence[seq_id]["reverse_targets"][targ_id]["coord_target"] = unique_targets[seq_id][targ_id]["coordinates"]

					
			
			# print the visual view of the unique targets
			highlight_targets_single(current_sequence)
			
			
			# first print the colored heading
			print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Table view of unique targets</h3></div>")
			
			# obtain a list of data for unique targets both in the forward and reverse strands
			unique_forward = []
			unique_reverse = []

			for targ_id in sorted(unique_targets[seq_id].keys()):
				# targets in the forward strand
				if unique_targets[seq_id][targ_id]["orientation"] == "forward":			
					unique_forward.append(unique_targets[seq_id][targ_id]["sequence"])
				# targets in the reverse strand
				if unique_targets[seq_id][targ_id]["orientation"] == "reverse":			
					unique_reverse.append(unique_targets[seq_id][targ_id]["sequence"])

			# test if there are any targets in the forward strand
			if unique_forward:
				print """
				<h4><a href="javascript:Toggle%d();">Table of the unique targets in the forward strand - expand or hide</a></h4>
				<div id="ToggleTarget%d">
				""" %(id_count, id_count)
		
				id_count += 1

			
				# print the <div> element defining the pretty table style
				print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
				# header for the forward strand consensus sequence target data
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

				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "forward":
						if "exon_number" in unique_targets[seq_id][targ_id].keys():	
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

				# target data output for this table goes here
				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "forward":
				
						print("<tr>")
						print('<td>%d</td>' %targ_id)
						print('<td>%s</td>' %unique_targets[seq_id][targ_id]["sequence"])

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								print('<td>%s</td>' %unique_targets[seq_id][targ_id]["seq_forscore"])
							else:
								print('<td>%s</td>' % "")

						print('<td>%d</td>' %(unique_targets[seq_id][targ_id]["coordinates"][0] + 1))
						print('<td>%d</td>' %(unique_targets[seq_id][targ_id]["coordinates"][1] + 1))

						print('<td>')
						print "{0:2.2f}".format(GCpercent(unique_targets[seq_id][targ_id]["sequence"]))
						print('</td>')

						print('<td>')
						print "{0:2.2f}".format(DNA_RNA_Tm(unique_targets[seq_id][targ_id]["sequence"]))
						print('</td>')

						if "exon_number" in unique_targets[seq_id][targ_id].keys():
							print('<td>%s</td>' %unique_targets[seq_id][targ_id]["exon_number"])	

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								print('<td>')
								print "{0:3.2f}".format(calcDoenchScore(unique_targets[seq_id][targ_id]["seq_forscore"]))
								print('</td>')
							else:
								print('<td>%s</td>' % "")	
					
						print("</tr>")
			
				print """
				</table>
				</div>
				<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>		
				"""

				#output the textarea with the same information as above	
				print """
				<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
				(<i>e.g.</i> Excel) or text editors.</p>
				<textarea rows=10 cols=175 class=mono>
				"""

				sys.stdout.write('%9s' % "Target ID")
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

				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "forward":
						if "exon_number" in unique_targets[seq_id][targ_id].keys():	
							exon_info_available = True
				
				if exon_info_available:
					sys.stdout.write('\t')
					sys.stdout.write('%20s' % "Exon")

				if (pam_seq == "NGG") and (length == 20):
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "Score")
		
				sys.stdout.write('\n')


				# target data output for this table goes here
				for targ_id in sorted(unique_targets[seq_id].keys()):

					if unique_targets[seq_id][targ_id]["orientation"] == "forward":
						sys.stdout.write('%9d' %targ_id)
						sys.stdout.write('\t')
						sys.stdout.write('%25s' %unique_targets[seq_id][targ_id]["sequence"])
						sys.stdout.write('\t')

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								sys.stdout.write('%35s' % unique_targets[seq_id][targ_id]["seq_forscore"])
								sys.stdout.write('\t')
							else:
								sys.stdout.write('%35s' % "")
								sys.stdout.write('\t')

						sys.stdout.write('%5d' %(unique_targets[seq_id][targ_id]["coordinates"][0] + 1))
						sys.stdout.write('\t')
						sys.stdout.write('%4d' %(unique_targets[seq_id][targ_id]["coordinates"][1] + 1))
						sys.stdout.write('\t')
						sys.stdout.write("{0:2.2f}".format(GCpercent(unique_targets[seq_id][targ_id]["sequence"])))
						sys.stdout.write('\t')
						sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(unique_targets[seq_id][targ_id]["sequence"])))

						if "exon_number" in unique_targets[seq_id][targ_id].keys():
							sys.stdout.write('\t')
							sys.stdout.write('%20s' %unique_targets[seq_id][targ_id]["exon_number"])

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								sys.stdout.write('\t')
								sys.stdout.write("{0:3.2f}".format(calcDoenchScore(unique_targets[seq_id][targ_id]["seq_forscore"])))
							else:
								sys.stdout.write('\t')
								sys.stdout.write('%5s' % "")

						sys.stdout.write('\n')

				print("</textarea>")	
				print("</div>")


			if unique_reverse:
				print """
				<h4><a href="javascript:Toggle%d();">Table of the unique targets in the reverse strand - expand or hide</a></h4>
				<div id="ToggleTarget%d">
				""" %(id_count, id_count)
		
				id_count += 1

			
				# print the <div> element defining the pretty table style
				print('<div class="CSSTableGenerator" style="width:1050px; margin:10px;">')
		
				# header for the forward strand consensus sequence target data
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

				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "reverse":
						if "exon_number" in unique_targets[seq_id][targ_id].keys():	
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


		
				# target data output for this table goes here
				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "reverse":
				
						print("<tr>")
						print('<td>%d</td>' %targ_id)
						print('<td>%s</td>' %unique_targets[seq_id][targ_id]["sequence"])

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								print('<td>%s</td>' %unique_targets[seq_id][targ_id]["seq_forscore"])
							else:
								print('<td>%s</td>' % "")

						print('<td>%d</td>' %(unique_targets[seq_id][targ_id]["coordinates"][0] + 1))
						print('<td>%d</td>' %(unique_targets[seq_id][targ_id]["coordinates"][1] + 1))

						print('<td>')
						print "{0:2.2f}".format(GCpercent(unique_targets[seq_id][targ_id]["sequence"]))
						print('</td>')

						print('<td>')
						print "{0:2.2f}".format(DNA_RNA_Tm(unique_targets[seq_id][targ_id]["sequence"]))
						print('</td>')

						if "exon_number" in unique_targets[seq_id][targ_id].keys():
							print('<td>%s</td>' %unique_targets[seq_id][targ_id]["exon_number"])	

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								print('<td>')
								print "{0:3.2f}".format(calcDoenchScore(unique_targets[seq_id][targ_id]["seq_forscore"]))
								print('</td>')
							else:
								print('<td>%s</td>' % "")
							
						print("</tr>")
			
				print """
				</table>
				</div>
				<p class='normalp'>You can use the target sequences to generate your sgRNA expression constructs.</p>		
				"""

				#output the textarea with the same information as above	
				print """
				<p class ="normalp">Tab-delimited text can be copy-pasted into spreadsheet softwares
				(<i>e.g.</i> Excel) or text editors.</p>
				<textarea rows=10 cols=175 class=mono>
				"""


				sys.stdout.write('%9s' % "Target ID")
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

				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "reverse":
						if "exon_number" in unique_targets[seq_id][targ_id].keys():	
							exon_info_available = True
				
				if exon_info_available:
					sys.stdout.write('\t')
					sys.stdout.write('%20s' % "Exon")

				if (pam_seq == "NGG") and (length == 20):
					sys.stdout.write('\t')
					sys.stdout.write('%5s' % "Score")
		
				sys.stdout.write('\n')



				# target data output for this table goes here
				for targ_id in sorted(unique_targets[seq_id].keys()):
					if unique_targets[seq_id][targ_id]["orientation"] == "reverse":

						sys.stdout.write('%9d' %targ_id)
						sys.stdout.write('\t')
						sys.stdout.write('%25s' %unique_targets[seq_id][targ_id]["sequence"])
						sys.stdout.write('\t')

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								sys.stdout.write('%35s' % unique_targets[seq_id][targ_id]["seq_forscore"])
								sys.stdout.write('\t')
							else:
								sys.stdout.write('%35s' % "")
								sys.stdout.write('\t')

						sys.stdout.write('%5d' %(unique_targets[seq_id][targ_id]["coordinates"][0] + 1))
						sys.stdout.write('\t')
						sys.stdout.write('%4d' %(unique_targets[seq_id][targ_id]["coordinates"][1] + 1))
						sys.stdout.write('\t')
						sys.stdout.write("{0:2.2f}".format(GCpercent(unique_targets[seq_id][targ_id]["sequence"])))
						sys.stdout.write('\t')
						sys.stdout.write("{0:2.2f}".format(DNA_RNA_Tm(unique_targets[seq_id][targ_id]["sequence"])))

						if "exon_number" in unique_targets[seq_id][targ_id].keys():
							sys.stdout.write('\t')
							sys.stdout.write('%20s' %unique_targets[seq_id][targ_id]["exon_number"])

						if (pam_seq == "NGG") and (length == 20):
							if "seq_forscore" in unique_targets[seq_id][targ_id]:
								sys.stdout.write('\t')
								sys.stdout.write("{0:3.2f}".format(calcDoenchScore(unique_targets[seq_id][targ_id]["seq_forscore"])))
							else:
								sys.stdout.write('\t')
								sys.stdout.write('%5s' % "")

						sys.stdout.write('\n')

				print("</textarea>")	
				print("</div>")

			print """
			</div>
			<br />
			"""

def highlight_targets_single(sequences_targets): 
 	
	for seq in sequences_targets.keys():
		# print basic headers 
		
		print("<div style='background-color:#737CA1; width: 1050px; height: 50px; color: white; margin:10px;'><h3>Visual View of unique sgRNA Targets</h3></div>")

		################################
		# Forward strand
		################################
		print("<p style='font-family: Times New Roman, Times, serif; font-size: 18px;'><strong>Guide RNA targets have been highlighted. Long stretches of highlighted sequence are due to target overlap.<br />Please see tables for individual sequences.</strong></p>")
		sys.stdout.write("<div style='background-color:#0000DD; color: white; width: 500px; height: 25px; margin:10px;'><p class='normalp'><strong>Targets in the forward strand:</strong></p></div>")
		
		# print the <pre> element which will contain all of the sequence
		sys.stdout.write("<pre style='margin:10px; padding:5px; line-height: 1.6em;'>")
		
		# store the sequence from the input dictionary in a variable
		# to make it easier to work with these data
		sequence = sequences_targets[seq]["sequence"]
		sequence = sequence.rstrip()
		direct_targets = {}
		direct_targets = sequences_targets[seq]["direct_targets"]
	
		# do a test for the possibility that no targets of sgRNA have been identified
		if len(direct_targets) == 0:
			sys.stdout.write("This sequence does not contain any sgRNA targets in the sense orientation.\n\n")
	
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

		######################################################
		# End of start coordinates list of coordinates
		######################################################
		
		# initialize the positions to use while producing output
		curr_pos = 0
		curr_slice = 100

		# print markers
		print_markers_uniq(start_coords, curr_slice)
		
		# print the initial coordinate
		sys.stdout.write('{0:5d}'.format(curr_slice - 99))
		sys.stdout.write("  ")
		
		# looping over sgRNA targets:
		for targ_id in sorted(direct_targets.keys()):
						
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

					# print markers
					print_markers_uniq(start_coords, curr_slice)
					
					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")					
				
				elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = direct_targets[targ_id]["coord_target"][1]

					# print markers
					print_markers_uniq(start_coords, curr_slice)
					
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

				# print markers
				print_markers_uniq(start_coords, curr_slice)
	
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")	
				
				# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
				while first >= curr_slice:
					
					sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100

					# print markers
					print_markers_uniq(start_coords, curr_slice)
					
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

					# print markers
					print_markers_uniq(start_coords, curr_slice)

					sys.stdout.write('{0:5d}'.format(curr_slice - 99))
					sys.stdout.write("  ")
					
				elif direct_targets[targ_id]["coord_target"][1] > curr_slice:
					last = curr_slice
					sys.stdout.write("<strong><font style='BACKGROUND-COLOR: #0000CC; color: white'>%s</font></strong>" %sequence[curr_pos: last])     # mark with highlighted font
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100
					last = direct_targets[targ_id]["coord_target"][1]

					# print markers
					print_markers_uniq(start_coords, curr_slice)

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

			# print markers
			print_markers_uniq(start_coords, curr_slice)
			
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
				print_markers_uniq(start_coords, curr_slice)

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
		#CF5300
		#FFA500
		
		# print the <pre> element which will contain all of the sequence
		print("<pre style='margin:10px; padding:5px; line-height: 1.6em;'>")
		
		# store the sequence from the input dictionary in a variable
		# to make it easier to work with these data
		sequence = sequences_targets[seq]["rc_sequence"]
		sequence = sequence.lstrip()
		reverse_targets = {}
		reverse_targets = sequences_targets[seq]["reverse_targets"]
	
		# do a test for the possibility that no targets of sgRNA have been identified
		if len(reverse_targets) == 0:
			print("This sequence does not contain any sgRNA targets in the anti-sense orientation.\n")
	
		# iterate over all the targets of CRISPR sgRNAs in the current sequence
		# 
		# Use the following styles for the sgRNA targets and adjacent PAM sequences:
		# <pre>
		# <strong><font style="BACKGROUND-COLOR: #0000CC; color: white">sgRNA target</font></strong>
		# </pre>

		######################################################
		#  Start site marker labeling
		######################################################
		

		# generate a list of target start sites
		start_coords = []
	
		# looping over sgRNA targets:
		for targ_id in sorted(reverse_targets, key = lambda key: reverse_targets[key]["coord_target"][0]):
			begin = reverse_targets[targ_id]["coord_target"][0]
			start_coords.append(begin)

		######################################################
		# End of start coordinates list of coordinates
		######################################################
		
		# initialize the positions to use while producing output
		curr_pos = 0
		curr_slice = 100

		# print markers
		print_markers_uniq(start_coords, curr_slice)
		
		# print the initial coordinate
		sys.stdout.write('{0:5d}'.format(curr_slice - 99))
		sys.stdout.write("  ")	
		
		# looping over sgRNA targets:
		
		for targ_id in sorted(reverse_targets.keys()):
						
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
					print_markers_uniq(start_coords, curr_slice)

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
					print_markers_uniq(start_coords, curr_slice)

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
				print_markers_uniq(start_coords, curr_slice)
				
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")
				
				# a while loop to iterate in case of multiple lines   separating the sgRNA target and the current position
				while first >= curr_slice:
					
					sys.stdout.write("%s" %sequence[curr_pos: curr_slice])
					sys.stdout.write("\n")
					
					curr_pos = curr_slice
					curr_slice += 100

					# print markers
					print_markers_uniq(start_coords, curr_slice)
					
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
					print_markers_uniq(start_coords, curr_slice)
					
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
					print_markers_uniq(start_coords, curr_slice)
					
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
			print_markers_uniq(start_coords, curr_slice)

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
				print_markers_uniq(start_coords, curr_slice)
				
				sys.stdout.write('{0:5d}'.format(curr_slice - 99))
				sys.stdout.write("  ")
			
			if curr_slice >= len(sequence):
				# output the rest of the sequence
				sys.stdout.write("%s" %sequence[curr_pos: ])
				sys.stdout.write("\n")
		print("</pre>")

	
################################
# form data variable assignment 
################################

# obtain the information from the form
form = cgi.FieldStorage(keep_blank_values=1)

# form fields related to the properties of the CRISPR guide RNA

# first two nucleotides - of some importance in current sgRNA molecules
dinuc = form["5dinuc"].value

# overall length
length = form["target_length"].value

# error checking of the target length input
try:
	length = int(length)
except ValueError:
	badrequest_length(length)

# orientation of the PAM sequence
oriented = form["oriented"].value
oriented = int(oriented)

# There is no point in allowing mismatches in this script, but to keep functions the same
# a mismatch variable has been created and assigned "No" value
mismatch = "No"

# an option where an user can choose whether to use the fixed most common PAM sequence
pam_value = form["PAM"].value
pam_value = int(pam_value)

if pam_value == 1:
	pam_seq = "NGG"
elif pam_value == 2:
	pam_seq = form["PAM_seq"].value

# obtain the species value from the form
species = form["species"].value

# retrieve the molecule type
mol_type = form["moltype"].value

# obtain the input gene identifier from the user 
gene_id = form["gene_id"].value
entered_geneid = gene_id


#####################################
# error checking and data validation
#####################################


# check and produce an error if the user did not specify their own PAM sequence
if not validate_pam(pam_seq):
	badrequest_pam(pam_seq)

if len(pam_seq) == 0:
	if pam_value == 2:
		badrequest_pam(pam_seq)
		
# Still, we have to check if the species value is empty or the gene with this name does not exist in this species 
# Other options are valid in principle but may not work in some circumstances with an input identifier
# but this will be tested separately

if species == "default":
	badrequest_species(species)

# first validate the gene identifier entered by the user
# the basic validation should be checking if the field was not empty
# the next step will be to check that it is possible to retrieve data from a database using this identifier, 
# which can also be part of the error checking process
			
# case 1: the gene identifier field was submitted empty
if gene_id == "":
# generate an error	
	badrequest_seq(gene_id)


# case 2: the gene identifier is not valid	
	
########################################
# ADD MORE VALIDATION FOR AN IDENTIFIER	
########################################
# - pattern matching of an identifier to find out what kind of identifier was entered
# - database query validation of an identifier and the species name

# clean up the gene_id variable by removing all possible spaces from it
gene_id = ''.join(gene_id.split())

# make the identifier uppercase, but then use case-insensitive comparison
gene_id.upper()

# also, establish a database connection
os.chdir("/home/multicri/cgi-bin")

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
else:
	badrequest_species(geneids_input)

# first check if the gene identifier has a format of an Ensembl ID
if mol_type == "ensgenes":

	# try to select the data from Genes table on the assumption that the gene identifier is of Ensembl ID format
	cur.execute("SELECT * from Genes WHERE geneid = ? AND species = ? COLLATE NOCASE", (gene_id, species))
	row = cur.fetchone()
	
	if row is None:
		# generate an error	
		badrequest_seq(gene_id)
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
		# generate an error	
		badrequest_seq(gene_id)
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

# a dictionary containing information about both the transcripts for the input gene
# and their constitutive exons including their coordinates
transcr_exons = {}

# another data structure will only contain the sequence information for transcripts
# under the transcriptID keys
transcripts = {}


# 2. Retrieve all the possible data items from the database and verify that this has worked well.
#    Consider raising some error if the process does not go as planned.

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
	
# block subsequent execution of the program if only one sequence is available for alignment
if len(transcripts) < 2:
	badrequest_input(gene_id)

	
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

	if end is not None:
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

# 3. Update subsequent algorithms.
# Appropriate code will be inserted in subsequent sections
# the most important test will be to check whether any given target is located within an exon
		
# The basic algorithm:

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
		
		
		
######################################################################
# Below are the instructions from the original crispr_targeter script		
######################################################################		
		
###########################################################################################################
# 1. Write the sequences of transcripts retrieved from the database to a file to prepare them for alignment
###########################################################################################################

# Reading has been done from the html form in the previous part of the script

# Preparing the data for the alignment will involve writing them to a temporary text file

# use tempfile functions to achieve this:
# http://pymotw.com/2/tempfile/
# http://docs.python.org/2/library/tempfile.html

# create a temporary file with a random name and write all the sequences to that file

# create a dummy file which will be used as an input for the ClustalW2 program
tmp_path = "/home/multicri/tmp/"

###########################################################
# Random folder generation for storing each user's results
##########################################################

rand_folder = str(''.join(random.sample(string.letters*8,8)))
new_folder = tmp_path + rand_folder + '/'
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
sequences = transcripts


for key in sequences.keys():
	
	id_line = ">" + key + " transcript" 
	
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
clustalw_exe = "/home/multicri/cgi-bin/clustalw2"

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
		
		# define a boolean variable, which is True if the target does not span two exons in any of the sequences
		cons_inside = True

		for record in align:
		#initialization for all of the aligned sequences
			sequences_matches_forward[record.id]["aligned_matches"][match_count] = {}
			sequences_matches_forward[record.id]["orig_matches"][match_count] = {}
		
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
			exon_N = ""
			
			# obtain the sequence of the target and check that it is inside a particular exon
			target_seq = seq[start: end]
			
			# check if this target does not overlap two exons			

			for exonID in sorted(gene_exons, key = lambda key: gene_exons[key]["start"]):
				if target_seq in gene_exons[exonID]["exon_seq"]:
					inside = True
					
					# determine exon number where this match occurred
					exon_N = exonID
					break

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
			
			
			# assign an exon number to the current common target site
			if exon_N:
				sequences_matches_forward["consensus"]["aligned_matches"][match_count]["exon_number"] = exon_N
			exon_N = ""
		
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

					# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
					# together with its sequence
 
					if (pam_seq == "NGG") and len(seq[start: end]) == 20:
						# check that the target site is sufficiently far from the end of the sequence
						if end < len(seq) - 6:
							 sequences_matches_forward[record.id]["aligned_matches"][match_count]["seq_forscore"] = seq[start-4: end + 6]

				
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

		for record in align:
		#initialization for all of the aligned sequences
			sequences_matches_reverse[record.id]["aligned_matches"][match_count] = {}
			sequences_matches_reverse[record.id]["orig_matches"][match_count] = {}

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

			# define a boolean variable to check if the target is inside a particular exon
			# it is defined as True initially because it is a priori more likely not to overlap two exons
			# once one of the conditions is satisfied, it becomes False and subsequent code execution 
			# occurs differently

			inside = False
			exon_N = ""
			
			# obtain the sequence of the target and check that it is inside a particular exon
			target_seq = rc(seq[start: end])
			
			# check if this target does not overlap two exons			
			for exonID in sorted(gene_exons, key = lambda key: gene_exons[key]["start"]):
				
				if target_seq in gene_exons[exonID]["exon_seq"]:
					inside = True
					# determine exon number where this match occurred
					exon_N = exonID
					break
			
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

			# assign an exon number to the current common target site
			if exon_N:
				sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["exon_number"] = exon_N
			exon_N = ""

		
			if oriented == 5:
				sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["coord_match"] = [start-len(pam_seq), end]
			elif oriented == 3:
				sequences_matches_reverse["consensus"]["aligned_matches"][match_count]["coord_match"] = [start, end + len(pam_seq)]
				
			for record in align:
				# matches in the original sequences:
				orig_seq = rc(str(sequences[record.id]))
				seq = rc(str(record.seq))

				gaps = int(seq[:match.start()].count('-'))
					
				# matches in the original sequences:

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

					# Nov 2014: scoring system for type II sgRNAs is now available so we can store the immediate neighbourhood of the target site
					# together with its sequence
 
					if (pam_seq == "NGG") and len(seq[start: end]) == 20:
						# check that the target site is sufficiently far from the end of the sequence
						if end < len(seq) - 6:
							sequences_matches_reverse[record.id]["aligned_matches"][match_count]["seq_forscore"] = seq[start-4: end + 6]
			
			# increment the count of matches
			match_count += 1			
	
#############################################################
#  CLEAN THE FILESYSTEM
#############################################################

os.remove(input_file)
os.remove(output_file)
os.remove(dnd_file)
os.rmdir(new_folder)
			
########################################################################################
# Output Phase: a highlighted alignment and tables with common and unique sgRNA targets.
########################################################################################

# The output of this program will consist of some headers identifying the sequence and
# what the results are, an alignment with highlighted common targets according to
# the alignment generated above and the searches for the targets that have been performed.
# In addition, researchers may be strongly interested in targets unique targets for each
# sequence. This need will be addressed by a different algorithm and since the results
# are not based on the multiple sequence alignments, theys will be presented in separate tables.
# The overall sequence of presented elements will be the following: headers, alignment with
# highlighted common targets, a table detailing the common targets and tables with the info
# on unique targets.

print "Content-type: text/html\n\n"

crispr_targeter_header(sequences)

###################################
# Output Functions
################################### 

# Output CRISPR sgRNA targets in the multiple sequence alignment

highlight_targets_output(sequences_matches_forward, sequences_matches_reverse)

# compute unique sgRNA targets in each sequence

compute_output_unique_targets(sequences, sgRNA_regex)


# Write some code to do verification and error checking of the above code especially for the reverse complement sequences

crispr_targeter_footer()