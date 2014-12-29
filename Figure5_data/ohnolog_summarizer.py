# This script is necessary to generate data for the Figure4

# Figure4 - part A:

# 1. Proportion of ohnolog pairs with identifiable common target sites.
# 2. Proportion of sites with mismatches.
# 3. Proportion of sites without mismatches.
# 4. Proportion of sites in the forward strand.
# 5. Proportion of sites in the reverse strand.

# Note:
# 1,2,3 need to be calculated based on the data in the "details" file
# 4,5 need to be calculated based on the counts file


# Figure4 - part B:

#  histogram of total counts of common target sites. This should be based on the "counts" file 
#  of common target sites.

# open files necessary for calculations
counts = open("ohnolog_target_sites_counts.txt", 'r')
details = open("ohnolog_target_sites_details.txt", 'r')

# total number of pairs
Npairs = 6305
Npairs = float(Npairs)
############################################
# Calculate all the numbers
############################################

# define proportions variables for all pieces of data
Prop_ohnologs_with_sites = 0.000
Prop_sites_mismatches = 0.000
Prop_sites_exact = 0.000
Prop_sites_forward = 0.000
Prop_sites_reverse = 0.000

# 1. Proportion of ohnolog pairs with identifiable common target sites.

counts_list = counts.readlines()
N_with_targets = float(len(counts_list))
Prop_ohnologs_with_sites = N_with_targets/Npairs


# 2. Proportion of sites with mismatches.

N_with_mism = 0
details_list = details.readlines()


for line in details_list:
	if "X" in line:
		N_with_mism += 1

N_with_mism = float(N_with_mism)
Prop_sites_mismatches = N_with_mism/float(len(details_list))

# 3. Proportion of sites without mismatches
Prop_sites_exact = 1 - Prop_sites_mismatches

# 4. Proportion of sites in the forward strand.
N_sites_forward = 0

for line in details_list:
	if "forward" in line:
		N_sites_forward += 1
N_sites_forward = float(N_sites_forward)
Prop_sites_forward = N_sites_forward/float(len(details_list))
		
# 5. Proportion of sites in the reverse strand.
N_sites_reverse = 0

for line in details_list:
	if "reverse" in line:
		N_sites_reverse += 1
N_sites_reverse = float(N_sites_reverse)	
Prop_sites_reverse = N_sites_reverse/float(len(details_list))

# prepare data for the histogram
total_counts = open("ohnologs_total_sites.txt", "w")
total_counts.write("Total counts of target sites for ohnolog pairs")
total_counts.write("\n")

for line in counts_list:
	line = line.rstrip('\n')
	items = line.split('\t')
	forw_count = int(items[2])
	rev_count = int(items[3])
	total = forw_count + rev_count
	
	total_counts.write(str(total))
	total_counts.write("\n")

total_counts.close()


# write out the data for the proportions bar graph
proportions = open("proportions_data_ohnolog_analysis.txt", "w")

# print the header
proportions.write("Ohnolog pairs with common target sites")
proportions.write("\t")
proportions.write("Sites with mismatches")
proportions.write("\t")
proportions.write("Sites without mismatches")
proportions.write("\t")
proportions.write("Sites in the sense strand")
proportions.write("\t")
proportions.write("Sites in the anti-sense strand")
proportions.write("\t")
proportions.write("\n")
	
# print the data
proportions.write("%.3f" % Prop_ohnologs_with_sites)
proportions.write("\t")
proportions.write("%.3f" % Prop_sites_mismatches)
proportions.write("\t")
proportions.write("%.3f" % Prop_sites_exact)
proportions.write("\t")	
proportions.write("%.3f" % Prop_sites_forward)
proportions.write("\t")
proportions.write("%.3f" % Prop_sites_reverse)
proportions.write("\t")
proportions.write("\n")





