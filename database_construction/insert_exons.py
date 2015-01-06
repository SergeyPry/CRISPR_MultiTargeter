# The purpose of this script is to insert the data into the tables
# of genomes database
import os
import sqlite3 as db


# obtain all the files from a directory matching a certain naming scheme
# assign these filenames to the corresponding variables

for file in os.listdir("."):

################################
# Insertion of EXON data
################################
	if file.endswith("_exon_table.txt"):
	# open the file containing exons
		ex_fh = open(file, 'r')
		# also, establish a database connection
		conn = db.connect('zebrafish.db')
		cur = conn.cursor()
		line_count = 0
		insert_count = 0
		
		# go iteratively through each line of the file
		for line in ex_fh:
			# check if you reached the end of the file
			if line != "":
				line_count = line_count + 1
				# clean up the line of text 
				line.strip('\n')
				# read the data items for database insertion
				(exonID, geneid, sequence, strand, chrstart, chrend) = line.split('\t')
				chrend.rstrip()
				# execute the insertion statement
				cur.execute('''INSERT INTO Exons (exonID, geneid, sequence, strand, chrstart, chrend) VALUES(?,?,?,?,?,?)''', (exonID, geneid, sequence, strand, chrstart, chrend))
				insert_count = insert_count + cur.rowcount
				
		conn.commit()
		conn.close()
		ex_fh.close()
		
		print(file)
		print("The number of lines gone through is %d" %line_count)
		print("The number of rows inserted is %d" %insert_count)