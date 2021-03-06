# The purpose of this script is to insert the data into the tables
# of genomes database
import os
import sqlite3 as db


# obtain all the files from a directory matching a certain naming scheme
# assign these filenames to the corresponding variables

for file in os.listdir("."):

################################
# Insertion of GENE data
################################
		
	if file.endswith("genes_table.txt"):
        # open the file containing exons
		tr_fh = open(file, 'r')
		# also, establish a database connection
		conn = db.connect('zebrafish.db')
		cur = conn.cursor()
		line_count = 0
		insert_count = 0
		
		# go iteratively through each line of the file
		for line in tr_fh:
			# check if you reached the end of the file
			if line != "":
				line_count = line_count + 1
				results = line.split('\t')
				if len(results) > 3:
				
					# read the data items for database insertion
					(geneid, symbol, species, sequence) = line.split('\t')
					sequence.rstrip()
					# execute the insertion statement
					cur.execute('''INSERT INTO Genes (geneid, symbol, sequence, species) VALUES(?,?,?,?)''', (geneid, symbol, sequence, species))
					insert_count = insert_count + cur.rowcount
				
		conn.commit()
		conn.close()
		tr_fh.close()
		
		print(file)
		print("The number of lines gone through is %d" %line_count)
		print("The number of rows inserted is %d" %insert_count)		
