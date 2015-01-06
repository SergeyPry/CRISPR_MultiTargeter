### CRISPR MultiTargeter SQLite3 database construction

#### Genomic data download and preparation

#### Database creation

The following code was run to create a single-species database

```python
import sqlite3
## the database can have any filename, e.g. "zebrafish.db"
conn = sqlite3.connect('zebrafish.db')
print "Database created successfully";

conn.execute("""CREATE TABLE Genes 
		(geneid CHAR(30) PRIMARY KEY NOT NULL, 
		symbol CHAR(30), 
		sequence TEXT NOT NULL,
		species CHAR(100) NOT NULL);""")
  		
conn.execute("""CREATE TABLE Exons
       (exonID CHAR(30) PRIMARY KEY	NOT NULL,
		geneid CHAR(30) NOT NULL,
		sequence TEXT NOT NULL,
		strand SMALLINT,
		chrstart INT NOT NULL,
        chrend INT NOT NULL,
		FOREIGN KEY(geneid) REFERENCES Genes(geneid));""") 

conn.execute("""CREATE TABLE Transcripts 
		(transcriptID CHAR(30) PRIMARY KEY NOT NULL,
		geneid CHAR(30) NOT NULL,
		sequence TEXT NOT NULL,
		FOREIGN KEY(geneid) REFERENCES Genes(geneid));""")
		
# try this out and experiment with inserting and querying data to make sure
# it works in the real script
		
conn.commit()
conn.close()
```

Because the database definition was evolving, the coordinates of exons inside a particular merged gene had to be inserted later. The database structure had to be updated as well.
```python
import os
import sqlite3 as db

# the purpose of this script is to add a couple of columns to the Exons table
# connect to the database
conn = db.connect('zebrafish.db')
cur = conn.cursor()

cur.execute("""ALTER TABLE Exons ADD COLUMN genestart INT;""")
conn.commit()
cur.execute("""ALTER TABLE Exons ADD COLUMN geneend INT;""")
conn.commit()
```

#### Inserting the data

For inserting exon information, the user needs to run two scripts:

```python
## the following script requires a file such as "example_zebrafish_exon_table.txt"
## it will insert main data on exons: ExonID, EnsemblGeneID, sequence, strand and chromosome coordinates
python insert_exons.py

## the following script requires a file such as "example_zebrafish_gene-exon_table.txt"
## it will insert only the coordinates of an exon inside a merged gene model
python insert_exon_coords.py
```



