### CRISPR MultiTargeter SQLite3 database construction

#### Genomic data download and preparation

Ensembl BioMart data were downloaded for exons, transcripts and basic gene information.
The data items that need to be downloaded are listed in the following database definition. Processing of the resulting files from BioMart into a database can be done directly or through an intermediate step of generating text files containing tables. I used the second option to generate the database files.

Gene sequences were assembled from exon sequences according to their coordinates. This was accomplished by running the following script with files such as  "example_exons_new_BioMart.txt" and "genes-symbols_BioMart.txt":

```python
python gene_allexons_merger.py
```

Please see the input files and the script for details.

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

#### Inserting the Exon data

For inserting exon information, the user needs to run two scripts:

The following script requires a file such as "example_zebrafish_exon_table.txt".
It will insert main data on exons: ExonID, EnsemblGeneID, sequence, strand and chromosome coordinates.
```python
python insert_exons.py
```

The next script requires a file such as "example_zebrafish_gene-exon_table.txt". It will insert the coordinates of an exon inside a merged gene model.
```python
python insert_exon_coords.py
```

#### Inserting Transcript data

To insert the transcript information containing EnsemblTranscriptID, EnsemblGeneID and the transcript sequence, it is necessary to simply run the following script with a file such as "example_zebrafish_transcript_table.txt":
```python
python insert_transcripts.py
```

#### Inserting Genes data

To insert the gene information containing EnsemblGeneID, gene symbol, species and the merged gene sequence, it is necessary to simply run the following script with a file such as "example_zebrafish_genes_table.txt":
```python
python insert_genes.py
```

#### Final things

To make transactions faster, we also created a number of indices for most frequent transactions. The following script contains the required code.
```python
python indexer.py
```

**Disclaimer:** The code and example data in this repository folder is only provided for illustrating how we generated our database. Anyone interested in using this code for their databases will need complete datasets from their species of choice and will need to make many adaptions of the code presented here.




