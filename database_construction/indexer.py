import sqlite3 as db

# establish a database connection
conn = db.connect('zebrafish.db')
cur = conn.cursor()

# index geneid and species columns in the Genes table
cur.execute("CREATE INDEX genes_species on Genes (geneid, species)") 

# index symbol and species columns in the Genes table
cur.execute("CREATE INDEX symbol_species on Genes (symbol, species)") 

# index geneid column in the Transcripts table
cur.execute("CREATE INDEX geneid_Transcripts on Transcripts (geneid)")

# index geneid column in the Exons table
cur.execute("CREATE INDEX geneid_Exons on Exons (geneid)")

# index geneid column in the Genes table
cur.execute("CREATE INDEX gene_ids on Genes (geneid)")

# index transcriptID column in the Transcripts table
cur.execute("CREATE INDEX trid_Transcripts on Transcripts (transcriptID)")

# index geneid column in the Genes table
cur.execute("CREATE INDEX refseq_ids on Refseq (refseqid)")


conn.commit()
conn.close()