import os
import sqlite3 as db

# the purpose of this script is to add a couple of columns to the Exons table
# and remove the Genes2Exons and TranscriptsExons tables

# connect to the database
conn = db.connect('zebrafish.db')
cur = conn.cursor()

cur.execute("""ALTER TABLE Exons ADD COLUMN genestart INT;""")
conn.commit()
cur.execute("""ALTER TABLE Exons ADD COLUMN geneend INT;""")
conn.commit()


cur.execute('select * from Exons')
r=cur.fetchone()
print cur.description
conn.close()

