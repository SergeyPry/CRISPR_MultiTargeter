import sqlite3

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

