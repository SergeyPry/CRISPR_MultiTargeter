###Workflow for generating Figure 4 of "CRISPR MultiTargeter: a web tool to find common and unique CRISPR single guide RNA targets in a set of similar sequences" manuscript

#### Step 1. Making a list of alternative transcripts.

This step makes a list of all alternative transcripts in zebrafish together with their corresponding genes based on the data available in Ensembl BioMart database. The **genes_selecter.py** script when run on **zebrafish_transcript_table.txt** accomplishes this task and produces **genes_with_alternative_transcripts.txt** file as output. However, Github has restrictions on the file sizes so the entire file was not possible to upload. We therefore opted to include **zebrafish_transcript_table_example.txt** file which has a similar structure. The complete file can be generated using BioMart.

#### Step 2. Finding the sites in alternative transcripts and counting them.

For this step you need to run **transcripts_cmdline.py** on the **zebrafish.db** SQLite3 database using the **genes_with_alternative_transcripts.txt** file as input. The database was too big to include into this repository but its construction is described in the [Database construction section](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/database_construction).

The output from this step are two files: **total_unique_sites.txt** and **detailed_unique_sites.txt**. The **total_unique_sites.txt** contains total counts of unique CRISPR/Cas9 sgRNA sites for all alternative transcripts analyzed in this workflow. The **detailed_unique_sites.txt** file contains the gene and transcript Ensemble identifiers of these alternative transcripts as well as counts of unique sgRNA sites in both sense and anti-sense strands of these transcripts, for example
**ENSDARG00000091715	ENSDART00000152509	14	6**.

#### Step 3. Generating data for unique sgRNA target site overall statistics and making the Figure 4A plot

The proportions of "Genes with isoform-specific sgRNA sites", "Transcripts with isoform-specific sgRNA sites", "Sites in the sense strand", "Sites in the anti-sense strand" were calculated based on the **detailed_unique_sites.txt** file in Excel and the resulting values were pasted into the **Figure4A.R** R script.

**Figure4A.R** script was run to produce the plot in Figure 4A. Please see the Figure 4 included in this repository.

#### Step 4. Plotting the distribution of unique sgRNA counts in alternative transcripts for the Figure 4B.

**Figure4B.R** script was run with the **total_unique_sites.txt** file to produce the plot in Figure 4B. Please see the Figure 4 included in this repository.






