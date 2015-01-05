###Workflow for generating Figure 5 of "CRISPR MultiTargeter: a web tool to find common and unique CRISPR single guide RNA targets in a set of similar sequences" manuscript

#### Step 1. Generating all possible ohnolog gene pairs in zebrafish.

#### Step 2. Finding common sites in pairs of ohnologs and counting them.

For this step you need to run **multigene_search_cmdline.py** on the **zebrafish.db** SQLite3 database using the **ohnologs_pairs_output.txt** file as input. The database was too big to include into this repository but its construction is described in the [Database construction section](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/database_construction). This script is designed to run on Windows using the Windows version of the ClustalW2 executable program. To run it, you will need to install ClustalW2 and configure the **clustalw_exe** variable, for example:
clustalw_exe = r"C:\\Program Files\\ClustalW2\\clustalw2.exe".

The output from this step are two files: **ohnologs_total_sites.txt** and **ohnolog_target_sites_details.txt**. The **ohnologs_total_sites.txt** contains total counts of common CRISPR/Cas9 sgRNA sites for all ohnolog pairs analyzed in this workflow. The **ohnolog_target_sites_details.txt** file contains the Ensemble identifiers of both ohnologs in a pair as well as counts of common sgRNA sites in both sense and anti-sense strands of these genes, for example
**ENSDARG00000071685	ENSDARG00000086104	4	4**.

#### Step 3. Generating data for common sgRNA target sites overall statistics and making the Figure 5A plot.

The proportions of "Ohnolog pairs with common target sites", "Sites with mismatches", "Sites without mismatches", "Sites in the sense strand", "Sites in the anti-sense strand" were calculated using **ohnolog_summarizer.py** and the resulting values were pasted into the **Figure5A.R** R script.

**Figure5A.R** script was run to produce the plot in Figure 5A. Please see the Figure 5 included in this repository.


#### Step 4. Plotting the distribution of common sgRNA counts in ohnolog pairs for the Figure 5B.