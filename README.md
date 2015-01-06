### Description of the CRISPR MultiTargeter repository

####  Software and availability

CRISPR MultiTargeter is a web-based tool for automatic searches of CRISPR guide RNA targets. It can find highly similar or identical target sites in multiple genes or transcripts or design targets unique to particular genes or transcripts. The search for common targets is based on generating a multiple sequence alignment, and unique targets are found using a string comparison algorithm among all possible targets for each sequence. The basic algorithm implemented in CRISPR MultiTargeter is versatile and can accommodate almost any possible target specificity of CRISPR/Cas system, which is important because new target specificities are discovered and will be used in genetic studies in different model systems.

CRISPR MultiTargeter is available for free use online at [http://multicrispr.net/](http://multicrispr.net/). Standalone versions for working with multiple genes or transcripts of a single gene have been developed and are available at [Figure4_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure4_data) and [Figure5_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure5_data), respectively. These programs require functional SQLite3 databases, which could not be uploaded but can be constructed as described in the [database_construction](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/database_construction) folder of this repository.

#### Repository structure

This repository is a repository for both the CRISPR MultiTargeter website and the corresponding manuscript currently in revision. It consists of the following parts. 

* [Figure4_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure4_data) - Data and scripts to generate the Figure 4 of the paper. The detailed description of this folder is available in its own README.md file.
* [Figure5_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure5_data) - Data and scripts to generate the Figure 5 of the paper. The detailed description of this folder is available in its own README.md file.
* [cgi-bin](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/cgi-bin) - This folder contains all the scripts used on the [CRISPR MultiTargeter website](http://multicrispr.net/) as well as the Linux/Unix version of the ClustalW2 software required for finding common sgRNA target sites in multiple sequences.
* [database_construction](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/database_construction) - This folder contains the database schema for a single-species database as well as a detailed description on how such a database can be generated using publicly-available genomic data. Please note that the current implementation of the CRISPR MultiTargeter website requires that the databases be located in the same folder as the Python scripts.  
