### Description of the CRISPR MultiTargeter repository

####  Software and availability

CRISPR MultiTargeter is a web-based tool for automatic searches of CRISPR guide RNA targets. It can find highly similar or identical target sites in multiple genes or transcripts or design targets unique to particular genes or transcripts. The search for common targets is based on generating a multiple sequence alignment, and unique targets are found using a string comparison algorithm among all possible targets for each sequence. The basic algorithm implemented in CRISPR MultiTargeter is versatile and can accommodate almost any possible target specificity of CRISPR/Cas system, which is important because new target specificities are discovered and will be used in genetic studies in different model systems.

CRISPR MultiTargeter is available for free use online at [http://multicrispr.net/](http://multicrispr.net/). Standalone versions for working with multiple genes or transcripts of a single gene have been developed and are available at [Figure4_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure4_data) and [Figure5_data](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/Figure5_data), respectively. These programs require functional SQLite3 databases, which could not be uploaded but can be constructed as described in the [database_construction](https://github.com/SergeyPry/CRISPR_MultiTargeter/tree/master/database_construction) folder of this repository.

