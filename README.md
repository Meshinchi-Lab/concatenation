# General R Utilities 
### for concatenation


**Important notes on the RNA-seq fusions data:**
1) Each of the RNA-seq fusion datasets (CICERO, STAR-Fusion, and TransAbyss) are initially cleaned individually. 

For example, the scripts `TARGET_AML_1031_Cat_StarFusion.Rmd`, `TARGET_AML_1031_Cat_TransAbyss_Fusions_5.7.18.Rmd` were each run independely as soon as the raw patient level data was avaialble. The CICERO data was provided in a single cleaned file from St. Jude collaborators, and this dataset was cleaned seperately in `2020.04.13_CICERO_St.Jude/TARGET_AML_RBD_Clean_CICERO_Fusions.Rmd`.

The cleaned output files from these data concatentation/cleaning scripts follow the approximate naming convention `TARGET_AML_{COHORT}_{FUSIONALGORITHM}_reformatted_FilteredForNBM_PrimaryFusions_{DATE}.csv`. 

2) The cleaned and concatenated fusion files are then merged (full join) per patient sample in the script `TARGET_AML_1031_0531_Relapse_Combine_STAR_TransAbyss_CICERO_Fusion_Calls.Rmd`. 

The concensus between the three fusion algorithms is noted by presence/absence and the specific breakpoint junction coordinates. The level of concensus between the fusion algorithms, and the karyotype data where available, is then summarized as a "confidence level" column. The confidence level 1 are the highest confidence of being a true positive, and confidence level 4 are one-offs that are considered the lowest confidence of being a true positive. 




Author: Jenny Leopoldina Smith<br>
ORCID: [0000-0003-0402-2779](https://orcid.org/0000-0003-0402-2779)
<br>
