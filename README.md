# General R Utilities 
### for concatenation


#### RNA-seq fusions data

1) Each of the RNA-seq fusion datasets (CICERO, STAR-Fusion, and TransAbyss) are initially cleaned individually. 

For example, the scripts `TARGET_AML_1031_Cat_StarFusion.Rmd`, `TARGET_AML_1031_Cat_TransAbyss_Fusions_5.7.18.Rmd` were each run independely as soon as the raw patient level data was avaialble. The CICERO data was provided in a single cleaned file from St. Jude collaborators, and this dataset was cleaned seperately in `2020.04.13_CICERO_St.Jude/TARGET_AML_RBD_Clean_CICERO_Fusions.Rmd`.

The cleaned output files from these data concatentation/cleaning scripts follow the approximate naming convention `TARGET_AML_{COHORT}_{FUSIONALGORITHM}_reformatted_FilteredForNBM_PrimaryFusions_{DATE}.csv`. 

2) The cleaned and concatenated fusion files are then merged (full join) per patient sample in the script `TARGET_AML_1031_0531_Relapse_Combine_STAR_TransAbyss_CICERO_Fusion_Calls.Rmd`. 

The concensus between the three fusion algorithms is noted by presence/absence and the specific breakpoint junction coordinates. The level of concensus between the fusion algorithms, and the karyotype data where available, is then summarized as a "confidence level" column. The confidence level 1 are the highest confidence of being a true positive, and confidence level 4 are one-offs that are considered the lowest confidence of being a true positive. 

#### RNA-seq Quantification Data 

1) Each batch of RNA-seq quantification data are first processed individually. For example, the down syndrome AML RNA-seq data was downloaded from S3 and then concatenated using `tximport` R package in `concat_Kallisto/TARGET_AML_RBD_Kallisto_Quant_DSAML_Concatenation.Rmd`. 

The output matrices are then saved in two locations A) fast drive under `/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3` in gene level or transcript level direcories, and B) `/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/0000.00.03_ExpressionMatrices`. The second directory is the main "working copy" and if anything damages it, the clean copies can be found in level 3 directories. 

2) After each new batch of RNA-seq is concatenated and save in the `0000.00.03_ExpressionMatrices` directory, run the `TARGET_AML_Combine_Ribodepleted_RNAseq_Matrices.Rmd` rmarkdown. 

This will combine the multiple large matrices (.csv/.txt or .RDS files), which will also run the `kallisto_rmDups()` function to remove duplicate gene names. That final de-duped matrix is the main shareable counts matrix - currently > 3000 samples. 

The Sample manifests were then updated for the concatenated final matrix using the `TARGET_AML_Sequencing_Manifests_05.05.21.Rmd` file under the `SequencingDataMatrix` github repository. This script is a data cleaning script, and can be used a template. But it probably should only be edited and modified for ~1 year. Then start a fresh data matrix generation script using it as a template. 



Author: Jenny Leopoldina Smith<br>
ORCID: [0000-0003-0402-2779](https://orcid.org/0000-0003-0402-2779)
<br>
