# predHBD
[![License: GPL-2.0](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)

This folder contains the scripts used in the study:
"Genomic prediction of individual inbreeding levels for the management of genetic diversity in populations with small effective size"

 General workflow :

 1) Generate the moderately inbred and highly inbred simulated populations
    Codes can be found in the /simulation directory

 2) Predict levels of homozygosity-by-descent (HBD) in a future offspring of a genotyped pair of parents using 16 methods
    Codes are in the /predict directory 
	   
	  * predict/input: to prepare the genomic data in the format required to run each prediction method
	  * predict/[prediction_method_name]: to run each prediction method

 3) Compute genome-wide and locus-specific reference inbreeding in the offspring 
    Codes can be found in /referenceF directory

 4) Tranform the raw output of each method into genome-wide and locus-specific inbreeding predictions for each couple, and merge the predictions with the reference inbreeding
    Codes can be found in the /merge directory
 
 5) Compute genome-wide and locus-specific performance measures (correlations, bias, ROC curves)
    Codes are located in /performance directory

      * for performance evaluation (eval_performance.R) and visualization (visualize_performance.R)
	  * for the simulated data, test_significance.R computes the mean performance and CI across replicates and tests the significance of the difference in performance
	  
This workflow is illustrated using 1 replicate of the simulated data (from the moderate inbreeding scenario), but can be applied to real data as long as the input files required to run the programs are given in the same format

Considerations for real data are commented throughout the code where appropiate
