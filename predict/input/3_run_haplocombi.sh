#!/bin/bash
#SBATCH --job-name=combi
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --partition=batch
#SBATCH --output=combi.out

# generates the four 2-by-2 combinations of parental haplotypes to be used by either RZooRoH or PLINK ROH for prediction 
# load required modules
module load R

phasedvcf='beagle.vcf.gz'
pairs='sim_ped'  

Rscript haplocombi.R $phasedvcf $pairs
