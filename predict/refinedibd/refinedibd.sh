#!/bin/bash
#SBATCH --job-name=refined
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4092
#SBATCH --partition=batch
#SBATCH --output=refined.out

progdir=/CECI/home/ulg/genan/forneris/programs/beagle
# Create folder temporal files (optional)
ScratchDir=TMP
mkdir -p ${ScratchDir}
# Load modules required to run Refined IBD
module load Java/1.8.0_281

vcffi=../input/beagle.vcf.gz
line=`zcat $vcffi | awk '$1=="#CHROM"{print NR}'`
# make map file
zcat $vcffi | awk -v line=$line 'NR>line{print $1,$1"_"$2,$2/1000000,$2}' > map

# run Refined IBD
java -Xss5m -Xmx4g -Djava.io.tmpdir=${ScratchDir} -jar $progdir/refined-ibd.17Jan20.102.jar nthreads=1 gt=$vcffi map=map out=refinedibd

# Merge reported IBD segments
zcat refinedibd.ibd.gz | java -jar $progdir/merge-ibd-segments.17Jan20.102.jar $vcffi map 0.5 1 > refinedibd.ibd.merged
gzip -9 refinedibd.ibd.merged

# Delete temporal files
rm -r ${ScratchDir}
rm map

#module load R
#Rscript /scratch/ulg/genan/forneris/SLiMtest/predict_modf/refinedibd/frefinedibd.R /scratch/ulg/genan/forneris/SLiMtest/modf$repli/sim_ped /scratch/ulg/genan/forneris/SLiMtest/modf$repli/sim_targets /scratch/ulg/genan/forneris/SLiMtest/modf$repli/predict sim_nontargets

