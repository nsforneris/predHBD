#!/bin/bash
#SBATCH --job-name=beagle
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4092
#SBATCH --partition=batch
#SBATCH --output=hapibd.out

progdir=/CECI/home/ulg/genan/forneris/programs/beagle

vcffi=../input/beagle.vcf.gz
line=`zcat $vcffi | awk '$1=="#CHROM"{print NR}'`
# make map file
zcat $vcffi | awk -v line=$line 'NR>line{print $1,$1"_"$2,$2/1000000,$2}' > map

# Create folder temporal files (optional)
ScratchDir='TMP'
mkdir -p ${ScratchDir}
# Load required modules to run hap-IBD
module load Java/1.8.0_281

# run hap-IBD
java -Xmx4g -Djava.io.tmpdir=${ScratchDir} -jar $progdir/hap-ibd.jar nthreads=1 gt=$vcffi map=map out=hapibd

# Delete temporal files
rm -r ${ScratchDir}
rm map

