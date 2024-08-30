#!/bin/bash
#SBATCH --job-name=beagle
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --partition=batch
#SBATCH --output=beagle.out

progdir='/CECI/home/ulg/genan/forneris/programs/beagle'

# load required modules to run Beagle
module load Java/1.8.0_281

vcf='plink.vcf.gz'
list='sim_nontargets' # list without the offspring 

# Make a new vcf called vcf.vcf excluding the offspring (target animals) (an alternative of how to do it)  
line=`zcat $vcf | awk '$1=="#CHROM"{print NR}'`
zcat $vcf | awk -v line=$line 'NR<line' > first
UNIQUE=`(awk '{print $1+9}' $list | tr '\n' ' ' | awk '{$1=$1};1' | sed 's/ /,/g')` #Array with desired field numbers (named $1 + 1st 9 fields in vcf file)
zcat $vcf | awk -v line=$line 'NR>=line' | awk '{$1=$1}1' | cut -d " " -f1-9,"${UNIQUE[@]}" | awk '{if(NR>1)$3=$1"_"$2};1' | sed 's/ /\t/g' > last 
cat first last > vcf.vcf
rm first last

# Create folder temporal files:
ScratchDir='TMP'
mkdir -p ${ScratchDir}

# Run Beagle
java -Xmx4000m -Djava.io.tmpdir=${ScratchDir} -jar $progdir/beagle.25Mar22.4f6.jar seed=050483 nthreads=2 gt=vcf.vcf out=beagle

# Delete temporal files:
rm -r ${ScratchDir}
rm vcf.vcf
