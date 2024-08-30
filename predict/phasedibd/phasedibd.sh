#!/bin/bash
#SBATCH --job-name=phasedibd
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --partition=batch
#SBATCH --output=phasedibd.out

# load required packages to run (and installed) phasedibd
module load releases/2021b
module load Python/3.9.6-GCCcore-11.2.0
module load SciPy-bundle/2021.10-foss-2021b

# git clone https://github.com/23andMe/phasedibd.git
# mkdir ~/my_venv
# virtualenv --system-site-packages ~/my_venv
# source ~/my_venv/bin/activate
#
# to compile: 
# cd phasedibd/
# make
# python setup.py install

# create sub-directory 'tests' in your working directory and place unit_tests.py with the arguments
# then, in your working directory
cd tests
mkdir -p resu data
# before the next step, copy shared resu directory if testing
rm -rf resu/* data/*

export PATH=/CECI/home/ulg/genan/forneris/programs/phasedibd:$PATH
export PYTHONPATH=/CECI/home/ulg/genan/forneris/programs/phasedibd/:$PYTHONPATH

vcffi='../../input/beagle.vcf.gz'
line=`zcat $vcffi | awk '$1=="#CHROM"{print NR}'`
zcat $vcffi | awk -v line=$line 'NR<=line' > first
zcat $vcffi | awk -v line=$line 'NR>line'  > last

# run phasedibd by chromosome (n=25)
for (( i=1; i<=25; i++ ))
do
  cat last | awk -v i=$i '$1==i' > lines
  cat first lines > data/myvcf.vcf
  python unit_tests.py
  cat resu/$i.csv >> phasedibd.segs
done

# recover original ID from the individuals for final output file
tail -1 first | awk '{$1=$1}1' | cut -d " " -f10- | awk '{ for(i=1; i<=NF; i++) {print i-1, $i}1 }' > targets
join -1 1 -2 1 -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2 <(cat phasedibd.segs | sed 's/,/ /g' | sort -k1,1) <(sort -k1,1 targets) > ff
join -1 2 -2 1 -o1.12,2.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11 <(sort -k2,2 ff) <(sort -k1,1 targets) | sort -u  > ../loc_phasedibd.txt

# delete temporal files (partial outputs still inside the resu/ directory)
rm -r first last lines targets ff phasedibd.segs data

module load R
Rscript post_run.R 



