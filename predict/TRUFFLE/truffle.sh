#!/bin/bash
#SBATCH --job-name=truffle
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --partition=batch
#SBATCH --output=truffle.out

progdir=/CECI/home/ulg/genan/forneris/programs/truffle

vcf=../input/plink.vcf.gz
list=../input/sim_nontargets

# Make a new vcf called truffle.vcf excluding the offspring (target animals) 
line=`zcat $vcf | awk '$1=="#CHROM"{print NR}'`
UNIQUE=`(awk '{print $1+9}' $list | tr '\n' ' ' | awk '{$1=$1};1' | sed 's/ /,/g')` #Array with desired field numbers (named $1 + 1st 9 fields in vcf file)
zcat $vcf | awk -v line=$line 'NR>=line' | awk '{$1=$1}1' | cut -d " " -f1-9,"${UNIQUE[@]}" | awk '{if(NR>1)$3=$1"_"$2};1' > last 

# extra cleaning specific for truffle.vcf 
zcat $vcf | awk -v line=$line 'NR<line' | awk '!($1 ~ /##contig/)' | awk '!($1 ~ /##INFO/)' > first
# for the simulated vcf
cat last | awk '{if(NR>1)$8="."};1' | awk '{if(NR>1)$9="GT"};1' | sed 's/|/\//g' > ff
 # for the real vcf (just keep the genotype (first three characters)
 # cat last | awk '{if(NR>1)$8="."};1' | awk '{if(NR>1)$9="GT"};1' | awk '{for(i=10;i<=NF;i++) {$i=substr($i,1,3)}OFS=" "; print $0}' | sed 's/|/\//g' > ff
cat first ff > truffle.vcf

# Run truffle 
$progdir/truffle --cpu 1 --vcf truffle.vcf --segments --maf -1 --missing 1 --mindist -1 --L 2 --ibs1markers 100 --ibs2markers 50 --out truffle

# Delete temporal files
rm first ff last 




