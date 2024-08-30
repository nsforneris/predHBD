#!/bin/bash
#SBATCH --job-name=locngsrel
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --partition=batch
#SBATCH --output=locngsrel.out

progdir=/CECI/home/ulg/genan/forneris/programs/LocalNgsRelate/src/cpp

vcf='../input/plink.vcf.gz'
list='../input/sim_nontargets'
coup='../input/sim_couples'
freqfi='../input/allele_frequencies.txt'

# 1.Make a genotype likelihood file without targets (offspring) 
line=`zcat $vcf | awk '$1=="#CHROM"{print NR}'`
UNIQUE=`(awk '{print $1+9}' $list | tr '\n' ' ' | awk '{$1=$1};1' | sed 's/ /,/g')` #Array with desired field numbers (named $1 + 1st 9 fields in vcf file)
zcat $vcf | awk -v line=$line 'NR>=line' | awk '{$1=$1}1' | cut -d " " -f1-9,"${UNIQUE[@]}"  > last 

 head -1 last | sed -e 's/_//g' -e 's/#CHROM/marker/g' -e 's/REF ALT/allele1 allele2/g' | awk '{ for(i=10; i<=NF; i++) $i=$i"_"$i"_"$i } 1' | sed 's/\_/ /g' | cut -d " " -f1,4,5,10- > header
 # for the simulated data 
 cat last | awk '{$1=$1"_"$2}1' | cut -d " " -f1,4,5,10- | sed 's/|/\//g' | sed 's/0\/0/1 0 0/g' | sed 's/0\/1/0 1 0/g' |  sed 's/1\/0/0 1 0/g' |  sed 's/1\/1/0 0 1/g' | sed 's/\.\/\./0.333 0.333 0.333/g'  > ff
 # for the real data
 # cat last | awk '{$1=$1"_"$2}1' | cut -d " " -f1,4,5,10- | awk '{for(i=4;i<=NF;i++) {$i=substr($i,1,3)}OFS=" "; print $0}' | sed 's/|/\//g' | sed 's/0\/0/1 0 0/g' | sed 's/0\/1/0 1 0/g' |  sed 's/1\/0/0 1 0/g' |  sed 's/1\/1/0 0 1/g' | sed 's/\.\/\./0.333 0.333 0.333/g'  > ff

cat header <(awk 'NR>1' ff) > genolike.beagle; gzip genolike.beagle

# 2.Run LocalNgsRelate (it will remove monomorphic SNP)
head -1 last | cut -d " " -f10- | awk '{ for(i=1; i<=NF; i++) {print i-1, $i}1 }' > inds
join -1 1 -2 2 -o1.1,1.2,2.1 <(sort -k1,1 $coup) <(sort -k2,2 inds) > ff
join -1 2 -2 2 -o1.1,1.2,1.3,2.1 <(sort -k2,2 ff) <(sort -k2,2 inds) > coups
ncoup=`wc -l coups | awk '{print $1}'`
nnt=`wc -l $list | awk '{print $1}'`

for (( i=1; i<=$ncoup; i++ ))
do
  awk -v i=$i 'NR==i' coups > line
  id1=`cat line | awk '{print $1}'` 
  id2=`cat line | awk '{print $2}'`
  indx1=`cat line | awk '{print $3}'`
  indx2=`cat line | awk '{print $4}'`
  $progdir/localngsrelate -gbeagle genolike -f $freqfi -n $nnt -O localngsrelate -l 0 -r 05041983 -a $indx1 -b $indx2 
  zcat localngsrelate.IBDtractinference.gz | awk -v id1=$id1 -v id2=$id2 'NR>1{print id1, id2, $1"_"$2, $4, 0.25*$6+0.5*$7}' >> loc_locngsrelate.txt
done

rm line coups inds last ff header localngsrelate.IBDtractinference.gz


