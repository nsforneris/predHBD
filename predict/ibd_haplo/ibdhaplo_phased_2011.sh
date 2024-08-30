#!/bin/bash
#SBATCH --job-name=phatho11
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4092
#SBATCH --partition=batch
#SBATCH --output=ibdhaplo11p.out

ibdhaplodir="/CECI/home/ulg/genan/forneris/programs/MORGAN_V34_Release/IBD_Haplo"
cp $ibdhaplodir/Gold/phased_2011.par .

cat ../input/sim_couples            > coup
cp ../input/beagle.vcf.gz           phfi
tl=`zcat phfi | awk '$1=="#CHROM"{print NR}'`
zcat phfi | awk -v tl=$tl 'NR==tl' | ../input/transpose.sh | sed 1,9d > sfi

# 1. modify in .par file the name of the marker file and the gametes to be compared
  sed -i "s/input marker data file.*/input marker data file       \"\.\/data.markers\"/g" phased_2011.par
  # print proband gametes to compare
  cat coup | awk '{print NR, $1, $2}' > progam
  cat progam | awk '{printf "%s %s %s %s %s %s %s %s %s %s %s\n","set scoreset",$1, \
"proband gametes",$2,"0",$2,"1",$3,"0",$3,"1"}' > p2
  rm progam
  # replace gametes in par file
  cat phased_2011.par | sed '/^set scoreset/d' > tmp1
  awk '//; /# proband gametes/{while(getline<"p2"){print}}' tmp1 > tmp2
  rm tmp1 p2; mv tmp2 phased_2011.par

  sed -i "s/select markers.*/select all markers/g" phased_2011.par
  sed -i '/select all markers/q' phased_2011.par
  echo par file modified

# 2. set data.markers file for each chromosome and run ibd_haplo
for chr in {1..25}
 do
   zcat phfi | awk -v tl=$tl 'NR>tl' | awk -v v1=$chr '$1==v1 {$1=$1;print}' > t1
   # number of markers in chromosome
   nm=`wc -l t1 | awk '{print $1}'`
   # number of characters required for number of markers
   len=`wc -l t1 | awk '{print $1}' | awk '{print length}'`
   # marker positions (in cM)
   cat <(echo "map marker positions") <(cat t1 | awk '{printf "%.6f\n", $2/1000000}' | xargs -n20) > t2

   # marker frequencies
   cat t1 | awk '{$1=$1};1' | cut -d " " -f10- | ../input/transpose.sh > t3
   mv t3 t1
   for ((isnp=1;isnp<=$nm;isnp++))
    do
      cat t1 | awk -v col=$isnp '$col!=".|."{print $col}' | sed 's/|/ /g' | awk '{print ($1+$2)*(-1)+2}' > 1snp
      awk -v col=$isnp -v v1=$len 'BEGIN{sum=0}; {sum=sum+$1}; END{a=sum/(2*NR); \
      printf "%s %" v1 "s %s %.6f %.6f\n","set marker",col,"allele frequencies", a, 1-a}' 1snp >> t3
    done

   # marker genotypes
   cat t1 | sed 's/|/ /g' | sed 's/1/2/g' | sed 's/0/1/g' | sed 's/\./0/g' > t5
   cat t2 t3 <(echo "set markers" $nm "data") <(paste sfi t5) > data.markers

   rm t1 t2 t3 t5 1snp

   # run IBD_haplo by chromosome
   $ibdhaplodir/ibd_haplo phased_2011.par >> out_ibdhaplo_2011.txt

   # replace name of output files
   mv phased_2011.qibd  $chr'_phased_2011.qibd'
   mv phased_2011.ids   $chr'_phased_2011.ids'
   mv data.markers      $chr'_data.markers'
   mv $chr'_'*          output_phased/.

 done

echo finished running ibdhaplo
rm coup sfi phfi
mv out_ibdhaplo_2011.txt output_phased/.


