#!/bin/bash
#SBATCH --job-name=tho11
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4092
#SBATCH --partition=batch
#SBATCH --output=ibdhaplo11.out

ibdhaplodir="/CECI/home/ulg/genan/forneris/programs/MORGAN_V34_Release/IBD_Haplo"
cp $ibdhaplodir/Gold/unphased_2011.par .

cat "../input/sim_couples"                 > coup
cat "../input/sim_nontargets_gen.txt"      > mfi
cat "../input/sim_nontargets_samples.txt"  > sfi

# 1. modify in .par file the name of the marker file and the gametes to be compared
  sed -i "s/input marker data file.*/input marker data file       \"\.\/data.markers\"/g" unphased_2011.par
  # print proband gametes to compare
  cat coup | awk '{print NR, $1, $2}' > progam
  cat progam | awk '{printf "%s %s %s %s %s %s %s %s %s %s %s\n","set scoreset",$1, \
"proband gametes",$2,"0",$2,"1",$3,"0",$3,"1"}' > p2
  rm progam
  # replace gametes in par file
  cat unphased_2011.par | sed '/^set scoreset/d' > tmp1
  awk '//; /# proband gametes/{while(getline<"p2"){print}}' tmp1 > tmp2
  rm tmp1 p2; mv tmp2 unphased_2011.par

  sed -i "s/select markers.*/select all markers/g" unphased_2011.par
  sed -i '/select all markers/q' unphased_2011.par
  echo par file modified

# 2. set data.markers file for each chromosome and run ibd_haplo
for chr in {1..25}
 do
   cat mfi | awk -v v1=$chr '$1==v1' > t1
   # number of markers in chromosome
   nm=`wc -l t1 | awk '{print $1}'`
   # number of characters required for number of markers
   len=`wc -l t1 | awk '{print $1}' | awk '{print length}'`
   # marker positions (in cM)
   cat <(echo "map marker positions") <(cat t1 | awk '{printf "%.6f\n", $3/1000000}' | xargs -n20) > t2
   # marker allele frequencies
   cat t1 | awk '{$1=$1};1' | cut -d " " -f6- | ../input/transpose.sh > t3
   mv t3 t1
   for ((isnp=1;isnp<=$nm;isnp++))
    do
      awk -v col=$isnp '$col!=9{print $col}' t1 > 1snp
      awk -v col=$isnp -v v1=$len 'BEGIN{sum=0}; {sum=sum+$1}; END{a=sum/(2*NR); \
      printf "%s %" v1 "s %s %.6f %.6f\n","set marker",col,"allele frequencies", a, 1-a}' 1snp >> t3
    done
   # marker genotypes
   cat t1 | sed 's/2/AA/g' | sed 's/1/AB/g' | sed 's/0/BB/g' \
   | sed 's/AA/1 1/g' | sed 's/AB/1 2/g' | sed 's/BB/2 2/g' | sed 's/9/0 0/g' > t5
   cat t2 t3 <(echo "set markers" $nm "data") <(paste sfi t5) > data.markers

   rm t1 t2 t3 t5 1snp

   # run IBD_haplo by chromosome
   $ibdhaplodir/ibd_haplo unphased_2011.par >> out_ibdhaplo_2011.txt

   # replace name of output files
   mv unphased_2011.qibd $chr'_unphased_2011.qibd'
   mv unphased_2011.ids  $chr'_unphased_2011.ids'
   mv data.markers       $chr'_data.markers'
   mv $chr'_'*           output_unphased/.

 done

echo finished running ibdhaplo
rm coup mfi sfi
mv out_ibdhaplo_2011.txt output_unphased/.




