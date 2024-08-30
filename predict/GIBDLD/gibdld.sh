#!/bin/bash
#SBATCH --job-name=gibdld
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --partition=batch
#SBATCH --output=gibdld.out

progdir='/CECI/home/ulg/genan/forneris/programs/IBDLDv3.38.1'

# prepare input files to run GIBDLD method 
mfi="../input/sim_nontargets_gen.txt"
sfi="../input/sim_nontargets_samples.txt"
coup="../input/sim_couples"

# setup STUDY file with couples
cat $coup | awk '{print "coup", $1, $2}' > hya.study

# setup MAP file (sorted by chr and position)
cat $mfi | cut -d " " -f1-5 | awk '{print $1, $1"_"$3, $3/1000000, $3}' > hya.map

# setup GENOTYPE file
cat $mfi | cut -d " " -f6- > t1
cat $sfi | awk 'NR{print "coup",$1,"0","0","0"}' > t2
cat t1 | ../input/transpose.sh | sed 's/2/AA/g' | sed 's/1/AB/g' | sed 's/0/BB/g' | sed 's/AA/1 1/g' | sed 's/AB/1 2/g' | sed 's/BB/2 2/g' |sed 's/9/0 0/g' > t3
paste t2 t3 | awk '{$1=$1};1' > hya.ped
rm t1 t2 t3

# RUN GIBDLD
$progdir/ibdld -p hya.ped -m hya.map -s hya.study -method GIBDLD -ploci 10 -dist 2 -ibd 90 --ibdtxt -segment -ibd2segment -r 4000

rm *.ibd2segment *.segment *.primal *.gtype *.reg

gzip -9 prefix*

