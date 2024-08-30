#!/bin/bash

# script to make the input files required to run the prediction methods compared in our study
# it uses the output files from replicate #100 of the moderate inbreeding scenario 
   /out_nodefamily.txt.gz 
   /plink.vcf.gz   #genotypes of the animals from the last 2 generations (250 offspring + 250 parents) 

# 1. make RZooRoH input files (gen and sample files) for the 500 genotyped animals
   # sim_gen.txt counts A1 alleles (0(ALT ALT) 1(REF ALT) 2(REF REF))
   # sim_samples.txt contains individual's ID (a character string) as 'tsk_'val, where val is (the order of appearance in plink.vcf.gz - 1)

   vcffile='plink.vcf'
   num=`cat $vcffile | awk '$1=="#CHROM"{print NR}'`
   cat $vcffile | awk -v num=$num 'NR==num' | awk '{$1=$1}1' | cut -d " " -f10- | ./transpose.sh > 'sim_samples.txt'
   cat $vcffile | awk -v num=$num 'NR>num' | awk '{$6="";$7="";$8="";$9=""}1' | tr -s " " | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' |  sed 's/1\/0/1/g' |  sed 's/1\/1/0/g' | sed 's/\.\/\./9/g' | sed 's/chr//'| awk '{t=$2; $2=$3;$3=t; print;}' | awk '{$2=$1"_"$3; print;}'  > 'sim_gen.txt'
   gzip -9 sim_gen.txt
   
# 2. create useful files:
   # a 'targets' list (here the first 100 offspring from last generation)
   # a couples' list
   # a 'non-targets' list (parents of the targets)
   # a 2-generation-pedigree file linking targets to parents 
   
   seq 0 499 | awk '{print NR, "tsk_"$1}' > sim_alive
   # animals with 0<=val<=99 correspond to targets used in our study
   head -n 100 sim_alive > sim_targets 
   # the parents of the targets are extracted from the nodes' info
   zcat out_nodefamily.txt.gz | awk '$2>=0 && $2<=99{print $2,$3,$4}' | sort -u | sort -k1,1n | awk '{print "tsk_"$1, "tsk_"$2, "tsk_"$3}' > sim_ped
   cat sim_ped | awk '{print $2, $3}' | sort -u > sim_couples
   cat <(awk '{print $2}' sim_ped) <(awk '{print $3}' sim_ped) | sort -u > ff
   join -1 1 -2 2 -o2.1,2.2 <(sort -k1,1 ff) <(sort -k2,2 sim_alive) | sort -k1,1n > sim_nontargets
   
# 3. create RZooroh "gen" and "sample" files for the 'non-targets'
   fullfile='sim_gen.txt.gz'
   list='sim_nontargets'
   smallfile='sim_nontargets_gen.txt'
   cat $list | awk '{print $2}' > 'sim_nontargets_samples.txt'
   UNIQUE=`(awk '{print $1+5}' $list | tr '\n' ' ' | awk '{$1=$1};1' | sed 's/ /,/g')`  # Array containing fields (named $1 + first 5 cols) to be printed
   zcat $fullfile | cut -d " " -f1-5,"${UNIQUE[@]}"  >  $smallfile

   # you can now estimate allele frequencies in the parents by running ZR_allelefreq.R inside ../RZooRoH directory

# 4. make PLINK .bim .bed and .fam files (with same REF-ALT assignment of alleles as plink.vcf.gz file)
   # by default, '--vcf --make-bed' sets REF alleles to "A2" by the autoconverter which are assigned to 2nd column in .bim
   # to force the REF allele to be "A1":
   progdir='/CECI/home/ulg/genan/forneris/programs/' 

   # a. for all animals
   vcffile='plink.vcf.gz'
   line=`zcat $vcffile | awk '$1=="#CHROM"{print NR}'`
   zcat $vcffile | awk -v line=$line 'NR>line {print $1"_"$2,$4,$5}'  > ref_alt
   $progdir/plink --memory 4000 --vcf $vcffile --const-fid --chr-set 25 no-xy --make-bed --out vcf
   cat vcf.bim | awk '{$2=$1"_"$4}1' | column -t > tmp
   mv tmp vcf.bim
   $progdir/plink --memory 4000 --chr-set 25 --bfile vcf --make-bed --keep-allele-order --a1-allele ref_alt 2 1 --out sim

   # b. for the non targets 
   bigfile='sim'
   short='sim_nontargets_samples.txt'
   awk '{print 0, $1}' $short > keepfile
   $progdir/plink --memory 4000 --chr-set 25 --bfile $bigfile --keep-allele-order --keep keepfile --make-bed --out "sim_nontargets"

# 5. remove temporary files
   rm ff sim_alive keepfile *nosex vcf.bed vcf.bim vcf.fam vcf.log
   


