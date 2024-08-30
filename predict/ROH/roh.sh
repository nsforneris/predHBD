#!/bin/bash
#SBATCH --job-name=roh
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4092
#SBATCH --partition=batch
#SBATCH --output=roh.out

progdir='/CECI/home/ulg/genan/forneris/programs'
samplefi='../input/haplocombi_samples.txt'

# tped and tfam format
cat $samplefi | awk '{print $1, $1, "0", "0", "0", "-9"}' > fake.tfam

# pat-pat combination
 genfi='../input/haplocombi_pat_pat_gen.txt'
 cat $genfi | awk '{print $1, $2, "0", $3}' > f0
 cat $genfi | cut -d " " -f6-  | sed 's/0/B B/g' | sed 's/1/A B/g' | sed 's/2/A A/g' > f1
 paste f0 f1 | awk '{$1=$1}1' > fake.tped
 rm f0 f1
 # ROH with 2 mb of min length, >40 snp, inverse density of 100 kb/snp, 500 kb gap and no heterozigous snps
 $progdir/plink --memory 4000 --chr-set 25 --tfile fake --keep-allele-order --homozyg-kb 2000 --homozyg-snp 40 --homozyg-density 100 --homozyg-gap 500 --homozyg-het 0 --homozyg-window-snp 40 --homozyg-window-het 0 --out fake_proh_pp

# pat-mat combination
 genfi='../input/haplocombi_pat_mat_gen.txt'
 cat $genfi | awk '{print $1, $2, "0", $3}' > f0
 cat $genfi | cut -d " " -f6-  | sed 's/0/B B/g' | sed 's/1/A B/g' | sed 's/2/A A/g' > f1
 paste f0 f1 | awk '{$1=$1}1' > fake.tped
 rm f0 f1
 $progdir/plink --memory 4000 --chr-set 25 --tfile fake --keep-allele-order --homozyg-kb 2000 --homozyg-snp 40 --homozyg-density 100 --homozyg-gap 500 --homozyg-het 0 --homozyg-window-snp 40 --homozyg-window-het 0 --out fake_proh_pm

# mat-pat combination
 genfi='../input/haplocombi_mat_pat_gen.txt'
 cat $genfi | awk '{print $1, $2, "0", $3}' > f0
 cat $genfi | cut -d " " -f6-  | sed 's/0/B B/g' | sed 's/1/A B/g' | sed 's/2/A A/g' > f1
 paste f0 f1 | awk '{$1=$1}1' > fake.tped
 rm f0 f1
 $progdir/plink --memory 4000 --chr-set 25 --tfile fake --keep-allele-order --homozyg-kb 2000 --homozyg-snp 40 --homozyg-density 100 --homozyg-gap 500 --homozyg-het 0 --homozyg-window-snp 40 --homozyg-window-het 0 --out fake_proh_mp

# mat-mat combination
 genfi='../input/haplocombi_mat_mat_gen.txt'
 cat $genfi | awk '{print $1, $2, "0", $3}' > f0
 cat $genfi | cut -d " " -f6-  | sed 's/0/B B/g' | sed 's/1/A B/g' | sed 's/2/A A/g' > f1
 paste f0 f1 | awk '{$1=$1}1' > fake.tped
 rm f0 f1
 $progdir/plink --memory 4000 --chr-set 25 --tfile fake --keep-allele-order --homozyg-kb 2000 --homozyg-snp 40 --homozyg-density 100 --homozyg-gap 500 --homozyg-het 0 --homozyg-window-snp 40 --homozyg-window-het 0 --out fake_proh_mm

rm *hom.summary *.nosex

