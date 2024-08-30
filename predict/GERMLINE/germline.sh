#!/bin/bash
#SBATCH --job-name=germline
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --partition=batch
#SBATCH --output=germline.out

export TMP=/tmpscratch
export TEMP=/tmpscratch
export TMPDIR=/tmpscratch
progdir=/CECI/home/ulg/genan/forneris/programs/germline-1-5-3
 
vcffi=../input/beagle.vcf.gz
line=`zcat $vcffi | awk '$1=="#CHROM"{print NR}'`
zcat $vcffi > vcffi

# for 25 autosomes
for (( i=1; i<=25; i++ ))
  do
    # make basic .map and .ped files per chromosome
    cat vcffi | awk -v line=$line 'NR>line{print $1,$1"_"$2,$2/1000000,$2}' | awk -v i=$i '$1==i' > germline.map
    cat vcffi | awk -v line=$line 'NR>=line' | awk '{$1=$1}1' | awk -v i=$i '$1==i || $1=="#CHROM"' | cut -d " " -f10- | ../input/transpose.sh | sed 's/0|0/1 1/g' | sed 's/0|1/1 2/g' | sed 's/1|0/2 1/g' | sed 's/1|1/2 2/g' | awk '{$1="FAM "$1" 0 0 0 -9"}1' > germline.ped

    # run germline
    $progdir/germline -input germline.ped germline.map -min_m 2 -haploid -output germline_$i

    # split field 2 into "individual1's ID" and "individual1's shared haplotype"; same for field 4
    # and concatenate results for each chromosome 
    cat germline_$i.match | awk 'BEGIN{}{gsub(/\./, " ", $2)} 1' | awk 'BEGIN{}{gsub(/\./, " ", $5)} 1' >> loc_germline_haplo.txt
  done

# Delete temporal files
rm vcffi germline.map germline.ped

