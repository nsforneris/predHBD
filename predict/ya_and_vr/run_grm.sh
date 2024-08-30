#!/bin/bash

# directory of the program
progdir='/CECI/home/ulg/genan/forneris/programs'

# run GCTA for methods UNI (grm_ya equal to 2*avFUNI) & GRM (grm_vr) using the couple's binary files generated with plink
$progdir/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta64 --bfile '../input/sim_nontargets' --maf 0 --autosome-num 25 --autosome --make-grm --make-grm-alg 0 --make-grm-gz --out 'ya'
$progdir/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta64 --bfile '../input/sim_nontargets' --maf 0 --autosome-num 25 --autosome --make-grm --make-grm-alg 1 --make-grm-gz --out 'vr'


