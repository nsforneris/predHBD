import numpy as np
import pandas as pd
import os
import sys
import unittest
import phasedibd as ibd

TEST_DATA_PATH = os.path.dirname(os.path.abspath(__file__)) + '/data/'


haplotypes = ibd.VcfHaplotypeAlignment(TEST_DATA_PATH + 'myvcf.vcf')
tpbwt = ibd.TPBWTAnalysis(template=[[1, 0, 1, 0],
 [0, 1, 0, 1],
 [1, 1, 0, 0],
 [0, 0, 1, 1],
 [1, 0, 0, 1],
 [0, 1, 1, 0]])
ibd_results = tpbwt.compute_ibd(haplotypes, L_m=50, L_f=2, use_phase_correction=False, segments_out_path='resu/')



