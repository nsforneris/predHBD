import tskit
import msprime
import numpy as np
import sys

semilla = sys.argv[1]
print(semilla)

ts = tskit.load("recap_w_pyslim.trees")

# save age of every node to an output
with open('out_nodeage.txt', 'w') as node_age_file:
 for nod_id in range(ts.num_nodes):
  nod_age = ts.nodes_time[nod_id]
  print(nod_id, nod_age, file = node_age_file)

# save individual ID and parents ID of every sampled node to an output
samp_ids = ts.samples()
with open('out_nodefamily.txt', 'w') as node_fam_file:
 for s in samp_ids:
   indnod = ts.node(s).individual
   parents = ts.individual(indnod).parents
   print(s, indnod, parents[0], parents[1], file = node_fam_file)

with open('out_segments.txt', 'w') as ibd_file:
 individuals = ts.individuals()

## print HBD segments for the first 20 individual IDs
 for s in range(20):
    ind=individuals[s]
    bb0=ind.nodes[0]
    bb1=ind.nodes[1]
    segments = ts.ibd_segments(within=[bb0, bb1], store_segments=True, store_pairs=True)
    for pair, segment_list in segments.items():
        for ibdseg in segment_list:
          print(s, pair[0], pair[1], ibdseg.left, ibdseg.right, ibdseg.node, file = ibd_file)

## overlay mutations
mutated = msprime.sim_mutations(ts, rate=1e-8, random_seed=semilla, keep=True) 
mutated.dump("recap_overlaid.trees")

# save variant info (position and ID start from 0 to genome length -1 and nvariants -1, respectively) 
with open('out_variant_info.txt', 'w') as variants_file:
 for var in mutated.variants():
  print(var.site.position, var.site.id, var.alleles, file = variants_file)

with open("output.vcf", "w") as vcf_file:
    mutated.write_vcf(vcf_file)







