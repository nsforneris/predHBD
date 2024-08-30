import tskit
import msprime
import pyslim

ts = tskit.load("fwrd_parents.trees")

with open('out_pedigree.txt', 'w') as pedigree_file:
 for s in range(ts.num_individuals):
   ind = ts.tables.individuals[s]
   print(s, ind.parents[0], ind.parents[1], file = pedigree_file)


