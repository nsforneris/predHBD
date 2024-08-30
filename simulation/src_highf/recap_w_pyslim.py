import tskit
import msprime
import pyslim
import sys

semilla = sys.argv[1]
print(semilla)

# load the .trees file generated with SLiM
orig_ts = tskit.load("fwrd.trees")

# adapt recombination map to msprime
positions = []
rates = []
with open('recomb_rates.tsv', 'r') as file:
  header = file.readline().strip().split("\t")
  assert(header[0] == "end_position" and header[1] == "rate(cM/Mb)")
  for line in file:
     components = line.split("\t")
     positions.append(float(components[0]))
     rates.append(1e-8 * float(components[1]))

# step 1
positions.insert(0, 0)
# step 2
positions[-1] += 1
assert positions[-1] == orig_ts.sequence_length

recomb_map = msprime.RateMap(position=positions, rate=rates)

# recapitate
rts = pyslim.recapitate(orig_ts,
                ancestral_Ne=10000,
                recombination_rate=recomb_map,
                model=[
                    msprime.DiscreteTimeWrightFisher(duration=20),
                    msprime.StandardCoalescent(),
                    ],
                random_seed=semilla)

print("finish recap")

assert(max([t.num_roots for t in rts.trees()]) == 1)
orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

rts.dump("recap_w_pyslim.trees")

