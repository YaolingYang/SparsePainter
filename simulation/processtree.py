import os
os.chdir("C:/Ubuntu/SLiM/sim_true/new5pop_1000ref2")

import subprocess, msprime
import matplotlib.pyplot as plt 
import numpy as np 
import tskit
import pandas as pd

def get_local_ancestry_pop(ts, pop, admixture_time, max_samples=100, per_rep=12):
	"""Return df describing local ancestry.

	Local ancestry is defined by the location of ancestors at admixture_time
	and does this by considering per_rep samples at once.
	"""
	ancestors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]
	pop_samples = ts.samples(population=pop)
	if max_samples:
		nsample = np.min([max_samples, pop_samples.size])  # number of samples to report
	else:
		nsample = pop_samples.size
	target_samples = pop_samples[:nsample]  # could allow random sampling of inds

	L = [x for x in range(0, nsample, per_rep)]
	R = [x for x in range(per_rep, nsample, per_rep)] + [nsample]
	assert(len(L) == len(R))
	dfs = []

	for i in range(len(L)):
		local = ts.tables.link_ancestors(
			samples=target_samples[L[i]:R[i]],
			ancestors=ancestors
		)
		local_df = pd.DataFrame({
			'left': local.left,
			'right': local.right,
			'parent': local.parent,
			'child': local.child
		})
		dfs.append(local_df)
		if i % 100 == 0:
			print(f'Done with local ancestry batch {i} out of {len(L)}!')

	local_ancestry_df = pd.concat(dfs)
	pop_of_node = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	# record the population of the ancestor, ie the local ancestry (localpop)
	local_ancestry_df['localpop'] = [pop_of_node[x] for x in local_ancestry_df['parent']]
	# record the population of the sample (samplepop)
	local_ancestry_df['samplepop'] = [pop_of_node[x] for x in local_ancestry_df['child']]
	local_ancestry_df = local_ancestry_df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	return(local_ancestry_df)

##ts = tskit.load("./ts.trees")
##a=get_local_ancestry_pop(ts, 3, 13, max_samples=1, per_rep=12)

## get true ancestries from 13-300 generations ago

ts = tskit.load("./ts.trees")

dfs = [] 

for i in range(13, 94, 1):
    print(i)
    a = get_local_ancestry_pop(ts, 5, i, max_samples=100, per_rep=12)
    
    a = a[["left", "right", "localpop","child"]]
    
    a["generation"] = i
    
    dfs.append(a)

a_combine = pd.concat(dfs, ignore_index=True)
a_combine.to_csv("sim_true.txt", sep="\t", index=False)
