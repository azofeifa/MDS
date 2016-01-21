import load
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
def compute_null(A,MC=False, N=10000, background=[0.25,0.25,0.25,0.25]):
	length 	= A.shape[0]
	lls 	= list()
	G 		= {}
	for i in range(N):
		vls 	= [np.random.multinomial(1,[0.25,0.25,0.25,0.25],1).argmax() for k in range(length)] 
		seq 	= "".join([str(v) for v in vls ])
		if seq not in G:
			lls.append(2*((sum([  math.log(A[k,j])  for k,j in enumerate(vls) ]) )-sum([ math.log(background[j])  for j in vls ]  )))
			G[seq] 	= 0
	return lls
def iterate_over_all(P):
	for motif in P:
		A 	= P[motif]
		compute_null(A)
def get_pvals(A,background=[0.25,0.25,0.25,0.25],bins=100 ):
	vls 	= list()
	counts,edges 	= None, None
	min_score 		= 2*sum([ min([math.log(A[i,j]) - math.log(background[j]) for j in range(A.shape[1])]) for i in range(A.shape[0])  ])
	max_score 		= 2*sum([ max([math.log(A[i,j]) - math.log(background[j]) for j in range(A.shape[1])]) for i in range(A.shape[0])  ])

	for i in range(A.shape[0]):
		nvls 	= list()
		weights = list()
		for j in range(A.shape[1]):
			if i==0:
				nvls.append( 2*(math.log(A[i,j])-math.log(background[j]))   )
			else:
				for k,vl in enumerate(edges):
					if counts[k]:
						nvls.append( vl+ 2*(math.log(A[i,j])-math.log(background[j])) )
						weights.append(counts[k])
		if i==0:
			counts,edges 	= np.histogram(nvls,bins=bins,range=(min_score, max_score))
			edges 			= edges[:-1]
		else:
			counts,edges 	= np.histogram(nvls,bins=bins, weights=weights,range=(min_score, max_score))
			edges 			= edges[:-1]
	return edges, counts
def across_all(P,background=[0.25,0.25,0.25,0.25]):
	#check to see if background alread exists
	src 	=  os.path.dirname(os.path.realpath(__file__))
	pval_dir = "/".join(src.split("/")[:-1])+ "/exact_pvalues/"
	FILE 	=pval_dir + "_".join([str(i) for i in background ]) + ".tsv"
	if not os.path.exists(FILE):
		FHW 	= open(FILE, "w")
		for m in P:
			print m
			A 			= P[m]
			vls,counts 	= get_pvals(A,background=background)
			counts 		= counts + [counts[-1]]
			N 			= float(sum(counts))
			FHW.write(m+"\t" + ",".join([str(x) + "_" + str(float(sum(counts[:i]))/N )for i,(x,y) in enumerate(zip(vls, counts))]) + "\n")
	pass


if __name__ == "__main__":
	PSSM 	= "/Users/joazofeifa/Lab/gTFI/files/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
	P 		= load.load_PSSMs(PSSM)
	across_all(P)
	#iterate_over_all(P)
	# get_pvals(A)
