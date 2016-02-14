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
		vls 	= [np.random.multinomial(1,background,1).argmax() for k in range(length)] 
		seq 	= "".join([str(v) for v in vls ])
		if seq not in G:
			lls.append(2*((sum([  math.log(A[k,j])  for k,j in enumerate(vls) ]) )))
			G[seq] 	= 0
	counts,edges 	= np.histogram(lls,bins=500)
	edges 			= edges[1:]
	return edges,counts
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
def find_5(A, background, edges, counts,N=50000):
	NN 	= float(sum(counts))
	print "-----------------------------"
	print background, 
	length 	= A.shape[0]
	t 		= 0
	MAP 	= {0:"A",1:"C",2:"G",3:"T"}
	scores 	= list()

	for i in range(N):
		vls 	= [np.random.multinomial(1,[0.25,0.25,0.25,0.25],1).argmax() for k in range(length)] 
		seq 	= "".join([MAP[v] for v in vls ])
		score 	= 2*((sum([  math.log(A[k,j])  for k,j in enumerate(vls) ]) ))
		S 		= 0.0
		for j in range(0,len(edges)):
			if score < edges[j]:
				break
			S+=counts[j]
	
		if S/NN > 0.99:
			scores.append(sum([  math.log(A[k,j])  for k,j in enumerate(vls) ]) )
#			print seq, (sum([  math.log(A[k,j])  for k,j in enumerate(vls) ]) )
			t+=1
	print len(scores)
	print 
	return scores

def across_all(P,background=[0.25,0.25,0.25,0.25]):
	#check to see if background alread exists
	src 	=  os.path.dirname(os.path.realpath(__file__))
	backgrounds  	= ([0.25,0.25,0.25,0.25], [0.15,0.35,0.35,0.15],[0.35,0.15,0.15,0.35])
	for m in P:
	#	print m
		A 			= P[m]
		if "CTCF" in m:
			F 			= plt.figure()
			ax 			= F.add_subplot(1,2,1)
			ax2 		= F.add_subplot(1,2,2)
			colors 		= ("r", "g", "b")

			for  i,background in enumerate(backgrounds):
#				vls,counts 	= get_pvals(A,background=background)
				vls, counts = compute_null(A, background=background)
				scores 		= find_5(A,background, vls,counts)
				ax.bar(vls, counts,label=str(background),alpha=0.2, color=colors[i])
				ax2.hist(scores, bins=30, color=colors[i],alpha=0.3)
			ax.legend(loc="best")
			plt.show()
			counts 		= counts + [counts[-1]]
			N 			= float(sum(counts))


if __name__ == "__main__":
	PSSM 	= "/Users/joazofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
	P 		= load.load_PSSMs(PSSM)
	across_all(P)
	#iterate_over_all(P)
	# get_pvals(A)
