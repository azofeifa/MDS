import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import linalg
def load_markov_chain_data(FILE):
	with open(FILE) as FH:
		i 	= -1
		L 	= [list(), list()]
		A 	= list()
		collect=False		
		for line in FH:
			if "Markov" in line:
				collect=True
				if A:
					L[i].append(np.array(A))
				A 		= list()
				i+=1
			elif collect and "$" in line:
				if A:
					L[i].append(np.array(A))
				A 	= list()
			elif collect:
				x 	= [float(u) for u in line.strip("\n").split("\t")[1].split(",")]
				A.append(x)
				pass
		if A:
			L[i].append(np.array(A))

	return L
def get_sums(L):
	X 	= list()
	for A in L:
		x 	= np.sum(A,axis=0)
		x 	/= sum(x)
		X.append(x)
	return X
def get_stationary(L):
	pi 	= np.array([0.25,0.25,0.25,0.25])
	X 	= list()
	for A in L:
		v = linalg.eig(A,left=True,right=False)[1][:,0]
		v = np.real(v / sum(v))
		if max(v) < 1.0:
			X.append(v)
		else:
			X.append(pi)

	return X




def plot_sum(L):
	F 		= plt.figure()
	ax1 	= F.add_subplot(1,2,1)
	ax2 	= F.add_subplot(1,2,2)
	axes 	= (ax1, ax2)
	for i,ax in enumerate(axes):
		X 	= get_stationary(L[i])
		for i in range(4):
			ax.plot(range(len(X)), [x[i] for x in X])

	plt.show()

if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/new_motif_distances/PSSM_DBS/new_Allen2014_generated_DB_2_MM.txt"
	L 		= load_markov_chain_data(FILE)
	plot_sum(L)