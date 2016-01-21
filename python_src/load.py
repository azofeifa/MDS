import re
import numpy as np
def load_PSSMs(FILE):
	P 	= {}
	MOTIF 	= None
	with open(FILE) as FH:
		for line in FH:
			if line[:5] == "MOTIF":
				line_array 	= re.split("\s+",line)
				MOTIF 		= line_array[1].split("_")[0]
				P[MOTIF] 	= [0,[]]
			elif line[:6] == "letter":
				line_array 	= re.split("\s+",line)
				P[MOTIF][0] = float(line_array[7])
			elif MOTIF is not None and line[:3]!= "URL":
				line_array 	= re.split("\s+", line.strip("\n"))
				P[MOTIF][1].append([float(x) for x in line_array])
			elif line[:3]=="URL":
				MOTIF 	= None
	U 	= {}
	for MOTIF in P:
		A 				= np.array(P[MOTIF][1])
		A[:,:]*=P[MOTIF][0]
		A[:,:]+=1
		for i in range(A.shape[0]):
			A[i,:]/=sum(A[i,:])

		U[MOTIF]=A
	return U
def load_bed_file(FILE, union=0):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if line[0]!= "#":
				chrom,start, stop 	= line.strip("\n").split("\t")[:3]
				if chrom not in G:
					G[chrom] 		= list()
				if union == 0:
					G[chrom].append((start, stop))
				else:
					x 		= (float(start) + float(stop))/2.
					G[chrom].append([x-union, x+union, ""])
	for chrom in G:
		G[chrom].sort()
	return G








if __name__ == "__main__"	:
	PSSM 	= "/Users/joazofeifa/Lab/gTFI/files/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
	P 		= load_PSSMs(PSSM)
