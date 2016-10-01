import sys
import numpy as np
import math
from scipy.special import erf
import matplotlib.pyplot as plt
from scipy.spatial import KDTree 

class GRID:
	def __init__(self, x1, x2, y1,y2, bins):
		self.max_y 	= y2
		self.max_x 	= x2
		self.grid 	= np.array([[x/x2,y/y2,0] for x in np.linspace(x1, x2,bins) for y in np.linspace(y1,y2,bins)])
		self.bins 	= bins
		self.tree 	= KDTree(self.grid)
		self.res_x 	= 1.0 #math.sqrt((x2-x1))
		self.res_y 	= 1.0 #math.sqrt((y2-y1)/bins)/1000.0
		self.T 		= 100
	def find_nearest_available(self, x):
		x[0,0]/=self.max_x
		x[0,1]/=self.max_y
		val, arg 	=  self.tree.query(x[:])
		if self.grid[arg,2] < 1:
			self.grid[arg,2]=1
			x1,y1 	=self.grid[arg][0,:2] 
			return x1*self.max_x,y1*self.max_y
		else:
			t 		= 0
			while t < self.T:
				x[0,0]=np.random.normal(x[0,0],self.res_x)
				x[0,1]=np.random.normal(x[0,1],self.res_y)
				val, arg 	=  self.tree.query(x)
				if self.grid[arg,2]==0:
					self.grid[arg,2]=1

					x1,y1 	=self.grid[arg][0,:2]
					return x1*self.max_x,y1*self.max_y
				t+=1
		x1,y1 	=self.grid[arg][0,:2]
		return x1*self.max_x,y1*self.max_y
	def insert_data_points(self, X):
		args 	= {}
		for x,y in X:
			xx 			= np.array([(x/self.max_x,y/self.max_y,0)])	
			val, arg 	=  self.tree.query(xx[:])
			args[arg[0]] 	= 1
		for arg in args:
			self.grid[arg, 2]=1
class PSSM:

	def __init__(self,motif,name, N, MD_score,  pv_s,pv_ns, null ):
		self.full 		= motif
		self.name 		= name
		self.N 			= N
		self.MD_score 	= MD_score
		self.pv_ns 		= pv_ns
		self.pv_s 		= pv_s
		self.null 		= null

def load_enrichment_file(FILE):
	P 	= {}
	with open(FILE) as FH:
		for line in FH:
			if "#" == line[0] and "Binned" in line:
				break
			if line[0]!="#":
				line_array 						= line.strip("\n").split("\t")
				n,MD, pv_s, pv_ns, null 		= [ [float(n) for n in line_array[i].split(",")]  for i in range(1,len(line_array))]
				p 								= PSSM( line_array[0],line_array[0].split("_")[0], n, MD, pv_s, pv_ns, null )
				if p.full in P:
					P[p.full] 					= max((P[p.name].MD_score[0],P[p.name]), (p.MD_score[0], p))[1]
				else:
					P[p.full] 					= p 
	return P
def get_SE(p1,p2,n1, n2, PSSM, M, cutoff=pow(10,-5.25), threshold=0.049):
	p 		= (p1*n1 + p2*n2) / float( n1 + n2  )

	SE 		= math.sqrt( p*(1-p)*( (1.0/ n1 ) + (1.0/n2 ) ))

	Z 		= (p2-p1 - M) / SE
	pval 	= 0.5*(1.0 + erf(Z / math.sqrt(2.0) ))
	size 	= 10.0
	color 	= "blue"
	if pval < cutoff and abs(p2-p1) > threshold:
		color 	= "green"
		size 	= 70.0
	elif pval > (1.0 - cutoff) and abs(p2-p1) > threshold :
		color 	= "red"
		size 	= 70.0
	return  p2-p1-M, SE, min(n1,n2),pval, color,size,PSSM.split("_")[0]
def compare(F1, F2, title, ax=None, cutoff=pow(10,-5), threshold=0.049, u=0):
	
	PSSMS 	= [p for p in F1 if p in F2 and F1[p].N[u] > 10 and F2[p].N[u]>10 and "CPEB1" not in p and "PPARA" not in p ]
	M 		= np.mean([ F2[p].MD_score[u] - F1[p].MD_score[u] for p in F1 if p in F2])
	
	stats 	= [ get_SE(F1[p].MD_score[u],F2[p].MD_score[u],F1[p].N[u], F2[p].N[u],p, M,cutoff=cutoff,threshold=threshold) for p in PSSMS]
	SS 		= False
	if ax is None:
		ax 	= plt.gca()
		SS 	= True

	
	ax.scatter( [ s[2] for s in stats], [s[0] for s in stats], color=[s[4] for s in stats], s=[s[5] for s in stats])
	x,y 	= [ s[2] for s in stats if s[4]=="red"],[ s[0] for s in stats if s[4]=="red"]
		
	if len(x):
		specs 	= list(set([s[6] for s in stats if s[4]=="red" ]))
		w 		= 2
		ax.scatter(x,y,label="\n".join([ ",".join(specs[i:i+w]) for i in np.arange(0,len(specs),w)])+"\n"  ,color="red",s=70)


	x,y 	= [ s[2] for s in stats if s[4]=="green"],[ s[0] for s in stats if s[4]=="green"]
	if len(x):
		specs 	= list(set([s[6] for s in stats if s[4]=="green" ]))
		w 		= 2
		ax.scatter(x,y,label="\n".join([ ",".join(specs[i:i+w]) for i in np.arange(0,len(specs),w)])+"\n"  ,color="green",s=70)

	ax.set_title(title)
	ax.legend(loc="upper right",scatterpoints=3)
	ax.set_xscale("log")
	ax.set_ylim((-max(ax.get_ylim())-0.03,max(ax.get_ylim())+0.03))
	ax.set_ylim(-0.25,0.25)
	ax.set_xlim(0,pow(10,7))
	ax.set_xticks([pow(10,1),pow(10,2),pow(10,3),pow(10,4)])
	T 		= GRID(math.log(ax.get_xlim()[0],10),math.log(ax.get_xlim()[1],10),ax.get_ylim()[0], ax.get_ylim()[1] ,15)  
	T.insert_data_points(zip([ math.log(s[2],10) for s in stats], [s[0] for s in stats] ))
	PS 		= {}
	
	plt.show()
	return [s[0] for s in stats],[s[3] for s in stats]

if __name__ == "__main__":
	F1 	= "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_human/SRR1105737_enrichment_stats.tsv"
	F2 	= "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_human/SRR1105739_enrichment_stats.tsv"

	F1 	= load_enrichment_file(F1)
	F2 	= load_enrichment_file(F2)

	compare(F1,F2, "Differential MD Score")

