import matplotlib.pyplot as plt
import numpy as np
import math
import os
import load
from matplotlib import rc
def hist(PSSMS,window=0.25):
	F 	= plt.figure()
	#ax1 = F.add_subplot(1,2,1)
	ax2 = F.add_subplot(1,1,1)

	#ax1.hist([P.pv[1]  for P in PSSMS ],bins=35, edgecolor="white", color="green")
	colors 	= ["red" if p.pv[1] <0.001 and p.pv2[1]<0.001 else "blue" for p in PSSMS  ]
	EXS 	= [p.null[0]  for p in PSSMS]
	OBS 	= [p.MD_score  for p in PSSMS]
	VARS 	= [ math.log( 1.0 / p.null[1],10) for p in PSSMS]
	FCS 	= [p.MD_score - (p.null[0] )  for p in PSSMS]

	ax2.scatter([VARS[i] for i in range(len(colors)) if colors[i]=="blue"] , [FCS[i] for i in range(len(colors)) if colors[i]=="blue"] , color="blue", alpha=0.5)
	ax2.scatter([VARS[i] for i in range(len(colors)) if colors[i]=="red"] , [FCS[i] for i in range(len(colors)) if colors[i]=="red"] , color="red", alpha=1.0,label=r"$p_{s} < 0.001$" +"\n" r'$p_{ns}<0.001$')
	ax2.grid()
	ax2.set_ylim(-max(FCS)-0.1,max(FCS)+0.1)
	ax2.set_ylabel(r"$MD_{obs} - E[MD_{null}]$",fontsize=20)
	ax2.set_xlabel(r"$1/Var(MD_{null})$",fontsize=20)
	ax2.legend(loc="best")
	plt.show()
def bootstrapped(PSSMS,IND=False):
	if (IND):
		for P in PSSMS:
			F 	= plt.figure(figsize=(15,10))
			ax1 	= F.add_subplot(111)
			ax1.set_title(P.name+ ":" + str(P.MD_score))
			ax1.hist( P.bt , bins=30, edgecolor="white")

			plt.show()

	F 	= plt.figure()
	ax1 = F.add_subplot(111)
	ax1.hist([ np.mean(P.bt["MD"]) for P in PSSMS],bins=30)

	plt.show()

if __name__ == "__main__":
	D 			= "/Users/joazofeifa/Lab/new_motif_distances/out_files/"
	for DMSO in os.listdir(D):
		print DMSO
		DMSO 	= load.load_stats_file(D+DMSO,BOOT=True)
		hist(DMSO, window=0.25)
		#observe_pvalue_dist.bootstrapped(DMSO, IND=True)




