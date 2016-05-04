import matplotlib.pyplot as plt
import numpy as np
def hist(PSSMS):
	F 	= plt.figure()
	ax1 = F.add_subplot(1,2,1)
	ax2 = F.add_subplot(1,2,2)

	ax1.hist([P.pv["MD"][0]  for P in PSSMS],bins=30, edgecolor="white")
	ax2.hist([P.pv["MD"][1]  for P in PSSMS],bins=30, edgecolor="white")

	plt.show()
def bootstrapped(PSSMS,IND=False):
	if (IND):
		for P in PSSMS:
			F 	= plt.figure(figsize=(15,10))
			ax1 	= F.add_subplot(121)
			ax1.set_title(P.name+ ":" + str(P.MD_score))
			ax1.hist( P.bt["MD"] , bins=30, edgecolor="white")

			ax2 	= F.add_subplot(122)
			ax2.set_title(P.name+ ":" + str(P.E_score))
			ax2.hist( P.bt["E"] , bins=30, edgecolor="white")
			plt.show()

	F 	= plt.figure()
	ax1 = F.add_subplot(121)
	ax1.hist([ np.mean(P.bt["MD"]) for P in PSSMS],bins=30)
	ax2 = F.add_subplot(122)
	ax2.hist([ np.mean(P.bt["E"]) for P in PSSMS],bins=30)

	plt.show()



