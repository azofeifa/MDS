
import load
import observe_pvalue_dist
import os
if __name__ == "__main__":
	D 			= "/Users/joazofeifa/Lab/new_motif_distances/out_files/"
	for DMSO in os.listdir(D):
		print DMSO
		DMSO 	= load.load_stats_file(D+DMSO,BOOT=True)
		observe_pvalue_dist.hist(DMSO, window=0.25)
		#observe_pvalue_dist.bootstrapped(DMSO, IND=True)

