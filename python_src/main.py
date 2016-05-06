
import load
import observe_pvalue_dist

if __name__ == "__main__":
	DMSO 	= "/Users/joazofeifa/Lab/new_motif_distances/out_files/Allen2014_DMSO2_3_enrichment_stats.txt"
	Nutlin 	= "/Users/joazofeifa/Lab/new_motif_distances/out_files/Allen2014_Nutlin_3_enrichment_stats.txt"

	DMSO 	= load.load_stats_file(Nutlin,BOOT=True)
	observe_pvalue_dist.hist(DMSO)
	observe_pvalue_dist.bootstrapped(DMSO, IND=True)