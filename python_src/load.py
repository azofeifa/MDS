from scipy.stats import norm
class PSSM:
	def __init__(self, name, MD,N):
		self.name 		= name
		self.MD_score 	= MD
		self.N 			= N
		self.pv 		= None
		self.pv2 		= None
		self.bt 		= list()

def load_stats_file(FILE, BOOT=False):
	PSSMS 	= list()
	i 			= 0
	begin 	= True
	boot 	= False
	OBS 	= False
	with open(FILE) as FH:
		for line in FH:
			if begin and  "#Binned" == line[:7]:
				begin 	= False
				OBS 	= True
			elif OBS and "#Empiracle" in line:
				begin 	= False
				boot 	= BOOT
				OBS 	= False
	
			elif begin and line[0]!="#":

				line_arrray 	= line.strip("\n").split("\t")
				name,N, MD 	= line_arrray[:3]
				MD,N 			= float(MD), float(N)
				P 				= PSSM(name, MD, N)
				P.pv 			= [float(x) for x in line_arrray[3].split(",")]
				P.pv2 			= [float(x) for x in line_arrray[4].split(",")]
				P.null 			= [float(x) for x in line_arrray[5].split(",")]

				PSSMS.append(P)
			elif boot and line[0]!="#":
				line_arrray 		= line.strip("\n").split("\t")
				PSSMS[i].bt 		= [float(x) for x in line_arrray[1].split(",")]
				i+=1
	return PSSMS










if __name__ == "__main__":
	DMSO 	= "/Users/joazofeifa/Lab/new_motif_distances/out_files/Allen2014_DMSO_3_enrichment_stats_250.txt"

	DMSO 	= load_stats_file(DMSO)
	#Nutlin 	= load_stats_file(Nutlin)