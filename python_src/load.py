class PSSM:
	def __init__(self, name, MD, E,N):
		self.name 		= name
		self.MD_score 	= MD
		self.E_score 	= E
		self.N 			= N
		self.pv 		= {"MD":[1.0,1.0], "E":[1.0,1.0]}
		self.bt 		= {"MD":list(), "E":list()}
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
				name,N, MD, E 	= line_arrray[:4]
				MD,E,N 			= float(MD), float(E),float(N)
				P 				= PSSM(name, MD, E,N)
				P.pv["MD"] 		= [float(x) for x in line_arrray[4].split(",")]
				P.pv["E"] 		= [float(x) for x in line_arrray[5].split(",")]
				PSSMS.append(P)
			elif boot and line[0]!="#":
				line_arrray 		= line.strip("\n").split("\t")
				PSSMS[i].bt['MD'] 	= [float(x) for x in line_arrray[1].split(",")]
				PSSMS[i].bt['E'] 	= [float(x) for x in line_arrray[2].split(",")]
				i+=1
	return PSSMS










if __name__ == "__main__":
	DMSO 	= "/Users/joazofeifa/Lab/new_motif_distances/out_files/Allen2014_DMSO2_3_enrichment_stats.txt"
	Nutlin 	= "/Users/joazofeifa/Lab/new_motif_distances/out_files/Allen2014_Nutlin_3_enrichment_stats.txt"

	DMSO 	= load_stats_file(DMSO)
	#Nutlin 	= load_stats_file(Nutlin)