import matplotlib.pyplot as plt
def load_for_freq(FILE):
	X 	= list()
	collect=False 
	with open(FILE) as FH:
		for line in FH:
			if collect:
				x 	= [float(i) for i in line.strip("\n").split("\t")[1].split(",")]
				X.append(x)
			elif "#Est" in line:
				collect 	= True
	for i in range(4):

		plt.plot(range(len(X)), [x[i] for x in X])
	plt.ylim(0.2,0.35)
	plt.show()
def look_at_BS(FILE, G):
	with open(FILE) as FH:
		for line in FH:
			if line[0]== ">":
				M 	= line[1:].strip("\n")
				if M not in G:
					G[M] 	= list()
			elif line[0]=="~":
				x 	= [int(i) for i in line.split("|")[1].split(",")]
				G[M].append(x)
	return G

def show(G):
	for m in G:
		F 		= plt.figure()
		print m
		ax1 	= F.add_subplot(1,2,1)
		ax1.hist(range(2000),  weights=G[m][0], bins=50, edgecolor="white")

		ax2 	= F.add_subplot(1,2,2)
		ax2.hist(range(2000), weights=G[m][1], bins=50, edgecolor="white")
		#ax1.set_ylim(ax2.get_ylim())
		plt.show()

if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/new_motif_distances/PSSM_DBS/new_Allen2014_generated_DB_3_MM.txt"
	FILE2 	= "/Users/joazofeifa/Lab/new_motif_distances/PSSM_DBS/new_Allen2014_generated_DB_3.txt"
	#load_for_freq(FILE)
	G 		= {}
	G 		= look_at_BS(FILE,G)
	G 		= look_at_BS(FILE2, G)
	show(G)