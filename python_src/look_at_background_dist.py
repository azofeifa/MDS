import matplotlib.pyplot as plt
def load_for_freq(FILE):
	X 	= list()
	collect=False 
	with open(FILE) as FH:
		for line in FH:
			if collect:
				x 	= [float("0." + i) for i in line.strip("\n").split("\t")[1].split(".")[1:]]
				X.append(x)
			elif "#Est" in line:
				collect 	= True
	for i in range(4):

		plt.plot(range(len(X)), [x[i] for x in X])
	plt.ylim(0.2,0.3)
	plt.show()
def look_at_BS(FILE):
	with open(FILE) as FH:
		for line in FH:
			if line[0]== ">":
				M 	= line[1:].strip("\n")
			elif line[0]=="~":
				x 	= [int(i) for i in line.split("|")[1].split(",")]
				plt.title(M)
				plt.bar(range(len(x)), x)
				plt.show()



if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/new_motif_distances/PSSM_DBS/new_Allen2014_generated_DB_3.txt"
	look_at_BS(FILE)
