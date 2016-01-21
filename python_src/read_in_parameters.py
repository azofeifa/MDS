import os
def checkFile(f):
	vl, text 	= os.path.isfile(f), ""

	if not vl:
		text 	= " does not exist..."
	return vl, text

def checkInt(vl):
	text 	= ""
	try:
		b 	= int(vl)
	except:
		text = " is not an integer..."
		return False, text
	return True,text
def checkFloat(vl):
	text 	= ""
	try:
		b 	= float(vl)
	except:
		text = " is not a floating point..."
		return False, text
	return True,text


	
def checkDir(path):
	vl, text 	= os.path.isdir(path), ""
	if not vl:
		text 	= " is not a valid path..."
	return vl, text
def checkChrom(vl):
	chroms 	=["all"] + ["chr" + str(i) for i in range(1, 24)] + ["chrX", "chrY", "chrM"]
	if vl in chroms:
		return True, ""
	return False, "this not a valid chromosome identifier"

class argv_wrapper():
	def __init__(self):
		self.PSSM_DB 		= ""
		self.fasta_file 	= ""
		self.bed_file 		= ""
		self.background 	= [0.25,0.25,0.25,0.25]
		self.chromosome 	= "all"
		self.out_directory 	= ""
		self.job_ID 		= 0

		


def run(argv):
	ax 	= 

	return aw
