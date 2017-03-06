#Converte mutation
from Bio.PDB import *

data = open("data3.txt").readlines()

ref = ["H.119", "W.120", "N.163", "E.164", "W.166", "C.167", "L.171", "H.178", "N.220", "L.221", "V.222", "F.240", "N.297", "Y.298", "Y.299", "S.300", "A.302", "E.357", "W.404", "E.411", "W.412", "F.420"]

for line in data:

	pos = line[3:].strip()
	amino = line[:3]

	amino1 = Polypeptide.three_to_one(amino)

	for r in ref:
		r = r.split(".")
		pos2 = r[1].strip()

		if pos == pos2:
			print "\""+r[0]+pos2+amino1+"\","