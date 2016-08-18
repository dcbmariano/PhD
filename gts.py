# GTS-Score: glucose tolerance signature score
# Author: Diego Mariano
# Date: July 31, 2016

import sys
from math import sqrt
import numpy

print "Running GTS-Score"
separator = ","

# File 
try:
	tsv = sys.argv[1]
except:
	tsv = "gt6.5_c10.csv"

# Tolerant beta-glucosidases
try:
	tolerants_list = sys.argv[2]
except:
	#tolerants_list = "3,4,5,8,10,17,19,20,21,22,23" #6A
	tolerants_list = "3,4,5,8,12,17,19,20,21,22,23" #6.5A	
	#tolerants_list = "4,5,8,19,20,21,22,23" #6.5A TESTE SEM 17, 3, 12

	#tolerants_list = "4,8,9,19,20,21" #7A
	#tolerants_list = "4,8,9,11,19,20,21" #7.5A
	#tolerants_list = "3,4,8,9,10,11,19,20,21"

tolerants = tolerants_list.split(",")


# CONSTRUCTING THE DISTANCE MATRIX ************************************
print "Constructing the distance matrix"
dist = []
tsv_file = open(tsv).readlines()
len_tsv = len(tsv_file)
cont_tsv = 1

for row1 in tsv_file:

	print str(cont_tsv)+"/"+str(len_tsv)
	cont_tsv+=1

	pos1 = row1.split(separator)

	dist_row = []

	for row2 in tsv_file:

		pos2 = row2.split(separator)

		# the distance is given by sqrt of all (pos2 - pos1)**2

		possum = 0
		for k in range(len(pos1)):
			p1 = float(pos1[k])
			p2 = float(pos2[k])
			dif = (p2 - p1)**2
			possum += dif  

		dist_aux = sqrt(possum)

		dist_row.append(dist_aux)

	dist.append(dist_row)

#print dist

# CALCULATING CENTROID *************************************************
print "Calculating glucose-tolerant centroid"
i = 1
centroid = []

# Fill centroid list with zeros
for row in tsv_file:
	pos = row.split(separator)
	for cell in range(len(pos)):
		centroid.append(0)
	break


for row in tsv_file:

	# Verify if the actual line is in tolerants list
	if str(i) in tolerants:

		pos = row.split(separator)
		tam_pos = len(pos)

		# Analyze each cell (k) in pos1 and pos2
		for cell in range(tam_pos):
			centroid[cell] += float(pos[cell])

	i += 1

#average
for cell in range(len(centroid)):
	centroid[cell] = centroid[cell]/len(tolerants)

#print centroid

# CALCULATING GTS-SCORE ***************************************
print "Calculating GTS-SCORE"
w = open("gts-score.txt","w")

for row in tsv_file:
	# euclidian distance between a line and the centroid

	pos = row.split(separator)

	aux = 0
	dist_centroid = 0

	for cell in range(len(pos)):
		aux += (float(pos[cell]) - float(centroid[cell]))**2

	dist_centroid = sqrt(aux)
	w.write(str(dist_centroid)+"\n")

w.close()

print "The results were saved in the file gts-score.txt"


# DETERMINING THE NEAREST RESULTS *****************************

# First determining the most distant tolerant from the centroid
gts = open("gts-score.txt")
gts_lines = gts.readlines()

i = 1
max_gts = 0
max_global = 0
for line in gts_lines:

	# Determining the highest tolerant
	if str(i) in tolerants:
		if float(line) > max_gts:
			max_gts = float(line)

	# Determining the highest global
	if float(line) > max_global:
		max_global = float(line)

	i+=1

if max_gts < 0:
	max_gts = max_gts * -1

gts.close()

gts = open("gts-score.txt")
gts_lines = gts.readlines()

#max_gts is the ray
# Determining all with distance less than the ray
score_aux = max_gts/max_global
print "GTS-Ray: "+str(max_gts)
print "Max-Ray: "+str(max_global)
print "Score: "+str(score_aux)
cont = 0
i = 1
new_tolerants = ""

for line in gts_lines:

	if float(line) <= float(max_gts) and str(i) not in tolerants:
		new_tolerants = new_tolerants+str(i)+","
		cont += 1
	i += 1

print "We detected "+str(cont)+" enzymes possibly tolerants to glucose inhibition."
print "New tolerants: "+new_tolerants
gts.close()

# New: testing tolerants
'''
t1 = open("original.tsv").readlines()
t2 = open("mutada.tsv").readlines()
aux1 = 0
aux2 = 0

for i in range(len(centroid)):
	aux1 += (float(t1[i]) - float(centroid[i]))**2
	aux2 += (float(t2[i]) - float(centroid[i]))**2

t1_dist = sqrt(aux1)
t2_dist = sqrt(aux2)

gts_t1 = t1_dist / 1
gts_t2 = t2_dist / 1

print "GTS MUTADA: "+str(gts_t1)
print "GTS ORIGINAL: "+str(gts_t2)

dgts_parcial = t2_dist - t1_dist

dgts = dgts_parcial / 1
print "DGTS: "+str(dgts)



mut_20 = open("mutantes_20.tsv").readlines()
aux = 0
gts_mut_20 = []
delta_gts = []

for bgl in mut_20:
	fbgl = bgl.split("\t")

	for i in range(len(centroid)):
		aux += (float(fbgl[i]) - float(centroid[i]))**2
	
	aux = sqrt(aux)
	
	gts_atual = aux / max_global
	gts_mut_20.append(gts_atual)

print "GTS: "
for g in gts_mut_20:
	print g

#Calculating deltaGTS
for i in range(82):
	dgts_parcial = gts_mut_20[i] - gts_mut_20[i+82]
	delta_gts.append(dgts_parcial)

print "DeltaGTS: "

for g in delta_gts:
	print g

'''

print "Successful"