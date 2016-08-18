# calculateDist.py
# Author: Diego Mariano
# Date: August 11, 2016


from math import sqrt
import sys
import glob

folders = glob.glob("data/*.csv")
print "* CALCULATE DIST *\n"
major = 0

for folder in folders:

	try:
		r = open(sys.argv[1]).readlines()
	except:
		r = open(folder).readlines()

	num_results_acsm = len(r)

	tolerants = [19] 

	len_tolerants = len(tolerants)



	################################################################ Determine centroid ################################################################

	# Filling centroid with 0
	centroid = []

	for i in range(num_results_acsm):
		len_line_1 = len(r[i].split(","))

		for value in range(len_line_1):
			centroid.append(0)
		break

	len_centroid = len(centroid)

	# For each line of an aCSM result in position i
	# CENTROID (n elements) (x1,y1)(x2,y2) => ( (x1+x2)/n,(y2+y2)/n )

	for i in range(num_results_acsm):
		
		# If the line represent a tolerant 
		if i+1 in tolerants:

			acsm_line = r[i].split(",")
			len_acsm_line = len(acsm_line)

			# Sum each centroid
			for j in range(len_acsm_line):
				centroid[j] += float(acsm_line[j])


	#average
	for i in range(len_centroid):
		centroid[i] = centroid[i]/len_tolerants

	#print "CENTROID: "
	#print centroid


	# Calculating euclidian distance #################################################

	aux_tolerant = 0
	dist_aux = 0
	dist_centroid = []

	
	for i in range(num_results_acsm):		
		acsm_line = r[i].split(",")
		for j in range(len(acsm_line)):
			dist_aux += (float(acsm_line[j]) - float(centroid[j]))**2		
		
		dist_aux = sqrt(dist_aux)
		dist_centroid.append(dist_aux)




	# Calculating (delta)GTS wild vs. mutant -#######################################################

	#gts_expected = [-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,-1,1,-1,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,1,1,1,1,-1,-1,1,1,1,1,1,0,1,1,1,1,1]
	#gts_expected = [-1,0,0,0,0,1,1,0,-1,-1,-1,-1,-1,0,0,-1,-1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,-1,1,-1,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0]
	gts_expected = [-1,0,0,1,-1,0,0,-1,-1,0,-1,-1,-1,0,0,-1,-1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,-1,1,-1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,-1,1,1,0,0,0]

	dgts = [] 

	#nearest ******************************************************************************************************************************
	nearest = []
	dist_all_nearest = []
	dist_nearest = 0


	for i in range(num_results_acsm):
		nearest.append(-1)
		dist_all_nearest.append(0)

	# para cada resultado acsm (todos)
	for i in range(187):

		nearest_min = 0
		dist_nearest_min = 9999999999999999999999999
		c1_acsm_line = r[i].split(",")

		

		# compara com as 23 tolerantes
		for j in range(23):
			#regra que impede que a bgl de cao seja comparada com ela mesmo (j == 11)
			if j!=11:
				c2_acsm_line = r[j].split(",")

				dist_nearest = 0

				# calcula distancia euclidiana
				for k in range(len_acsm_line):

					# distancia da tolerante j para a bgl mutante i
					dist_nearest += (float(c1_acsm_line[k]) - float(c2_acsm_line[k]))**2		
			

				dist_nearest = sqrt(dist_nearest)

				if dist_nearest < dist_nearest_min:
					dist_nearest_min = dist_nearest
					nearest_min = j


		nearest[i] = nearest_min
		dist_all_nearest[i] = dist_nearest_min


	#print nearest
	#print dist_all_nearest

	
	# fim nearest ******************************************************************************************************************************

	# Classificando com nearest ********************************
	for i in range(24,106):

		mutant = i - 1
		wild = (i + 82) - 1
		'''
		print mutant+1
		print str(nearest[mutant])+"x"+str(nearest[wild])
		print "M: "+str(dist_all_nearest[mutant])+" - W:"+str(dist_all_nearest[wild])
		print dist_all_nearest[mutant] - dist_all_nearest[wild]
		print "-"
'''
		dgts.append( dist_all_nearest[mutant] - dist_all_nearest[wild] )

	# FIM Classificando com nearest ********************************

	# Classificando com centroide *********************************
	'''
	for i in range(24,106):

		mutant = i - 1
		wild = (i + 82) - 1

		#print str(wild)+"x"+str(mutant)
		dgts.append( dist_centroid[mutant] - dist_centroid[wild] )
	'''

	# FIM Classificando com centroide *********************************

	positive = 0
	negative = 0
	neutral = 0
	accuracy = 0
	print dgts

	for i in range(len(dgts)):
		if dgts[i] * gts_expected[i] > 0:
			positive += 1
		elif dgts[i] * gts_expected[i] < 0:
			negative += 1
		elif dgts[i] * gts_expected[i] == 0:
			neutral += 1

	accuracy = 100 * positive / (positive + negative)
	if accuracy > major:
		major = accuracy

	print folder+": "+str(accuracy)+"%"

print "Higher: "+str(major)+"%"