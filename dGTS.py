# dGTS.py
# Author: Diego Mariano
# Date: August 11, 2016


from math import sqrt
import sys
import glob
import numpy as np

folders = glob.glob("data/*.csv")
print "* CALCULATING delta GTS *\n"
major = 0

for folder in folders:

	try:
		r = open(sys.argv[1]).readlines()
	except:
		r = open(folder).readlines()


	# SVD ************************************************************************************************
	print "Running SVD"

	matrix = []
	for row in range(len(r)):
		each_row = r[row].split(",")
		for col in range(len(each_row)):
			each_row[col] = int(each_row[col])
		matrix.append(each_row) 

	# transpose matrix
	matrix = map(list, zip(*matrix))
	
	U, s, V = np.linalg.svd(matrix)
	S = np.diag(s)
	V = V.transpose()

	print "Matrix S"
	print s

	factor = 20

	Sfactor = S[:factor,:factor]
	#Vfactor = V[:factor,:]
	#AUXfactor = np.dot(Sfactor, Vfactor)
	Vfactor = V[:,:factor] 
	AUXfactor = np.dot(Sfactor, Vfactor.transpose()) 
	
	R = AUXfactor.transpose()
	
	# Writing R in a csv file
	w = open(folder+"_svd.txt",'w')
	for row in R:
		for col in range(len(row)):
			if col != len(row)-1:
				w.write(str(row[col])+",")
			else:
				w.write(str(row[col]))
		w.write("\n")

	w.close()
	
	# r receives a new value
	r = open(folder+"_svd.txt").readlines()


	# FIM AREA SVD ***************************************************************************************


	num_results_acsm = len(r)
	tolerants = range(23)
	len_tolerants = len(tolerants)

	# number of cols
	len_acsm_line = len(r[0].split(","))

	#print num_results_acsm

	# Calculating (delta)GTS wild vs. mutant -#######################################################
	
	dgts = [] 
	dgts_m2 = []
	nearest_wild = []
	nearest_mutant = []

	#nearest ******************************************************************************************************************************
	nearest = []
	dist_all_nearest = []
	dist_nearest = 0


	for i in range(num_results_acsm):
		nearest.append(-1)
		dist_all_nearest.append(0)

	# para cada resultado acsm (todos)
	for i in range(num_results_acsm):

		nearest_min = 0
		dist_nearest_min = 1000000000000000
		c1_acsm_line = r[i].split(",")
		

		# compara com as 23 tolerantes
		for j in range(len_tolerants):
			#regra que impede que a bgl de cao seja comparada com ela mesmo (j == 11)
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


	print "How is the nearest?"
	print nearest

	print "\nDistance nearest: "
	print dist_all_nearest


	# All distances for wild -----------------------------------------------------

	all_dist_wild = []

	for i in range(len_tolerants):

		c1_acsm_line = r[num_results_acsm-1].split(",")

		# calcula distancia euclidiana
		dist_nearest = 0
		for k in range(len_acsm_line):

			c2_acsm_line = r[i].split(",")
			# distancia da tolerante j para a bgl mutante i
			dist_nearest += (float(c1_acsm_line[k]) - float(c2_acsm_line[k]))**2		
			

		dist_nearest = sqrt(dist_nearest)
		all_dist_wild.append(dist_nearest)

	print "\nAll distances of wild to tolerants"
	print all_dist_wild

	# END -- All distances for wild -----------------------------------------------------


	# fim nearest ******************************************************************************************************************************

	# Calculating dGTS ********************************

	print "\ndGTS: "

	wild = num_results_acsm - 1
	
	for i in range(num_results_acsm):

		mutant = i
		print "Mutant: "+str(mutant)+" (line "+str(mutant+1)+")"

		# METHOD 1: nearest tolerant
		print "Nearest mutant: "+str(nearest[mutant])+"\tNearest wild: "+str(nearest[wild])
		print "M: "+str(dist_all_nearest[mutant])+" - W:"+str(dist_all_nearest[wild])
		print dist_all_nearest[mutant] - dist_all_nearest[wild]
		nearest_wild.append(str(nearest[wild]))
		nearest_mutant.append(str(nearest[mutant]))
		dgts.append( dist_all_nearest[mutant] - dist_all_nearest[wild] )

		# METHOD 2: nearest tolerant for mutant (only one tolerant) 
		#print "Nearest mutant: "+str(nearest[mutant])		
		#print "M: "+str(dist_all_nearest[mutant])+" - W:"+str(all_dist_wild[nearest[mutant]])		
		#print dist_all_nearest[mutant] - all_dist_wild[nearest[mutant]]		
		dgts_m2.append( dist_all_nearest[mutant] - all_dist_wild[nearest[mutant]] )

		print ""

	# FIM Classificando com nearest ********************************

def translateMutation(m):
	
	#tradutions = ["L171A","F240E","L221F","N220A","L171G","Y298F","L221A","L221G","S300F","V222A","V222G","N220S","F240A","W166H","V222S","W412A","C167H","F240S","W166R","Y298H","N297S","F240R","A302H","A302S","A302R","L171I","E411S","W166N","L221I","C167T","L221N","F240I","N220T","V222N","Y298I","L221T","F240N","W166L","V222T","S300N","H178L","S300T","A302N","A302L","Y298W","L171D","F240K","W120Y","N220D","H178M","W166Y","V222D","F240M","F240Y","L221C","S300M","A302Y","W166Q","A302M","F420Y","L171Q","W120F","C167V","A302Q","W166F","L171V","V222E","H178F","S300V", "wild"]
	#tradutions = ["H119A","Y299E","L171F","W120A","S300E","H178F","N163A","A302E","N220F","E164A","W404E","L221F","W166A","W412E","V222F","C167A","F420E","N297F","L171A","H119G","Y298F","H178A","W120G","Y299F","N220A","N163G","S300F","L221A","E164G","A302F","V222A","W166G","E357F","F240A","C167G","W404F","N297A","L171G","E411F","Y298A","H178G","W412F","Y299A","N220G","H119P","S300A","L221G","W120P","E357A","V222G","N163P","W404A","F240G","E164P","E411A","N297G","W166P","W412A","Y298G","C167P","F420A","Y299G","L171P","H119R","S300G","H178P","W120R","A302G","N220P","N163R","E357G","L221P","E164R","W404G","V222P","W166R","E411G","F240P","C167R","W412G","N297P","L171R","F420G","Y298P","H178R","W120H","Y299P","N220R","N163H","S300P","L221R","E164H","A302P","V222R","W166H","E357P","F240R","C167H","W404P","N297R","L171H","E411P","Y298R","N220H","W412P","Y299R","L221H","F420P","S300R","V222H","H119S","A302R","F240H","W120S","E357R","N297H","N163S","W404R","Y298H","E164S","E411R","Y299H","W166S","W412R","S300H","C167S","F420R","A302H","L171S","H119N","E357H","H178S","W120N","W404H","N220S","E164N","E411H","L221S","W166N","W412H","V222S","C167N","F420H","F240S","L171N","H119I","N297S","H178N","W120I","Y298S","L221N","N163I","Y299S","V222N","E164I","A302S","F240N","W166I","E357S","Y298N","C167I","W404S","Y299N","L171I","E411S","S300N","H178I","W412S","A302N","N220I","F420S","E357N","L221I","H119T","W404N","V222I","W120T","E411N","F240I","N163T","W412N","N297I","E164T","F420N","Y298I","W166T","H119D","Y299I","C167T","W120D","S300I","L171T","N163D","A302I","H178T","E164D","E357I","N220T","W166D","W404I","L221T","C167D","E411I","V222T","L171D","W412I","F240T","H178D","F420I","N297T","N220D","H119L","Y298T","L221D","W120L","Y299T","V222D","N163L","S300T","F240D","E164L","A302T","N297D","W166L","E357T","Y298D","C167L","W404T","Y299D","H178L","E411T","S300D","N220L","W412T","A302D","V222L","F420T","E357D","F240L","H119W","W404D","N297L","N163W","E411D","Y298L","E164W","W412D","Y299L","C167W","F420D","S300L","L171W","H119C","A302L","H178W","W120C","E357L","N220W","N163C","W404L","L221W","E164C","E411L","V222W","W166C","W412L","F240W","L171C","F420L","N297W","H178C","H119K","Y298W","N220C","W120K","Y299W","L221C","N163K","S300W","V222C","E164K","A302W","F240C","W166K","E357W","N297C","C167K","E411W","Y298C","L171K","F420W","Y299C","H178K","H119Y","S300C","N220K","W120Y","A302C","L221K","N163Y","E357C","V222K","E164Y","W404C","F240K","W166Y","E411C","N297K","C167Y","W412C","Y298K","L171Y","F420C","Y299K","H178Y","H119Q","S300K","N220Y","W120Q","A302K","L221Y","N163Q","E357K","V222Y","E164Q","W404K","F240Y","W166Q","E411K","N297Y","C167Q","W412K","S300Y","L171Q","F420K","A302Y","H178Q","H119M","E357Y","N220Q","W120M","W404Y","L221Q","N163M","E411Y","V222Q","E164M","W412Y","F240Q","W166M","F420Y","N297Q","C167M","H119V","Y298Q","L171M","W120V","Y299Q","H178M","N163V","S300Q","N220M","E164V","A302Q","L221M","W166V","E357Q","V222M","C167V","W404Q","F240M","L171V","E411Q","N297M","H178V","W412Q","Y298M","N220V","F420Q","Y299M","L221V","H119E","S300M","F240V","W120E","A302M","N297V","N163E","E357M","Y298V","W166E","W404M","Y299V","C167E","E411M","S300V","L171E","W412M","A302V","H178E","F420M","E357V","N220E","H119F","W404V","L221E","W120F","E411V","V222E","N163F","W412V","F240E","E164F","F420V","N297E","W166F","Y298E","C167F", "wild"]
	#tradutions = ["tolerants/ambgl18_sitioALA240ASN302.pdb","tolerants/ambgl18_sitioTHR222ALA240THR300.pdb","tolerants/ambgl18_sitioALA240.pdb","tolerants/ambgl18_sitioTHR222ASN302.pdb","tolerants/ambgl18_sitioALA240THR300ASN302.pdb","tolerants/ambgl18_sitioTHR222.pdb","tolerants/ambgl18_sitioALA240THR300.pdb","tolerants/ambgl18_sitioTHR222THR300ASN302.pdb","tolerants/ambgl18_sitioASN302.pdb","tolerants/ambgl18_sitioTHR222THR300.pdb","tolerants/ambgl18_sitioTHR222ALA240ASN302.pdb","tolerants/ambgl18_sitioTHR300ASN302.pdb","tolerants/ambgl18_sitioTHR222ALA240.pdb","tolerants/ambgl18_sitioTHR300.pdb","tolerants/ambgl18_sitioTHR222ALA240THR300ASN302.pdb","thermo_tolerants/ambgl18_sitioALA240.pdb","thermo_tolerants/ambgl18_sitioSER222.pdb","thermo_tolerants/ambgl18_sitioALA240SER302.pdb","thermo_tolerants/ambgl18_sitioSER222SER302.pdb","thermo_tolerants/ambgl18_sitioALA240THR300.pdb","thermo_tolerants/ambgl18_sitioSER222THR300.pdb","thermo_tolerants/ambgl18_sitioALA240THR300SER302.pdb","thermo_tolerants/ambgl18_sitioSER222THR300SER302.pdb","thermo_tolerants/ambgl18_sitioSER222ALA240.pdb","thermo_tolerants/ambgl18_sitioSER302.pdb","thermo_tolerants/ambgl18_sitioSER222ALA240SER302.pdb","thermo_tolerants/ambgl18_sitioTHR300.pdb","thermo_tolerants/ambgl18_sitioSER222ALA240THR300.pdb","thermo_tolerants/ambgl18_sitioTHR300SER302.pdb","thermo_tolerants/ambgl18_sitioSER222ALA240THR300SER302.pdb","wild"]
	tradutions = ["yang/d0ALA125.pdb_temp.pdb", "yang/d0GLU298.pdb_temp.pdb", "yang/d0PHE226.pdb_temp.pdb", "yang/d0ALA126.pdb_temp.pdb", "yang/d0GLU299.pdb_temp.pdb", "yang/d0PHE227.pdb_temp.pdb", "yang/d0ALA169.pdb_temp.pdb", "yang/d0GLU301.pdb_temp.pdb", "yang/d0PHE228.pdb_temp.pdb", "yang/d0ALA170.pdb_temp.pdb", "yang/d0GLU399.pdb_temp.pdb", "yang/d0PHE246.pdb_temp.pdb", "yang/d0ALA172.pdb_temp.pdb", "yang/d0GLU407.pdb_temp.pdb", "yang/d0PHE296.pdb_temp.pdb", "yang/d0ALA173.pdb_temp.pdb", "yang/d0GLU415.pdb_temp.pdb", "yang/d0PHE297.pdb_temp.pdb", "yang/d0ALA177.pdb_temp.pdb", "yang/d0GLY125.pdb_temp.pdb", "yang/d0PHE298.pdb_temp.pdb", "yang/d0ALA184.pdb_temp.pdb", "yang/d0GLY126.pdb_temp.pdb", "yang/d0PHE299.pdb_temp.pdb", "yang/d0ALA226.pdb_temp.pdb", "yang/d0GLY169.pdb_temp.pdb", "yang/d0PHE301.pdb_temp.pdb", "yang/d0ALA227.pdb_temp.pdb", "yang/d0GLY170.pdb_temp.pdb", "yang/d0PHE353.pdb_temp.pdb", "yang/d0ALA228.pdb_temp.pdb", "yang/d0GLY172.pdb_temp.pdb", "yang/d0PHE399.pdb_temp.pdb", "yang/d0ALA246.pdb_temp.pdb", "yang/d0GLY173.pdb_temp.pdb", "yang/d0PHE406.pdb_temp.pdb", "yang/d0ALA296.pdb_temp.pdb", "yang/d0GLY177.pdb_temp.pdb", "yang/d0PHE407.pdb_temp.pdb", "yang/d0ALA297.pdb_temp.pdb", "yang/d0GLY184.pdb_temp.pdb", "yang/d0PRO125.pdb_temp.pdb", "yang/d0ALA298.pdb_temp.pdb", "yang/d0GLY226.pdb_temp.pdb", "yang/d0PRO126.pdb_temp.pdb", "yang/d0ALA299.pdb_temp.pdb", "yang/d0GLY227.pdb_temp.pdb", "yang/d0PRO169.pdb_temp.pdb", "yang/d0ALA301.pdb_temp.pdb", "yang/d0GLY228.pdb_temp.pdb", "yang/d0PRO170.pdb_temp.pdb", "yang/d0ALA353.pdb_temp.pdb", "yang/d0GLY296.pdb_temp.pdb", "yang/d0PRO172.pdb_temp.pdb", "yang/d0ALA399.pdb_temp.pdb", "yang/d0GLY297.pdb_temp.pdb", "yang/d0PRO173.pdb_temp.pdb", "yang/d0ALA406.pdb_temp.pdb", "yang/d0GLY298.pdb_temp.pdb", "yang/d0PRO177.pdb_temp.pdb", "yang/d0ALA407.pdb_temp.pdb", "yang/d0GLY299.pdb_temp.pdb", "yang/d0PRO184.pdb_temp.pdb", "yang/d0ALA415.pdb_temp.pdb", "yang/d0GLY301.pdb_temp.pdb", "yang/d0PRO226.pdb_temp.pdb", "yang/d0ARG125.pdb_temp.pdb", "yang/d0GLY353.pdb_temp.pdb", "yang/d0PRO227.pdb_temp.pdb", "yang/d0ARG126.pdb_temp.pdb", "yang/d0GLY399.pdb_temp.pdb", "yang/d0PRO228.pdb_temp.pdb", "yang/d0ARG169.pdb_temp.pdb", "yang/d0GLY406.pdb_temp.pdb", "yang/d0PRO246.pdb_temp.pdb", "yang/d0ARG170.pdb_temp.pdb", "yang/d0GLY407.pdb_temp.pdb", "yang/d0PRO296.pdb_temp.pdb", "yang/d0ARG172.pdb_temp.pdb", "yang/d0GLY415.pdb_temp.pdb", "yang/d0PRO297.pdb_temp.pdb", "yang/d0ARG173.pdb_temp.pdb", "yang/d0HIS126.pdb_temp.pdb", "yang/d0PRO298.pdb_temp.pdb", "yang/d0ARG177.pdb_temp.pdb", "yang/d0HIS169.pdb_temp.pdb", "yang/d0PRO299.pdb_temp.pdb", "yang/d0ARG184.pdb_temp.pdb", "yang/d0HIS170.pdb_temp.pdb", "yang/d0PRO301.pdb_temp.pdb", "yang/d0ARG226.pdb_temp.pdb", "yang/d0HIS172.pdb_temp.pdb", "yang/d0PRO353.pdb_temp.pdb", "yang/d0ARG227.pdb_temp.pdb", "yang/d0HIS173.pdb_temp.pdb", "yang/d0PRO399.pdb_temp.pdb", "yang/d0ARG228.pdb_temp.pdb", "yang/d0HIS177.pdb_temp.pdb", "yang/d0PRO406.pdb_temp.pdb", "yang/d0ARG246.pdb_temp.pdb", "yang/d0HIS226.pdb_temp.pdb", "yang/d0PRO407.pdb_temp.pdb", "yang/d0ARG296.pdb_temp.pdb", "yang/d0HIS227.pdb_temp.pdb", "yang/d0PRO415.pdb_temp.pdb", "yang/d0ARG297.pdb_temp.pdb", "yang/d0HIS246.pdb_temp.pdb", "yang/d0SER125.pdb_temp.pdb", "yang/d0ARG298.pdb_temp.pdb", "yang/d0HIS296.pdb_temp.pdb", "yang/d0SER126.pdb_temp.pdb", "yang/d0ARG299.pdb_temp.pdb", "yang/d0HIS297.pdb_temp.pdb", "yang/d0SER169.pdb_temp.pdb", "yang/d0ARG301.pdb_temp.pdb", "yang/d0HIS298.pdb_temp.pdb", "yang/d0SER170.pdb_temp.pdb", "yang/d0ARG353.pdb_temp.pdb", "yang/d0HIS299.pdb_temp.pdb", "yang/d0SER172.pdb_temp.pdb", "yang/d0ARG399.pdb_temp.pdb", "yang/d0HIS301.pdb_temp.pdb", "yang/d0SER173.pdb_temp.pdb", "yang/d0ARG406.pdb_temp.pdb", "yang/d0HIS353.pdb_temp.pdb", "yang/d0SER177.pdb_temp.pdb", "yang/d0ARG407.pdb_temp.pdb", "yang/d0HIS399.pdb_temp.pdb", "yang/d0SER184.pdb_temp.pdb", "yang/d0ARG415.pdb_temp.pdb", "yang/d0HIS406.pdb_temp.pdb", "yang/d0SER226.pdb_temp.pdb", "yang/d0ASN125.pdb_temp.pdb", "yang/d0HIS407.pdb_temp.pdb", "yang/d0SER227.pdb_temp.pdb", "yang/d0ASN126.pdb_temp.pdb", "yang/d0HIS415.pdb_temp.pdb", "yang/d0SER228.pdb_temp.pdb", "yang/d0ASN170.pdb_temp.pdb", "yang/d0ILE125.pdb_temp.pdb", "yang/d0SER246.pdb_temp.pdb", "yang/d0ASN172.pdb_temp.pdb", "yang/d0ILE126.pdb_temp.pdb", "yang/d0SER296.pdb_temp.pdb", "yang/d0ASN173.pdb_temp.pdb", "yang/d0ILE169.pdb_temp.pdb", "yang/d0SER297.pdb_temp.pdb", "yang/d0ASN177.pdb_temp.pdb", "yang/d0ILE170.pdb_temp.pdb", "yang/d0SER298.pdb_temp.pdb", "yang/d0ASN184.pdb_temp.pdb", "yang/d0ILE172.pdb_temp.pdb", "yang/d0SER299.pdb_temp.pdb", "yang/d0ASN227.pdb_temp.pdb", "yang/d0ILE173.pdb_temp.pdb", "yang/d0SER301.pdb_temp.pdb", "yang/d0ASN228.pdb_temp.pdb", "yang/d0ILE177.pdb_temp.pdb", "yang/d0SER353.pdb_temp.pdb", "yang/d0ASN246.pdb_temp.pdb", "yang/d0ILE184.pdb_temp.pdb", "yang/d0SER399.pdb_temp.pdb", "yang/d0ASN297.pdb_temp.pdb", "yang/d0ILE226.pdb_temp.pdb", "yang/d0SER406.pdb_temp.pdb", "yang/d0ASN298.pdb_temp.pdb", "yang/d0ILE227.pdb_temp.pdb", "yang/d0SER407.pdb_temp.pdb", "yang/d0ASN299.pdb_temp.pdb", "yang/d0ILE228.pdb_temp.pdb", "yang/d0SER415.pdb_temp.pdb", "yang/d0ASN353.pdb_temp.pdb", "yang/d0ILE246.pdb_temp.pdb", "yang/d0THR125.pdb_temp.pdb", "yang/d0ASN399.pdb_temp.pdb", "yang/d0ILE296.pdb_temp.pdb", "yang/d0THR126.pdb_temp.pdb", "yang/d0ASN406.pdb_temp.pdb", "yang/d0ILE297.pdb_temp.pdb", "yang/d0THR169.pdb_temp.pdb", "yang/d0ASN407.pdb_temp.pdb", "yang/d0ILE298.pdb_temp.pdb", "yang/d0THR170.pdb_temp.pdb", "yang/d0ASN415.pdb_temp.pdb", "yang/d0ILE299.pdb_temp.pdb", "yang/d0THR172.pdb_temp.pdb", "yang/d0ASP125.pdb_temp.pdb", "yang/d0ILE301.pdb_temp.pdb", "yang/d0THR173.pdb_temp.pdb", "yang/d0ASP126.pdb_temp.pdb", "yang/d0ILE353.pdb_temp.pdb", "yang/d0THR177.pdb_temp.pdb", "yang/d0ASP169.pdb_temp.pdb", "yang/d0ILE399.pdb_temp.pdb", "yang/d0THR184.pdb_temp.pdb", "yang/d0ASP170.pdb_temp.pdb", "yang/d0ILE406.pdb_temp.pdb", "yang/d0THR226.pdb_temp.pdb", "yang/d0ASP172.pdb_temp.pdb", "yang/d0ILE407.pdb_temp.pdb", "yang/d0THR227.pdb_temp.pdb", "yang/d0ASP173.pdb_temp.pdb", "yang/d0ILE415.pdb_temp.pdb", "yang/d0THR228.pdb_temp.pdb", "yang/d0ASP177.pdb_temp.pdb", "yang/d0LEU125.pdb_temp.pdb", "yang/d0THR246.pdb_temp.pdb", "yang/d0ASP184.pdb_temp.pdb", "yang/d0LEU126.pdb_temp.pdb", "yang/d0THR296.pdb_temp.pdb", "yang/d0ASP226.pdb_temp.pdb", "yang/d0LEU169.pdb_temp.pdb", "yang/d0THR297.pdb_temp.pdb", "yang/d0ASP227.pdb_temp.pdb", "yang/d0LEU170.pdb_temp.pdb", "yang/d0THR298.pdb_temp.pdb", "yang/d0ASP228.pdb_temp.pdb", "yang/d0LEU172.pdb_temp.pdb", "yang/d0THR301.pdb_temp.pdb", "yang/d0ASP246.pdb_temp.pdb", "yang/d0LEU173.pdb_temp.pdb", "yang/d0THR353.pdb_temp.pdb", "yang/d0ASP296.pdb_temp.pdb", "yang/d0LEU184.pdb_temp.pdb", "yang/d0THR399.pdb_temp.pdb", "yang/d0ASP297.pdb_temp.pdb", "yang/d0LEU226.pdb_temp.pdb", "yang/d0THR406.pdb_temp.pdb", "yang/d0ASP298.pdb_temp.pdb", "yang/d0LEU227.pdb_temp.pdb", "yang/d0THR407.pdb_temp.pdb", "yang/d0ASP299.pdb_temp.pdb", "yang/d0LEU228.pdb_temp.pdb", "yang/d0THR415.pdb_temp.pdb", "yang/d0ASP301.pdb_temp.pdb", "yang/d0LEU246.pdb_temp.pdb", "yang/d0TRP125.pdb_temp.pdb", "yang/d0ASP353.pdb_temp.pdb", "yang/d0LEU296.pdb_temp.pdb", "yang/d0TRP169.pdb_temp.pdb", "yang/d0ASP399.pdb_temp.pdb", "yang/d0LEU297.pdb_temp.pdb", "yang/d0TRP170.pdb_temp.pdb", "yang/d0ASP406.pdb_temp.pdb", "yang/d0LEU298.pdb_temp.pdb", "yang/d0TRP172.pdb_temp.pdb", "yang/d0ASP407.pdb_temp.pdb", "yang/d0LEU299.pdb_temp.pdb", "yang/d0TRP173.pdb_temp.pdb", "yang/d0ASP415.pdb_temp.pdb", "yang/d0LEU301.pdb_temp.pdb", "yang/d0TRP177.pdb_temp.pdb", "yang/d0CYS125.pdb_temp.pdb", "yang/d0LEU353.pdb_temp.pdb", "yang/d0TRP184.pdb_temp.pdb", "yang/d0CYS126.pdb_temp.pdb", "yang/d0LEU399.pdb_temp.pdb", "yang/d0TRP226.pdb_temp.pdb", "yang/d0CYS169.pdb_temp.pdb", "yang/d0LEU406.pdb_temp.pdb", "yang/d0TRP227.pdb_temp.pdb", "yang/d0CYS170.pdb_temp.pdb", "yang/d0LEU407.pdb_temp.pdb", "yang/d0TRP228.pdb_temp.pdb", "yang/d0CYS172.pdb_temp.pdb", "yang/d0LEU415.pdb_temp.pdb", "yang/d0TRP246.pdb_temp.pdb", "yang/d0CYS177.pdb_temp.pdb", "yang/d0LYS125.pdb_temp.pdb", "yang/d0TRP296.pdb_temp.pdb", "yang/d0CYS184.pdb_temp.pdb", "yang/d0LYS126.pdb_temp.pdb", "yang/d0TRP297.pdb_temp.pdb", "yang/d0CYS226.pdb_temp.pdb", "yang/d0LYS169.pdb_temp.pdb", "yang/d0TRP298.pdb_temp.pdb", "yang/d0CYS227.pdb_temp.pdb", "yang/d0LYS170.pdb_temp.pdb", "yang/d0TRP299.pdb_temp.pdb", "yang/d0CYS228.pdb_temp.pdb", "yang/d0LYS172.pdb_temp.pdb", "yang/d0TRP301.pdb_temp.pdb", "yang/d0CYS246.pdb_temp.pdb", "yang/d0LYS173.pdb_temp.pdb", "yang/d0TRP353.pdb_temp.pdb", "yang/d0CYS296.pdb_temp.pdb", "yang/d0LYS177.pdb_temp.pdb", "yang/d0TRP406.pdb_temp.pdb", "yang/d0CYS297.pdb_temp.pdb", "yang/d0LYS184.pdb_temp.pdb", "yang/d0TRP415.pdb_temp.pdb", "yang/d0CYS298.pdb_temp.pdb", "yang/d0LYS226.pdb_temp.pdb", "yang/d0TYR125.pdb_temp.pdb", "yang/d0CYS299.pdb_temp.pdb", "yang/d0LYS227.pdb_temp.pdb", "yang/d0TYR126.pdb_temp.pdb", "yang/d0CYS301.pdb_temp.pdb", "yang/d0LYS228.pdb_temp.pdb", "yang/d0TYR169.pdb_temp.pdb", "yang/d0CYS353.pdb_temp.pdb", "yang/d0LYS246.pdb_temp.pdb", "yang/d0TYR170.pdb_temp.pdb", "yang/d0CYS399.pdb_temp.pdb", "yang/d0LYS296.pdb_temp.pdb", "yang/d0TYR172.pdb_temp.pdb", "yang/d0CYS406.pdb_temp.pdb", "yang/d0LYS297.pdb_temp.pdb", "yang/d0TYR173.pdb_temp.pdb", "yang/d0CYS407.pdb_temp.pdb", "yang/d0LYS298.pdb_temp.pdb", "yang/d0TYR177.pdb_temp.pdb", "yang/d0CYS415.pdb_temp.pdb", "yang/d0LYS299.pdb_temp.pdb", "yang/d0TYR184.pdb_temp.pdb", "yang/d0GLN125.pdb_temp.pdb", "yang/d0LYS301.pdb_temp.pdb", "yang/d0TYR226.pdb_temp.pdb", "yang/d0GLN126.pdb_temp.pdb", "yang/d0LYS353.pdb_temp.pdb", "yang/d0TYR227.pdb_temp.pdb", "yang/d0GLN169.pdb_temp.pdb", "yang/d0LYS399.pdb_temp.pdb", "yang/d0TYR228.pdb_temp.pdb", "yang/d0GLN170.pdb_temp.pdb", "yang/d0LYS406.pdb_temp.pdb", "yang/d0TYR246.pdb_temp.pdb", "yang/d0GLN172.pdb_temp.pdb", "yang/d0LYS407.pdb_temp.pdb", "yang/d0TYR296.pdb_temp.pdb", "yang/d0GLN173.pdb_temp.pdb", "yang/d0LYS415.pdb_temp.pdb", "yang/d0TYR299.pdb_temp.pdb", "yang/d0GLN177.pdb_temp.pdb", "yang/d0MET125.pdb_temp.pdb", "yang/d0TYR301.pdb_temp.pdb", "yang/d0GLN184.pdb_temp.pdb", "yang/d0MET126.pdb_temp.pdb", "yang/d0TYR353.pdb_temp.pdb", "yang/d0GLN226.pdb_temp.pdb", "yang/d0MET169.pdb_temp.pdb", "yang/d0TYR399.pdb_temp.pdb", "yang/d0GLN227.pdb_temp.pdb", "yang/d0MET170.pdb_temp.pdb", "yang/d0TYR406.pdb_temp.pdb", "yang/d0GLN228.pdb_temp.pdb", "yang/d0MET172.pdb_temp.pdb", "yang/d0TYR407.pdb_temp.pdb", "yang/d0GLN246.pdb_temp.pdb", "yang/d0MET173.pdb_temp.pdb", "yang/d0TYR415.pdb_temp.pdb", "yang/d0GLN296.pdb_temp.pdb", "yang/d0MET177.pdb_temp.pdb", "yang/d0VAL125.pdb_temp.pdb", "yang/d0GLN297.pdb_temp.pdb", "yang/d0MET184.pdb_temp.pdb", "yang/d0VAL126.pdb_temp.pdb", "yang/d0GLN298.pdb_temp.pdb", "yang/d0MET226.pdb_temp.pdb", "yang/d0VAL169.pdb_temp.pdb", "yang/d0GLN299.pdb_temp.pdb", "yang/d0MET227.pdb_temp.pdb", "yang/d0VAL170.pdb_temp.pdb", "yang/d0GLN301.pdb_temp.pdb", "yang/d0MET228.pdb_temp.pdb", "yang/d0VAL172.pdb_temp.pdb", "yang/d0GLN353.pdb_temp.pdb", "yang/d0MET246.pdb_temp.pdb", "yang/d0VAL173.pdb_temp.pdb", "yang/d0GLN399.pdb_temp.pdb", "yang/d0MET296.pdb_temp.pdb", "yang/d0VAL177.pdb_temp.pdb", "yang/d0GLN406.pdb_temp.pdb", "yang/d0MET297.pdb_temp.pdb", "yang/d0VAL184.pdb_temp.pdb", "yang/d0GLN407.pdb_temp.pdb", "yang/d0MET298.pdb_temp.pdb", "yang/d0VAL226.pdb_temp.pdb", "yang/d0GLN415.pdb_temp.pdb", "yang/d0MET299.pdb_temp.pdb", "yang/d0VAL228.pdb_temp.pdb", "yang/d0GLU125.pdb_temp.pdb", "yang/d0MET301.pdb_temp.pdb", "yang/d0VAL246.pdb_temp.pdb", "yang/d0GLU126.pdb_temp.pdb", "yang/d0MET353.pdb_temp.pdb", "yang/d0VAL296.pdb_temp.pdb", "yang/d0GLU169.pdb_temp.pdb", "yang/d0MET399.pdb_temp.pdb", "yang/d0VAL297.pdb_temp.pdb", "yang/d0GLU172.pdb_temp.pdb", "yang/d0MET406.pdb_temp.pdb", "yang/d0VAL298.pdb_temp.pdb", "yang/d0GLU173.pdb_temp.pdb", "yang/d0MET407.pdb_temp.pdb", "yang/d0VAL299.pdb_temp.pdb", "yang/d0GLU177.pdb_temp.pdb", "yang/d0MET415.pdb_temp.pdb", "yang/d0VAL301.pdb_temp.pdb", "yang/d0GLU184.pdb_temp.pdb", "yang/d0PHE125.pdb_temp.pdb", "yang/d0VAL353.pdb_temp.pdb", "yang/d0GLU226.pdb_temp.pdb", "yang/d0PHE126.pdb_temp.pdb", "yang/d0VAL399.pdb_temp.pdb", "yang/d0GLU227.pdb_temp.pdb", "yang/d0PHE169.pdb_temp.pdb", "yang/d0VAL406.pdb_temp.pdb", "yang/d0GLU228.pdb_temp.pdb", "yang/d0PHE170.pdb_temp.pdb", "yang/d0VAL407.pdb_temp.pdb", "yang/d0GLU246.pdb_temp.pdb", "yang/d0PHE173.pdb_temp.pdb", "yang/d0VAL415.pdb_temp.pdb", "yang/d0GLU296.pdb_temp.pdb", "yang/d0PHE177.pdb_temp.pdb", "yang/d0GLU297.pdb_temp.pdb", "yang/d0PHE184.pdb_temp.pdb", "wild"]
	return tradutions[m]

def translateSpecie(s):
	s = int(s)
	specie = ["Acidilobus saccharovorans", "Bacillus subtilis", "Caldicellulosiruptor bescii", "Exiguobacterium antarcticum B7", "Fervidobacterium islandicum", "Humicola grisea var thermoidea", "Hypocrea jecorina Trichoderma reesei", "Metagenome China South Sea", "Metagenome hydrothermal spring", "Metagenome Kusaya gravy", "Metagenome soil", "Metagenome Turpan Depression", "Mucor circinelloides", "Nasutitermes takasagoensis", "Neotermes koshunensis", "Neurospora crassa", "Pyrococcus furiosus", "Talaromyces funiculosus Penicillium funiculosum", "Thermoanaerobacter brockii", "Thermoanaerobacterium aotearoense", "Thermoanaerobacterium thermosaccharolyticum", "Thermotoga naphthophila", "Thermotoga petrophila"]
	return specie[s]

# Saving results
w = open("result.txt",'w')
w.write("Mutation\tdGTS Method 1\tdGTS Method 2\tWild reference\tMutate reference\n")

for i in range(len(dgts)):

	if i >= len_tolerants:
		mutation_name = translateMutation(i-len_tolerants)
		mutant_name = translateSpecie(nearest_mutant[i])
		wild_name = translateSpecie(nearest_wild[i])

		w.write(str(mutation_name)+"\t"+str(dgts[i])+"\t"+str(dgts_m2[i])+"\t"+str(wild_name)+"\t"+str(mutant_name)+"\n")

w.close()

print "Successful | Results saved at result.txt"