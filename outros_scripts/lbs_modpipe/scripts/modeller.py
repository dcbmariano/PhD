# -*- coding: utf-8 -*-

from modeller import *         	 
from modeller.automodel import *    
import sys

log.verbose()	# request verbose output
env = environ()  # create a new MODELLER environment to build this model in


for i in len(sys.argv):
	if sys.argv[i] == "-d":
		diretorio = sys.argv[i+1]
	if sys.argv[i] == "-a":
		alinhamento = sys.argv[i+1]
	if sys.argv[i] == "-t":
		template = sys.argv[i+1]
	if sys.argv[i] == "-s":
		seq = sys.argv[i+1]
		

env.io.atom_files_directory = [diretorio]

env.io.hetatm = True

env.io.water = True

a = automodel(env,
       	alnfile  = alinhamento, 	# alignment filename
       	knowns   = template,          	# codes of the templates ## Molde ()
       	sequence = seq)          	# code of the target ## Modelo ()

a.starting_model= 1             	# index of the first model
a.ending_model  = 100             	# index of the last model
a.make()                        	# do the actual homology modeling
