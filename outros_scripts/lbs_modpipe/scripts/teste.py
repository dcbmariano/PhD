# -*- coding: utf-8 -*-

from modeller import *         	 
from modeller.automodel import *    

log.verbose()	# request verbose output
env = environ()  # create a new MODELLER environment to build this model in

env.io.atom_files_directory = ['/home/diego/bin/lbs_modpipe/scripts']

env.io.hetatm = True

env.io.water = True

a = automodel(env,
       	alnfile  = 'alinhamento.pir', 	# alignment filename
       	knowns   = '4ptv',          	# codes of the templates ## Molde ()
       	sequence = 'bgl')          	# code of the target ## Modelo ()

a.starting_model= 1             	# index of the first model
a.ending_model  = 100             	# index of the last model
a.make()                        	# do the actual homology modeling
