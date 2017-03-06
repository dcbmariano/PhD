import glob
import sys
from Bio.PDB import *
from math import sqrt


raio = 6 #default

for i in range(len(sys.argv)):
	
	if i % 2  == 1:

		# Ligante
		if sys.argv[i] == '-r':
			pasta_receptores = sys.argv[i+1]
			receptores = glob.glob(pasta_receptores+"/*")

		# Diretorio receptores
		if sys.argv[i] == '-l':
			pasta_ligantes = sys.argv[i+1]
			ligantes = glob.glob(pasta_ligantes+"/*")

		if sys.argv[i] == '-d':
			raio = int(sys.argv[i+1])

		# Help
		if sys.argv[i] == '-h' or sys.argv[i] == '--help':
			print "Run: python retorna_aminoacidos_proximos.py -l pasta_ligantes -d raio -r pasta_receptores\nNao use '/' ao final do nome da pasta"



for ligante in ligantes:
	aux = ligante.split(".")
	formato = aux[len(aux)-1]

	if formato == "pdbqt" or formato == "pdb":
		print "Lendo "+ligante

		arquivo_ligante = open(ligante).readlines()
		aux = ligante.split("docked_") # defini como padrao que a string 'docked_' eh o que difere um nome de outro
		nome_ligante = aux[len(aux)-1] #retira o nome 

		for linha in arquivo_ligante:

			if linha[0:4] == 'ATOM':
				x = float(linha[30:38])
				y = float(linha[38:46])
				z = float(linha[46:54])


				# CALCULA A DISTANCIA DESSE ATOMO PARA TODOS OS OUTROS
				#arquivo_receptor = open(pasta_receptores+"/"+nome_ligante[0:-2]) #-2 remove o qt de pdbqt
				try:
					arquivo_receptor = open(pasta_receptores+"/"+nome_ligante) #-2 remove o qt de pdbqt

					residuos = []

					for linha2 in arquivo_receptor:
						if linha2[0:4] == 'ATOM':
							x2 = float(linha2[30:38])
							y2 = float(linha2[38:46])
							z2 = float(linha2[46:54])

							#calcula distancia euclidiana
							dist = sqrt( (x2-x)**2 + (y2-y)**2 + (z2-z)**2)

							if dist < raio:
								residuo = linha2[17:20]+linha2[22:26]
								if residuo not in residuos:
									residuos.append(residuo)


						elif linha[0:3] == 'TER':
							break # pega apenas a primeira cadeia
				except:
					print "."

			elif linha == "MODEL 2\n":
				break
		print residuos