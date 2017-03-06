# Extrai aminoacidos a partir de centroide
# Regra: sao extraidos aminoacidos cuja distancia para o caborno alfa eh inferir a 6, 9, 12 ou 15 angstrons

from Bio.PDB import *
from math import sqrt
import sys
import warnings
warnings.filterwarnings("ignore") # esconde erros do pdb

parser = PDBParser()

raio = 15

# Resultado
w = open("aminoacidos_"+str(raio)+"A.txt",'w') #formato fasta
w2 = open("posicao_aminoacidos_"+str(raio)+"A.txt","w") #formato proprio arquivo (tab) E166(aminoacido, numero) (tab) E167... etc

# Le centroide
lista_centroide = open("centroide.txt").readlines()

for linha in lista_centroide:

	celula = linha.split("\t")

	arquivo = celula[0]
	xC = float(celula[11])
	yC = float(celula[12])
	zC = float(celula[13])

	info_amino = celula[1].split(".")
	cadeia = info_amino[0].strip()
	

	# Navegando em cada PDB
	estrutura = parser.get_structure('PDB',arquivo)

	print "Lendo arquivo "+arquivo
	w.write(">"+arquivo+"\n")
	w2.write(arquivo)
	w2.write("\t"+cadeia)

	try:
		for residuo in estrutura[0][cadeia]:
			if is_aa(residuo):
				try:
					rid = residuo.id[1]
					rname = Polypeptide.three_to_one(residuo.resname)

					# Modificacao para pegar maior distancia
					coord = residuo['CA'].get_coord()
					xR = coord[0]
					yR = coord[1]
					zR = coord[2]

					# Calculando distancia
					dist = sqrt( (xC-xR)**2 + (yC-yR)**2 + (zC-zR)**2 )

					if dist < raio:
						w.write(rname)
						w2.write("\t"+str(rname)+str(rid))
				except:
					print "Falha ao detectar um carbono alfa."
	except:
		# Em caso de soh ter uma cadeia, e essa nao ter nome
		for cadeia in estrutura[0]:
			for residuo in cadeia:
				if is_aa(residuo):
					try:
						rid = residuo.id[1]
						rname = Polypeptide.three_to_one(residuo.resname)
						coord = residuo['CA'].get_coord()
						xR = coord[0]
						yR = coord[1]
						zR = coord[2]

						# Calculando distancia
						dist = sqrt( (xC-xR)**2 + (yC-yR)**2 + (zC-zR)**2 )

						if dist < raio:
							w.write(rname)
							w2.write("\t"+str(rname)+str(rid))
					except:
						print "Falha ao detectar um carbono alfa."

	w.write("\n")
	w2.write("\n")

w.close()
w2.close()