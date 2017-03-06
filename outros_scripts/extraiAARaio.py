# Extrai aminoacidos a partir de centroide
# Regra: sao extraidos aminoacidos cuja distancia para o caborno alfa eh inferir a 6, 9, 12 ou 15 angstrons

from Bio.PDB import *
from math import sqrt
import sys
import warnings
warnings.filterwarnings("ignore") # esconde erros do pdb


raio = 9


parser = PDBParser()
matrix_atomos = []


# Resultado
w = open("aminoacidos_"+str(raio)+"A.txt",'w') #formato fasta
w2 = open("posicao_aminoacidos_"+str(raio)+"A.txt",'w') #formato proprio arquivo (tab) E166(aminoacido, numero) (tab) E167... etc


# Le centroide
lista_centroide = open("centroide.txt").readlines()

print "ETAPA 1: extrai atomos"

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

					# Modificacao para pegar o atomo mais proximo
					maior_dist = 0
					for atomo in residuo:
						coord = atomo.get_coord()
						xR = coord[0]
						yR = coord[1]
						zR = coord[2]

						# Calculando distancia
						dist = sqrt( (xC-xR)**2 + (yC-yR)**2 + (zC-zR)**2 )
						
						if dist > maior_dist:
							maior_dist = dist

					if maior_dist < raio:
						w.write(rname)
						w2.write("\t"+str(rname)+str(rid))
						# Se o raio eh menor, grave todos os atomos
						for atomo in residuo:
							coord = atomo.get_coord()
							x = coord[0]
							y = coord[1]
							z = coord[2]
							matrix_atomos.append([arquivo,cadeia,rname,x,y,z,rid])
				except:
					print "Falha ao detectar atomos."

	except:
		# Em caso de soh ter uma cadeia, e essa nao ter nome
		for cadeia in estrutura[0]:
			for residuo in cadeia:
				if is_aa(residuo):
					try:
						rid = residuo.id[1]
						rname = Polypeptide.three_to_one(residuo.resname)

						# Modificacao para pegar o atomo mais proximo
						maior_dist = 0
						for atomo in residuo:
							coord = atomo.get_coord()
							xR = coord[0]
							yR = coord[1]
							zR = coord[2]

							# Calculando distancia
							dist = sqrt( (xC-xR)**2 + (yC-yR)**2 + (zC-zR)**2 )

							if dist > maior_dist:
								maior_dist = dist

						if maior_dist < raio:
							w.write(rname)
							w2.write("\t"+str(rname)+str(rid))
							# Se o raio eh menor, grave todos os atomos
							for atomo in residuo:
								coord = atomo.get_coord()
								x = coord[0]
								y = coord[1]
								z = coord[2]
								matrix_atomos.append([arquivo,'-',rname,x,y,z,rid])
					except:
						print "Falha ao detectar atomos."

	w.write("\n")
	w2.write("\n")

	print "Gerando a matriz de distancia para: "+arquivo

	nome_final_arquivo = arquivo.split("/")
	
	w3 = open("output_matriz_distancia/"+nome_final_arquivo[1]+"_matriz_distancia_"+str(raio)+".txt",'w')

	# Devido a minha falta de capacidade tecnica com Biopython, irei implementar distancia euclidiana manualmente :(
	for a in matrix_atomos:

		linha_dist = ""

		for b in matrix_atomos:
			dist_a_b = sqrt( (a[3] - b[3])**2 + (a[4] - b[4])**2 + (a[5] - b[5])**2 )
			dist_a_b_str = "%0.2f" %dist_a_b
			linha_dist = linha_dist+dist_a_b_str+"\t"

		w3.write(linha_dist+"\n")

	# Apaga a matrix com atomos
	matrix_atomos = [] 
	w3.close()

w.close()
w2.close()
