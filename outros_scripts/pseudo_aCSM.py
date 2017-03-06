from Bio.PDB import *
import sys
import os
import warnings
warnings.filterwarnings("ignore") 
import numpy as np

from BioSVD import *


print "aCSM"

parser = PDBParser()

# Preenchendo matriz de assinatura
assinatura = []
for i in range(0,150):
	assinatura.append(0)

assinaturas = []

print "Lendo matrizes de distancia."

#pasta = os.listdir('output_matriz_distancia')
#pasta = os.listdir('teste')
pasta = os.listdir('output_9')


for arquivo in pasta:

	#pdb = open('output_matriz_distancia/'+arquivo).readlines()
	#pdb = open('teste/'+arquivo).readlines()
	pdb = open('output_9/'+arquivo).readlines()


	for linha in pdb:
		celulas = linha.split("\t")

		for celula in celulas:

			# "For" deve percorrer de 0 a 30 angstrons, somando 0,5 a cada rodada, ou seja 60 itens
			for i in range(0,150):
				if (celula != "\t" and celula != "\n") and float(celula) > i/5 and float(celula) <= (i/5)+0.2:
					assinatura[i] += 1
	
	print "Obtendo assinatura do arquivo: "+arquivo
	
	assinaturas.append(assinatura)

	# Zera assinatura
	assinatura = []
	for i in range(0,150):
		assinatura.append(0)


# Necessario transpor a matriz assinaturas
# Preenche auxiliar vazio
assinaturas_aux = []
for i in range(len(assinaturas[0])):
	assinaturas_aux.append(list())

# Transpoe matriz
for i in range(len(assinaturas[0])):
	for j in range(len(assinaturas)):
		assinaturas_aux[i].append(float(assinaturas[j][i]))

#print assinaturas_aux

# Rodando SVD
assinaturas_aux = np.matrix(assinaturas_aux)
[U,S,V] = BioSVD.svds(assinaturas_aux)
m2 = BioSVD.extractFactor(S,V,2)
# Testing function "factor"
s = BioSVD.factor(S,'plot')

BioSVD.plot2(m2,'plot',pasta)