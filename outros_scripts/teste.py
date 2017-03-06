
from Bio.PDB import *
from math import sqrt
import sys
import warnings
warnings.filterwarnings("ignore") # esconde erros do pdb


from BioSVD import *

assinaturas = [[1,2,3],[4,5,6],[7,8,9]]

# Preenche auxiliar vazio
assinaturas_aux = []
for i in range(len(assinaturas[0])):
	assinaturas_aux.append(list())

# Transpoe matriz
for i in range(len(assinaturas[0])):
	for j in range(len(assinaturas)):
		assinaturas_aux[i].append(assinaturas[j][i])

print assinaturas_aux

'''
parser = PDBParser()
estrutura = parser.get_structure('PDB','out_149/1bga.pdb_out_15.pdb')

for residuo in estrutura[0]['A']:
	print residuo.id[1]

teste = ".B.123"

amino = teste.split(".")

print amino[0]
print amino[1]
print amino[2]

w = open("149/1bga.pdb").readlines()

for linha in w:
	if linha[0:4] == "ATOM":
		print linha
		print linha[21:22]
		print linha[17:20]
		print linha[23:26]
		print linha[30:38]
		print linha[38:46]
		print linha[46:54]
		print linha[13:15]


ATOM   1677  O   PRO A 209     -24.199 -38.571 -31.632  1.00  9.07           O  

COMECO -1
FINAL  -2

atom[0:3]
cadeia[21:21]
amino[17:19]
num_amino[23:25]
x[30:37]
y[38:45]
z[46:53]
'''