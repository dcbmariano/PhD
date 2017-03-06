# Extrai aminoacidos a partir de ligante
# Regra: sao extraidos aminoacidos cuja distancia para o caborno alfa eh inferir a 6, 9, 12 ou 15 angstrons

from Bio.PDB import *
from math import sqrt
import sys
import warnings
warnings.filterwarnings("ignore") # esconde erros do pdb


try:
	raio = sys.argv[1]
except:
	raio = 8

ligante = "CBI"
atomos_ligante = []
aminoacidos_raio = []


pdb = open("3vik.pdb").readlines()

# capturando atomos ligante
for linha in pdb:

	busca = linha[17:20].strip()
	tipo = linha[0:6].strip()

	if busca == ligante and tipo == "HETATM":

		# capturando informacoes basicas
		atom = linha[13:15].strip()
		x = linha[30:38].strip()
		y = linha[38:46].strip()
		z = linha[46:54].strip()

		atomos_ligante.append([atom,x,y,z])


#para cada atomo do ligante, construa uma matriz de distancia

for atomo in atomos_ligante:

	# captura aminoacidos no raio
	for linha in pdb:

		tipo = linha[0:6].strip()

		if tipo == "ATOM":
		
			cadeia = linha[21:22].strip()
			amino3l = linha[17:20].strip()
			amino = Polypeptide.three_to_one(amino3l) # Converte o codigo de 3 letras para 1 letra
			num_amino = linha[23:26].strip()
			
			xr = linha[30:38].strip()
			yr = linha[38:46].strip()
			zr = linha[46:54].strip()

			# calcula distancia
			dist = sqrt( (float(atomo[1]) - float(xr) )**2 + (float(atomo[2]) - float(yr) )**2 + (float(atomo[3]) - float(zr) )**2 )

			if dist < raio:
				aa = cadeia+"."+amino+"."+num_amino

				if aa not in aminoacidos_raio:
					aminoacidos_raio.append(aa)

print "Aminoacidos num raio de %0.1f angstrons: \n" %raio
print aminoacidos_raio
print "\n"