from Bio.PDB import *
import sys
from math import sqrt


print "ETAPA 1"
print "Calculando centroide sitio"

# lendo prediction.txt
pred_txt = open("prediction.txt").readlines()
w = open("centroide.txt","w")
for linha in pred_txt:
	
	celulas = linha.split("\t") 

	nome_arquivo = celulas[1]
	
	GLU1 = celulas[2].strip()
	GLU2 = celulas[3].strip()

	# Valida se existe GLU1 e GLU2
	if GLU1 != "\n" and GLU2 != "\n" and GLU1 != "" and GLU2 != "":
		
		print nome_arquivo
		arquivo_atual = open(nome_arquivo).readlines()

		# Detectando cadeia, aminoacido e num_aminoacido que deve ser busca
		# GLU1
		info_amino1 = GLU1.split(".")
		cadeia1 = info_amino1[0].strip()
		amino1 = info_amino1[1].strip()
		num_amino1 = info_amino1[2].strip()

		# GLU2
		info_amino2 = GLU2.split(".")
		cadeia2 = info_amino2[0].strip()
		amino2 = info_amino2[1].strip()
		num_amino2 = info_amino2[2].strip()

		atom1 = ""
		atom2 = ""

		# Capturando X, Y, Z de GLU1 e GLU2
		for linha2 in arquivo_atual:
			if linha2[0:4] == "ATOM":

				# Extraindo informacoes de cada linha
				cadeia = linha2[21:22].strip()
				amino3l = linha2[17:20].strip()
				amino = Polypeptide.three_to_one(amino3l) # Converte o codigo de 3 letras para 1 letra
				num_amino = linha2[23:26].strip()
				atom = linha2[13:15].strip()
				x = linha2[30:38].strip()
				y = linha2[38:46].strip()
				z = linha2[46:54].strip()


				#GLU1
				if cadeia == cadeia1 and amino == amino1 and num_amino == num_amino1 and (atom == 'CA' or atom == 'OE'):
					print "GLU1 detectado."
					atom1 = atom
					x1 = x
					y1 = y
					z1 = z

				#GLU2
				if cadeia == cadeia2 and amino == amino2 and num_amino == num_amino2 and (atom == 'CA' or atom == 'OE'):
					print "GLU2 detectado."
					atom2 = atom
					x2 = x
					y2 = y
					z2 = z


		# Calcula o centroide
		xc = (float(x2)+float(x1))/2
		yc = (float(y2)+float(y1))/2
		zc = (float(z2)+float(z1))/2

		dist_euclidiana = sqrt( (float(x2) - float(x1))**2 + (float(y2) - float(y1))**2 + (float(z2) - float(z1))**2)

		# Grava resultado - padrao: nome_arquivo GLU1 x1 y1 z1 GLU2 xx y2 z2
		result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n" %(nome_arquivo,GLU1,atom1,x1,y1,z1,GLU2,atom2,x2,y2,z2,xc,yc,zc,dist_euclidiana)
		w.write(result)

	# Reseta os valores de GLU1 e GLU2
	GLU1 = GLU2 = ""

w.close()

print "\nConcluido com sucesso."
print "Confira o resultado em centroide.txt"
print "Consta em centroide.txt: nome_arquivo,GLU1,atom1,x1,y1,z1,GLU2,atom2,x2,y2,z2,xc,yc,zc,dist_euclidiana\n"