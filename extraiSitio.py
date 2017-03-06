import os
import glob
from Bio.PDB import *


print "Detectando atomos do sitio ativo baseado em 3VIK.\n"

# Declaracoes
template = "3vik.pdb"

#6 angstrons
#res_sitio = ['A.W.149', 'A.S.193', 'A.T.196', 'A.N.253', 'A.Y.337', 'A.W.374', 'A.E.402', 'A.N.192', 'A.N.335', 'A.H.148', 'A.W.452', 'A.Q.45', 'A.W.444', 'A.E.451', 'A.F.460', 'A.M.207', 'A.N.255', 'A.T.338', 'A.I.254', 'A.Y.273', 'A.L.340', 'A.L.195']

#6.5 angstrons
res_sitio = ['A.W.149', 'A.S.193', 'A.T.196', 'A.M.207', 'A.N.253', 'A.Y.337', 'A.W.374', 'A.E.402', 'A.E.451', 'A.H.148', 'A.N.192', 'A.N.335', 'A.W.452', 'A.W.444', 'A.Q.45', 'A.F.460', 'A.L.195', 'A.I.254', 'A.N.255', 'A.Y.273', 'A.T.338', 'A.F.336', 'A.L.340', 'A.D.199']

#7 angstrons
#res_sitio = ['A.W.149', 'A.S.193', 'A.T.196', 'A.M.207', 'A.N.253', 'A.N.255', 'A.N.335', 'A.Y.337', 'A.W.374', 'A.E.402', 'A.E.451', 'A.W.452', 'A.H.148', 'A.N.192', 'A.F.197', 'A.Q.45', 'A.W.444', 'A.F.460', 'A.L.375', 'A.S.372', 'A.E.458', 'A.P.194', 'A.L.195', 'A.L.252', 'A.I.254', 'A.R.103', 'A.F.450', 'A.D.199', 'A.Y.273', 'A.T.338', 'A.L.340', 'A.F.336', 'A.Q.363']

#7.5 angstrons
#res_sitio = ['A.W.149', 'A.N.192', 'A.S.193', 'A.T.196', 'A.M.207', 'A.N.253', 'A.N.255', 'A.N.335', 'A.Y.337', 'A.W.374', 'A.E.402', 'A.E.451', 'A.W.452', 'A.F.460', 'A.H.148', 'A.L.195', 'A.F.197', 'A.W.444', 'A.Q.45', 'A.F.336', 'A.L.375', 'A.S.372', 'A.E.458', 'A.P.194', 'A.L.252', 'A.I.254', 'A.R.103', 'A.S.251', 'A.A.42', 'A.N.449', 'A.F.450', 'A.R.454', 'A.D.199', 'A.Y.273', 'A.T.338', 'A.L.340', 'A.Q.363', 'A.S.373', 'A.W.256']

#8 angstrons
#res_sitio = ['A.H.148', 'A.W.149', 'A.N.192', 'A.S.193', 'A.T.196', 'A.F.197', 'A.M.207', 'A.N.253', 'A.I.254', 'A.N.255', 'A.N.335', 'A.F.336', 'A.Y.337', 'A.W.374', 'A.L.375', 'A.E.402', 'A.W.444', 'A.E.451', 'A.W.452', 'A.F.460', 'A.Q.45', 'A.P.194', 'A.L.195', 'A.L.252', 'A.R.103', 'A.S.372', 'A.G.404', 'A.S.406', 'A.N.449', 'A.F.450', 'A.R.454', 'A.E.458', 'A.S.251', 'A.A.42', 'A.T.401', 'A.N.403', 'A.S.445', 'A.L.453', 'A.D.199', 'A.Y.273', 'A.N.277', 'A.T.338', 'A.L.340', 'A.Q.363', 'A.S.373', 'A.A.339', 'A.R.353', 'A.L.361', 'A.W.256']

res_outros = []

#res_outros.append([template]+res_sitio)

teste = "149/1bga.pdb"
f149 = glob.glob("149/*.pdb") # posso usar tambem: os.listdir("149")
f21 = glob.glob("pdbs/*.pdb")
fsag = glob.glob("sagarana/*.pdb")
fmod = glob.glob("best/*.pdb")
fdelta = glob.glob("teste_delta/*.pdb")
mutantes = glob.glob("mutantes/*.pdb")
bgl18 = glob.glob("bgl18/*.pdb")
yang = glob.glob("yang/*.pdb")



########################### 21

print "PARTE 1: Realizando alinhamentos com Multiprot. Aguarde."

pcont = 1
for pdb in f21:

	pdb = pdb.replace("|","\|")
	print "%d. %s" %(pcont, pdb)
	pcont += 1

	comando = "export PATH=$PATH:/home/diego/MultiProtInstall \
		&& multiprot.Linux %s %s >> log.txt" %(template,pdb)

	os.system(comando)

	linhas = open("2_sol.res").readlines()

	# Detectando amino

	res_temp = []
	res_temp.append(pdb)

	for linha in linhas:

		if linha != "End of Match List\n":

			res_pdb = linha[0:7]

			# aqui faco um loop pra cada residuo no sitio
			for res in res_sitio:
				if res_pdb == res:
					res_temp.append(linha[8:].strip())
		else:
			break

	res_outros.append(res_temp)
	


print res_outros


print "PARTE 2: Criando arquivos PDB separados"

for residuos in res_outros:

	pdb_name = residuos[0]
	pdb_name = pdb_name.replace("\|","|") 
	print pdb_name

	# Pega nome do diretorio e cria pasta dentro de out
	diretorio_name = pdb_name.split("/")
	os.system("mkdir out/"+diretorio_name[0])

	pdb = open(pdb_name).readlines()

	w = open("out/"+pdb_name+"_temp.pdb","w")
	
	try:
		for linha in pdb:

			for residuo in residuos:
				if residuo != pdb_name:
					if linha[0:7].strip() == "ATOM":
						cadeia = linha[21:22].strip()
						amino3l = linha[17:20].strip()
						amino = Polypeptide.three_to_one(amino3l) # Converte o codigo de 3 letras para 1 letra
						num_amino = linha[23:26].strip()

						compara = residuo.split(".")

						if cadeia == compara[0] and amino == compara[1] and num_amino == compara[2]:
							w.write(linha)
	except:
		print "Falha ao detectar alguns aminoacidos."

	w.close()


print "PARTE 3: Gera arquivo FASTA"

w = open('f21.fasta','w')

for pdb in res_outros:

	i = 0 
	seq = ""

	for res in pdb:
		if i == 0:
			w.write(">"+res+"\n")

		else:
			resi = res.split(".")
			seq += resi[1]
		
		i += 1

	w.write(seq)
	w.write("\n")
