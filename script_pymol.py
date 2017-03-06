#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pymol
from pymol import cmd
from pymol import util
import glob
import os
import sys


# Default
diretorio = "/home/diego/Downloads/pdb_leonardo"
arquivos = os.listdir("/home/diego/Downloads/pdb_leonardo")

# Recebendo valores
for i in range(len(sys.argv)):
	if i % 2  == 1:
		# Diretorio
		if sys.argv[i] == '-d':
			diretorio = sys.argv[i+1]
			arquivos = os.listdir(diretorio)

		# Help
		if sys.argv[i] == '-h' or sys.argv[i] == '--help':
			print "Como executar: "
			print "opt/pymol/pymol -cpr $scripts_path/script_pymol.py -d endereco_diretorio_com_arquivos_pdb"
			print "Nao use '/' ao final do nome da pasta"


# Paranaues do pymol
#pymol.finish_launching()
global pdbPath, outputPath, outputLog, verbose


print "Calculando RMSD"
w = open("/home/diego/result.txt","w")
w.write("\t")


# Carregando PDBs
count = 0
for a in arquivos:

	print "Lendo arquivo: "+a

	if (os.path.exists(diretorio+'/'+a)):
		
		# Carrega arquivo PDB
		cmd.load(diretorio+'/'+a,a) #cmd.load(endereco_completo_mais_nome, apenas_nome)

		w.write(a)
		w.write("\t")

		# Seleciona hidrofobicos
		hidro = "(resn ala+gly+val+ile+leu+phe+met) and "+a
		cmd.select("s"+str(count), hidro)

		# Extrai 
		cmd.extract("h"+str(count),"s"+str(count))

		count += 1

	else:
		print "Arquivo nao encontrado: "+a


# Alinhamentos
print "Construindo alinhamentos"


for i in range(count):

	w.write("\n")
	w.write(arquivos[i]+"\t")

	for j in range(count):

		# Alinhando os mesmos apresenta um erro
		if i != j:
			string = cmd.align("h"+str(i), "h"+str(j))


			w.write(str(string[0]))
			w.write("\t")
		else:
			w.write("0.00")
			w.write("\t")

w.close()
print "Executado com sucesso."