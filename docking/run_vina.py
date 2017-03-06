#run_vina.py
#executa o vina para todo um diretorio

import sys
import glob
import subprocess


aresta = "15.00"

for i in range(len(sys.argv)):
	
	if i % 2  == 1:

		# Ligante
		if sys.argv[i] == '-l':
			ligante = sys.argv[i+1]

		# Diretorio receptores
		if sys.argv[i] == '-r':
			diretorio = sys.argv[i+1]

		# arquivo com centro 
		if sys.argv[i] == '-c':
			centroide = sys.argv[i+1]

		if sys.argv[i] == '-a':
			aresta = sys.argv[i+1]


print "Realizando docking\n"
arquivos = glob.glob(diretorio+"/*.pdbqt")


for a in arquivos:

	print "Lendo arquivo: "+a
	# retorna centroide
	centroide_arquivo = open(centroide).readlines()
	
	for linha in centroide_arquivo:

		celula = linha.split("\t")
		arquivo_comparacao = celula[0].split("/")
		arquivo_lido = a.split("/")
		arquivo_comparacao[1] += 'qt'

		if arquivo_comparacao[1] == arquivo_lido[1]:
			x = celula[11]
			y = celula[12]
			z = celula[13]

	
	# Cria um arquivo de configuracao temporario
	config = "receptor = "+a
	config += "\nligand = "+ligante
	config += "\ncenter_x = "+x
	config += "\ncenter_y = "+y
	config += "\ncenter_z = "+z
	config += "\nsize_x = "+aresta+"\nsize_y = "+aresta+"\nsize_z = "+aresta
	config += "\nout = docked/docked_"+arquivo_lido[1]
	config += "\nlog = docked/docked_"+arquivo_lido[1]+".log"
	config += "\nexhaustiveness = 150"
	config += "\nnum_modes = 10"

	# Grava resultado
	temp = open('config.txt','w')
	temp.write(config)
	temp.close()
	print "Gravado com sucesso."


	# Excutando o vina
	print "Execuntando docking para: "+a
	subprocess.call(["./dock/vina", "--config", "config.txt"])

print "Concluido com sucesso."