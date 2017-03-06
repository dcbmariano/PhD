# Definindo o melhor modelo (better.py adaptado) **********************************************************************************
# Usa a funcao DOPE para definir o melhor

import os

nomes_pastas = os.listdir("../out")

for nome_pasta in nomes_pastas:

	lista = {}

	print "Lendo diretorio "+nome_pasta

	if nome_pasta != "best" and nome_pasta != "templates_modelos.txt":
		arquivos = os.listdir("../out/"+nome_pasta)

		for nome in arquivos:
				
			n = nome.split(".")
			last = len(n)				

			if n[last-1] == "pdb":
				r = open("../out/"+nome_pasta+"/"+nome)
				arquivo = r.readlines()
				n_linha = 0

				for linha in arquivo:

					if n_linha == 1:
						try:
							lista[nome] = float(linha[40:])
						except:
							print 'Fail. File: '+nome
					if n_linha > 1:
						break

					n_linha += 1

		# SALVA O MELHOR MODELO NA PASTA "BEST"
		better = min(lista,key = lista.get)
		print "Melhor modelo: "+better	
		better = better.replace("|","\|")
		os.system("cd ../out && mkdir best") # Criando pasta best / soh funciona na primeira execucao
		os.system("cd ../out && cp "+nome_pasta+"/"+better+" best/.")
