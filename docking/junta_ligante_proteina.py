import glob
import sys


for i in range(len(sys.argv)):
	
	if i % 2  == 1:

		# Ligante
		if sys.argv[i] == '-r':
			pasta_receptores = sys.argv[i+1]
			receptores = glob.glob(pasta_receptores+"/*")

		# Diretorio receptores
		if sys.argv[i] == '-l':
			pasta_ligantes = sys.argv[i+1]
			ligantes = glob.glob(pasta_ligantes+"/*")

		# Help
		if sys.argv[i] == '-h' or sys.argv[i] == '--help':
			print "Run: python junta_ligante_proteina -l pasta_ligantes -r pasta_receptores\nNao use '/' ao final do nome da pasta"


for receptor in receptores:
	
	# abre arquivo pdb
	print "Lendo "+receptor
	arquivo_receptor = open(receptor,"a+")

	aux = receptor.split('/')
	r = aux[len(aux)-1]	#pega ultimo elemento do vetor
	r += 'qt' #simula arquivo pdbqt

	# abre resultado de docking
	for ligante in ligantes:

		aux = ligante.split("docked_") # defini como padrao que a string 'docked_' eh o que difere um nome de outro
		l = aux[len(aux)-1] #retira o nome 

		if l == r:
			print "OK"
			arquivo_ligante = open(ligante).readlines()

			for linha in arquivo_ligante:
				if linha != "MODEL 2\n": # grava apenas model 1
					arquivo_receptor.write(linha)
				else:
					break

print "Concluido com sucesso."


