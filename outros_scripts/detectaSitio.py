import os
import glob

print "Detectando sitios ativos baseado em 4MDO.\n"

# Declaracoes
template = "4mdo.pdb"
teste = "149/1bga.pdb"
f149 = glob.glob("149/*.pdb") # posso usar tambem: os.listdir("149")
f21 = glob.glob("21/*.pdb")

# Predicao 
w = open("prediction.txt","w")
i = 0
j = 0

# Executa
print "Lendo pasta 149. Aguarde."
for pdb in f149:

	i += 1

	comando = "export PATH=$PATH:/home/diego/MultiProtInstall \
		&& multiprot.Linux %s %s >> log.txt" %(template,pdb)

	os.system(comando)

	linhas = open("2_sol.res").readlines()

	# Detectando GLU1
	GLU1 = ""
	for linha in linhas:
		if linha[0:7] == "A.E.166": #linha[0:10] = 112

			GLU1 = linha[9:16]			
			break

	# Detectando GLU2
	GLU2 = ""
	for linha in linhas:
		if linha[0:7] == "A.E.377":
			GLU2 = linha[9:16]
			break

	result = "%d\t%s\t%s\t%s\n" %(i,pdb,GLU1,GLU2)
	print result
	w.write(result)

	if GLU1 != "" and GLU2 != "":
		j += 1


# Predicao
acertos = j/i
print "Total de acertos: %d / %d = " %(j,i)
print int(acertos)


########################### 21

i = 0
j = 0

# Executa
print "Lendo pasta 21. Aguarde."
for pdb in f21:

	i += 1

	comando = "export PATH=$PATH:/home/diego/MultiProtInstall \
		&& multiprot.Linux %s %s >> log.txt" %(template,pdb)

	os.system(comando)

	linhas = open("2_sol.res").readlines()

	# Detectando GLU1
	GLU1 = ""
	for linha in linhas:
		if linha[0:7] == "A.E.166": #linha[0:10] = 112

			GLU1 = linha[9:16]			
			break

	# Detectando GLU2
	GLU2 = ""
	for linha in linhas:
		if linha[0:7] == "A.E.377":
			GLU2 = linha[9:16]
			break

	result = "%d\t%s\t%s\t%s\n" %(i,pdb,GLU1,GLU2)
	print result
	w.write(result)

	if GLU1 != "" and GLU2 != "":
		j += 1

w.close()

# Predicao
acertos = j/i
print "Total de acertos: %d / %d = " %(j,i)
print int(acertos)


print "Concluido"

w.close()