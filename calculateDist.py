# calculateDist.py
# Author: Diego Mariano


from math import sqrt

r = open("result2.txt").readlines()

i = 1

for line in r:
	if i == 3:
		tolerante = line.split(",")
		print line
	if i == 1:
		wild = line.split(",")
		print line
	if i == 2:
		mutante = line.split(",")
		print line

	i += 1

# Calculando distancia euclidiana
tam = len(tolerante)
dist_aux = 0
# distancia mutante

cont = 0
for j in range(tam):
	if cont < 4000:
		dist_aux += (float(tolerante[j]) - float(mutante[j]))**2
	cont += 1


dist_aux = sqrt(dist_aux)

print "Distancia da mutante para tolerante eh: "+str(dist_aux)

# distancia wild

dist_aux = 0
cont = 0
for j in range(tam):
	if cont < 4000:
		dist_aux += (float(tolerante[j]) - float(wild[j]))**2
	cont += 1


dist_aux = sqrt(dist_aux)

print "Distancia da wild para tolerante eh: "+str(dist_aux)