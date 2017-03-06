c1 = ["222 THR", "240 ALA", "300 THR", "302 ASN"]
c1 = ["222 SER", "240 ALA", "300 THR", "302 SER"]


combinacoes = []

for i1 in c1:	
	print "python mutate.py ambgl18_sitio "+i1+" ' ' > log.txt"

	j1 = i1.split(" ")
	combinacoes.append([i1,"","",""])

	for i2 in c1:
		j2 = i2.split(" ")

		atual2 = [i1,i2,"",""]
		atual2.sort()

		if atual2 not in combinacoes:
			combinacoes.append(atual2)
			if i1 != i2:
				print "python mutate.py ambgl18_sitio"+j1[1]+j1[0]+" "+i2+" ' ' > log.txt"

		for i3 in c1:
			j3 = i3.split(" ")

			atual3 = [i1,i2,i3,""]
			atual3.sort()

			if atual3 not in combinacoes:
				combinacoes.append(atual3)
				if i1!=i2 and i1!=i2 and i1!=j3 and i2 != i3:
					print "python mutate.py ambgl18_sitio"+j1[1]+j1[0]+j2[1]+j2[0]+" "+i3+" ' ' > log.txt"

			for i4 in c1:
				j4 = i3.split(" ")

				atual4 = [i1,i2,i3,i4]
				atual4.sort()
				
				if atual4 not in combinacoes:
					combinacoes.append(atual4)
					if i1!=i2 and i1!=i2 and i1!=j3 and i1 != i4 and i2 != i3 and i2!=i4 and i3!=i4:
						print "python mutate.py ambgl18_sitio"+j1[1]+j1[0]+j2[1]+j2[0]+j3[1]+j3[0]+" "+i4+" ' ' > log.txt"

'''



#for i in c2:
#	print "python mutate.py ambgl18_sitio "+i+" ' ' > log.txt"
'''