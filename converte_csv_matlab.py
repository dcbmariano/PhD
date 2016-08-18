# Converte CSV para MATLAB
# Author: Diego Mariano
# Date: August 18, 2016

import sys

try:
	nome_arquivo = sys.argv[1]
except:
	nome_arquivo = "data3/bgl.csv"

arquivo = open(nome_arquivo).readlines()
nome = nome_arquivo.split(".csv")
out = open(nome[0]+"_out.m","w")

out.write("matrix = [")

for linha in arquivo:
	out.write(linha)
	out.write(";")

out.write("]; ")


print "Successful"