import glob
import sys
import subprocess
import os

try:
	arquivos = glob.glob(sys.argv[1]+"/*")
except:
	arquivos = glob.glob("receptores/*.pdb")

for a in arquivos:
	query = "grep -vwE \"HETATM\" %s > %s_" %(a,a)
	os.system(query)