# prepara_ligantes.py

import sys
import glob
import subprocess

try:
	diretorio = sys.argv[1]
except:
	diretorio = 'data/*'

arquivos = glob.glob(diretorio)


for a in arquivos:
	print 'Lendo '+a
	subprocess.call(['./pythonsh', 'scripts/prepare_receptor4.py','-r', a])
