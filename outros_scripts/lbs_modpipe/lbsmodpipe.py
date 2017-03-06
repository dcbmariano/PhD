# Modelagem em larga escala
# -*- coding: utf-8 -*-

import os
import sys
import glob
from Bio import *
from Bio import SeqIO
from Bio.Blast.Applications import *
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.PDB import *
import warnings
from random import randint
import psutil
from time import sleep
from multiprocessing import Process,Array

warnings.filterwarnings('ignore') # Reduz mensagens de erro do Bio.PDB


total_modelos = 100

print "***************************************************************************************************************"
print "LBS-MODELLER v0.3"
print "\nAVISO: [...] indica que o sistema esta aguardando liberacao do processador.\n"
print "***************************************************************************************************************\n"

# DECLARACOES ***************************************************************************************************************
# Databases
seqdb = "db/seq/db.fasta"
pdbdb = glob.glob("db/pdb/*.pdb")

#seqs
try:
	seqs_fasta = sys.argv[1]
except:
	seqs_fasta = "data/seqs.fasta"


# Funcoes de paralelizacao *****************************************************************************************************

procs =  6 #countCPU()
cont = 0 # contador para paralelizacao
controle = [False]*procs
arr = Array("i",controle)

def countCPU():
	#return: Número de CPUs do computador
    tCPUs = 0
    try:
        tCPUs = psutil.cpu_count()
    except:
        tCPUs = psutil.NUM_CPUS
    return tCPUs if tCPUs ==1 else tCPUs -1

def hasCoreFree(arrayControle):
    #return the first core free id.
    index = -2
    try:
        index = arrayControle[:].index(False)
    except:
        index = -1
    return index


# Funcao que executa o Modeller paralelizado ************************************************************************************
def modela(arrayControle,seqModel,coreId, nome_pasta):

	print "\nIniciando a modelagem de %s no core %s." %(seqModel,coreId)

	#comando = "cd out/"+nome_pasta+" && /home/diego/bin/modpipe-2.2.0/ext/mod/bin/modSVN.in modeller_tmp.py"

	#Atualizacao servidor BIOINFO
	comando = "cd out/"+nome_pasta+" && /programs/modeller/bin/mod9.17 modeller_tmp.py"

	#os.system(comando)
	print "\n[ Modelagem de %s foi concluida com sucesso ]\n" %(seqModel)

	arrayControle[coreId] = False


# Lendo cada sequencia fasta que queremos MODELAR ********************************************************************************
j = 1
diretorio = os.path.dirname(os.path.realpath(__file__)) #busca o diretorio atual
templates_modelos = open("out/templates_modelos.txt", "w")

for i in SeqIO.parse(seqs_fasta,'fasta'):
	
	templates_modelos.close() # Permite que o arquivo seja salvo antes da finalizacao do script
	templates_modelos = open("out/templates_modelos.txt", "a+") 

	print "[ ... ]\n"
	coreFree = hasCoreFree(arr)

	if coreFree > -1 :

		arr[coreFree] = True

		# Definindo o template *************************************************************************************************************************
		print "Analisando seq: %d" %j

		# Gerando id randomico
		rand_id = str(randint(0,999999))

		# Criando arquivo temporario com seq atual
		tmp_name = "tmp/tmp_%s.fasta" %rand_id
		tmp_blast = "tmp/blast_%s.out" %rand_id
		tmp_align = "tmp/align_%s.pir" %rand_id
		tmp = open(tmp_name,'w')
		tmp.write(">"+i.id+"\n"+str(i.seq))
		tmp.close()

		# blastp
		blastp = NcbiblastpCommandline( query = tmp_name, subject = seqdb, outfmt = "'6 qseqid sseqid pident score slen sstart send'" , out = tmp_blast)
		stdout, stderr = blastp()

		# blastp sem biopython
		#blastp = "blastp -out %s -outfmt '6 qseqid sseqid pident score slen sstart send' -query %s -subject %s" %(tmp_blast, tmp_name, seqdb)

		# SERVIDOR BIOINFO UPDATE
		#blastp = "/programs/ncbi-blast-2.2.31+/bin/blastp -out %s -outfmt '6 qseqid sseqid pident score slen sstart send' -query %s -subject %s" %(tmp_blast, tmp_name, seqdb)
		#os.system(blastp)
		

		blast_result = open(tmp_blast).readlines()

		max_score = 0

		linha_max = []

		for linha in blast_result:
			
			#Pegando o maior score
			info = linha.split("\t")
			score = int(info[3])
			if(score > max_score):
				max_score = score
				linha_max = info

		print "Template com maior score: "
		print linha_max

		seq = linha_max[0]
		template = linha_max[1]
		identidade = float(linha_max[2])
		split_template = template.split("|")
		split_pdbid  = split_template[0].split(":")
		pdbid_template = split_pdbid[0].lower()
		cadeia_template = split_pdbid[1]


		# Necessario estimar coverage ****************************************************************************
		# Formula: final do subject menos inicio do subject sobre o tamanho do subject 
		coverage = 100 * ( float(linha_max[6]) - float(linha_max[5]) ) / float(linha_max[4])


		# PARAR MODELAGEM QUANDO:
		if identidade < 25:
			print "IDENTIDADE"
			print identidade
			print "***************** AVISO: IDENTIDADE BAIXA - Ignorando modelagem *****************"
			arr[coreFree] = False # Libera o processador
			j += 1	#incrementa o contador de modelos
			continue

		if coverage < 25:
			print "COBERTURA: "
			print coverage
			print "***************** AVISO: COBERTURA BAIXA - Ignorando modelagem *****************"
			arr[coreFree] = False # Libera o processador
			j += 1	#incrementa o contador de modelos
			continue
		

		templates_modelos.write(linha_max[0]+"\t"+linha_max[1]+"\t"+linha_max[2]+"\t"+linha_max[3]+"\t"+str(coverage)+"\n")


		# Criando pasta do modelo
		nome_pasta = seq.replace("|","_")
		nome_pasta = nome_pasta.replace(":","_")
		os.system("mkdir out/"+nome_pasta)


		# Criando arquivo temporario para alinhamento **************************************************************************************************
		align_tmp = open(tmp_align,"w")
		align_tmp.write(">"+i.id+"\n"+str(i.seq)+"\n")

		# PEGAR SEQUENCIA PDB

		align_tmp.write(">"+pdbid_template+"\n")

		parser = PDBParser()
		estrutura = parser.get_structure("PROT","db/pdb/"+pdbid_template+".pdb")
		for residuo in estrutura[0][cadeia_template]:
			if is_aa(residuo):
				try:
					resname = Polypeptide.three_to_one(residuo.resname)
					align_tmp.write(resname)
				except:
					print "Falha ao detectar residuo: "+str(residuo.resname)+". Residuo sera ignorado."

		
		align_tmp.close()

		# Alinhamento com CLUSTALW
		clustalw = ClustalwCommandline("clustalw",infile = tmp_align, output = "pir")
		clustalw()

		# sem biopython
		#clustalw = "clustalw2 -infile=%s -output=pir" %tmp_align


		# SERVIDOR BIOINFO UPDATE
		#clustalw = "/home/diegomariano/bin/clustalw-2.1/src/clustalw2 -infile=%s -output=pir" %tmp_align
		#os.system(clustalw)


		# Configurando alinhamento para o formato MODELLER
		# Alinhamento salvo em: tmp/align_tmp.pir
		align = open(tmp_align).readlines()
		align_mod = open("out/"+nome_pasta+"/align_tmp_mod.pir","w")

		tipo = 0 # 0 sequencia; 1 pdb
		

		#obtendo o primeiro e ultimo residuo
		menor = 1000000
		maior = 0
		pdb_atual = open("db/pdb/"+pdbid_template+".pdb")
		for linha in pdb_atual:
			tipo_linha = linha[0:5].strip()
			cadeia = linha[21] 

			if cadeia == cadeia_template and tipo_linha == "ATOM":
				resid = int(linha[22:26])

				if resid > maior:
					maior = resid
				if resid < menor:
					menor = resid

		for a in align:
			
			if a[0] == ">" and tipo == 0:
				align_mod.write(a)
				align_mod.write("sequence:%s:1:A:%d:A::::" %(seq, len(i.seq) ) )			
				tipo += 1

			elif a[0] == ">" and tipo == 1:
				align_mod.write(">P1;%s\n" %pdbid_template)
				align_mod.write("structure:%s:%d:%s:%d:%s::::" %(pdbid_template, menor, cadeia_template, maior, cadeia_template) )

			else:
				align_mod.write(a)


		align_mod.close()


		#  Modeller

		# Criando script temporario
		script_tmp = open("out/"+nome_pasta+"/modeller_tmp.py","w")

		script_tmp.write("# -*- coding: utf-8 -*-\n")
		script_tmp.write("from modeller import *\n")
		script_tmp.write("from modeller.automodel import *\n")
		script_tmp.write("log.verbose()\n")
		script_tmp.write("env = environ()\n")
		#script_tmp.write("env.io.atom_files_directory = ['"+diretorio+"/out/"+nome_pasta+"']\n")
		script_tmp.write("env.io.atom_files_directory = ['/home/diegomariano/modelagem/lbs_modpipe/out/"+nome_pasta+"']\n")
		
		script_tmp.write("env.io.hetatm = True\n")
		script_tmp.write("env.io.water = True\n")
		script_tmp.write("a = automodel(env, alnfile = 'align_tmp_mod.pir', knowns = '"+pdbid_template+"', sequence = '"+seq+"')\n")
		script_tmp.write("a.starting_model= 1\n")
		script_tmp.write("a.ending_model  = "+str(total_modelos)+"\n")
		script_tmp.write("a.make()")

		script_tmp.close()

		#copiando pdb para diretorio
		os.system("cp db/pdb/"+pdbid_template+".pdb out/"+nome_pasta+"/.")


		# Excutando modeller ********************************************************************************************************************
		p = Process(target=modela, args=(arr, nome_pasta, coreFree, nome_pasta)) # corrigir essa redundancia depois
		p.start()		
		

		print "" #insere uma quebra de linha
		j += 1	#incrementa o contador de modelos

	else:
		while(coreFree < 0):
			print "[ ... ]"
			sleep(10) # Aguarde 10 seg e verifique novamente se o proc esta livre
			coreFree = hasCoreFree(arr)
	 



# Limpando diretorio tmp
#os.system("cd tmp && rm *")
print "Executado com sucesso."