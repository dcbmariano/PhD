#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os
import psutil
from time import sleep
from multiprocessing import Process,Array




def modela(arrayControle,inFilePDB,coreId):

	print "Realizando a modelagem de %s no core %s" %(inFilePDB,coreId)
	print "..."
	for i in range(10000000):
		j = i # soh pra travar um pouco
	print "Concluida"
	arrayControle[coreId] = False


def run(blockDir):
    
    filesToProcess = []
    for i in range(100000):
    	filesToProcess.append(i)

    tam =  len(filesToProcess)
    procs = 4 #countCPU()
    ct = 0    #contador usado para iterar na lista que será processada

    controle = [False]*procs
    arr = Array("i",controle)

    #while que irá iterar na lista a ser processada
    while ct < tam:
    	print "Aguardando"
        coreFree = hasCoreFree(arr)
        if coreFree > -1 :
            arr[coreFree] = True
            p = Process(target=modela, args=(arr, filesToProcess[ct], coreFree))
            p.start()
            print("Arquivo %s number %d of %d ")%(filesToProcess[ct],ct,tam)
            ct+=1
        else:
            sleep(10)


if __name__=="__main__":
    run("pdb/")

print "End"