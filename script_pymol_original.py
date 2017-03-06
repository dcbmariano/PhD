# como executar:
# opt/pymol/pymol -cpr $scripts_path/generate_cluster_images.py


#!/usr/bin/python2.7
from __future__ import with_statement
import __main__
import os
import os.path
import sys
import getopt
import pprint
import simplejson as json
import re
import pymol
from pymol import cmd
from pymol import util

pymol.finish_launching()


def generateSuperposition(cluster,clusterNumber):
  global pdbPath, outputPath, outputLog, verbose
    
  logResult = ""
  
  if (verbose == True): print "Generating the superposition to the cluster: " + clusterNumber + "." 
  
  ligandSelection = ""
  
  for entry in cluster:
    arrayAux = entry.split(":") 
         
    pdb = tryGetArrayElement(arrayAux,0)    
    chain = tryGetArrayElement(arrayAux,1)    
    ligand = tryGetArrayElement(arrayAux,2)
    ligandNumber = tryGetArrayElement(arrayAux,3)

    # If one of the three keys have not been passed, next.
    if (pdb == None or chain == None or ligand == None or ligandNumber == None):      
      if (verbose == True): print "A PDB or a Chain or a Ligand were not defined to the entry: " + entry + "."
      continue

    # Define the PDB filename.         
    filename = pdbPath + "/" + pdb + ".chain" + chain + ".align.pdb"    
  
    # Load the PDB file if it exist.
    if ( os.path.exists(filename) ):
      cmd.load(filename,entry)           
      if (ligandSelection != ""): ligandSelection += " or "
      ligandSelection += "(resn " + ligand + " and res " + ligandNumber + " and chain " + chain + ")"          
      logResult += entry + "\tOK\n"            
    else:       
      if (verbose == True): print "File " + filename + " do not exist."
      logResult += entry + '\t' + 'FILE_NOT_FOUND\n'           
      continue # if the file do not exist, next.
  
  if (ligandSelection == ""): return

  cmd.hide("everything")
  cmd.select("ligands",ligandSelection)
  cmd.show("sticks","ligands")
  # Show spheres. With a scale equal to 0.1, the spheres only appears when sticks do not exist, i.e., when there are Ions in the PDB file.
  cmd.show("spheres","ligands")  
  cmd.set("sphere_scale", 0.1,"ligands")
  
  cmd.bg_color("white")
  
  cmd.orient("ligands")
  cmd.zoom("ligands",0,0,1)
 
  util.cbag("all") # Color atoms by element.

  cmd.png(outputPath + "/cluster" + clusterNumber + ".png",2200,1500,300,1)  
  cmd.delete("all")
  
  # Write the results to the log file if it was required by the user.
  if (outputLog != None):
    with open(outputLog, "a") as f:
      f.write(logResult)