from copy import deepcopy
from os import listdir
from os.path import isfile, join, dirname
import os
import sys
import subprocess
from datetime import datetime
from collections import defaultdict
from collections import namedtuple
import functions
import math
import time
import psutil

## Get the current directory
filepath = os.path.dirname(os.path.realpath(__file__)) + "/"
hlapath = filepath+ "hla"
allfiles = [f for f in listdir(hlapath) if isfile(join(hlapath, f)) and f.endswith(".txt")]
allfiles = set(allfiles)
length = 11
tab = "\t"
#print allfiles

##############################
##############################
syfpeithiList= set(functions.readSyphethi())
IEDBList = set(functions.readIEDB())
netmhcpanList = set(functions.readnetmhcpan())

file_hlas_map =dict();
for file in allfiles:
     file_hlas_map[file] = list();
     with open(hlapath+ "/"+ file) as f:
          for line in f:
               hlas = line.strip().split("\t")
               file_hlas_map[file].extend(list(set(hlas[1:len(hlas)])))

outputFolder =  str(datetime.utcnow()).replace(" ","-").replace(":","-").replace(".","-") + "/"
outputFilePath = filepath + "output/" + outputFolder
isomers = ['15','17','19','21']
#print file_hlas_map



#############################################
# Core Processing
#############################################


procs = []
key_ = []
value_ = []

patientSet = set()
patientHlaMAP = defaultdict(list)
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     patientSet.add(patientID)
     hlaSet = set()
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]

          syfpeithiStr = "False"
          IEDBStr = "False"
          netmhcpanStr = "False"
          if hla in syfpeithiList:
               syfpeithiStr = "True"

          value_.append(hla)     

          if hla in IEDBList:
               IEDBStr  = "True"
          if hla in netmhcpanList:
                netmhcpanStr = "True"
          updatedHLA = ""
          
          if IEDBStr == "False":
               hla =  functions.getClosestHLA(hla, IEDBList)
               key_.append(hla)
          else:
               key_.append(hla)

          if hla not in hlaSet:
               hlaSet.add(hla);
               patientHlaMAP[patientID].append(hla)

          if not os.path.exists(outputFilePath + patientID):
               os.makedirs(outputFilePath + patientID)


for index, patient in enumerate(patientSet):
     processPatient = True
     while processPatient and index < len(patientSet)  :
          if psutil.cpu_percent() < 85.0:
               hlas = patientHlaMAP[patient]
               print patient + "-- " + str(hlas)
               for hla in hlas:
                    #if not os.path.exists(outputFilePath+ patient +"/" + hla.replace(":","-") + "/"):
                    os.makedirs(outputFilePath+ patient +"/" + hla.replace(":","-") + "/")          
                    for num in isomers:
                         filename = "TCGA-" + patient +"_Varscan_variants_filter.pass."+ str(num) +".peptide"
                         functions.writeInputFile(filepath + "/peptides/", filename, patient,str(num))                         
                         proc = subprocess.Popen([sys.executable, 'epitope_new.py', hla, num , filename,  patient, outputFolder, syfpeithiStr, IEDBStr, netmhcpanStr])
                         procs.append(proc)
               processPatient = False                 
               time.sleep(5)        
          else:
               time.sleep(60)
     print patient + tab + str(patientHlaMAP[patient])
          

newOldHLAMap = dict(zip(key_,value_))
for proc in procs:
     proc.wait()

for patient in patientSet:
     fwrite = open(outputFilePath+ patient + "/"+patient+ '.tsv','w')
     fwrite.write(functions.getHeaderText()) 
     hlas = patientHlaMAP[patient]

     for hla in hlas:
          IEDB_transcriptMap_MT = defaultdict(list)
          IEDB_transcriptMap_WT = defaultdict(list)
          IEDB_TranscriptMap_SameSeq = defaultdict(list)
          transcript_SET_G = set()
          dataTuple = namedtuple("dataTuple", ("mer","data"))
          
          for num in isomers:

               transcript_SET, mutant_Map, wildType_Map = functions.getMutantWildTypeData(filepath, patient, num)
               peptideLenStr = "."+ str(num) + ".txt"
               outputFolder = functions.getFilePath(outputFilePath,hla,patient)
               outputIEDBFile= outputFolder + "output_IEDB" + peptideLenStr
               IEDBMap, IEDBTuples =  functions.getMapwithValuesIEDB(outputIEDBFile)
               
               ###################################################
               # Transcript - peptide - min score 
               ###################################################
               
               peptideSet = set()  
               transcript_SET_G = transcript_SET  
               sameSeqTupleList= list()
               for t in transcript_SET:
                    mutantLowestScoreTuple = functions.getLowestScore(mutant_Map[t][0], IEDBTuples )
                    IEDB_transcriptMap_MT[t].append( dataTuple( mer = num, data = mutantLowestScoreTuple))
                    IEDB_transcriptMap_WT[t].append(dataTuple( mer= num, data = functions.getLowestScore(wildType_Map[t][0], IEDBTuples )))
                    IEDB_TranscriptMap_SameSeq[t].append(dataTuple( mer= num, data = functions.getSameSeqScore(mutantLowestScoreTuple, wildType_Map[t][0], IEDBTuples )))
                    peptideSet.update(functions.getPeptides(IEDB_transcriptMap_MT[t],num))
                    peptideSet.update(functions.getPeptides(IEDB_transcriptMap_WT[t],num) )
                    peptideSet.update(functions.getPeptides(IEDB_TranscriptMap_SameSeq[t],num) )
               
          for transcript in transcript_SET_G:
               fwrite.write(patient+  tab + hla+ tab + transcript )
               merList  = list()
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_transcriptMap_MT[transcript], num)
                    merList.append((num,float(bindingScore), peptide))
                    if num == "15":
                         fwrite.write(tab  + "MT") 
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))          
               fileData = functions.getFinalMer(merList)
               mer = fileData[0]
               score = fileData[1]
               peptide = fileData[2]
               fwrite.write(tab+ peptide + tab + mer + tab+ str(score))

               merList  = list()
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_transcriptMap_WT[transcript], num)
                    merList.append((num,float(bindingScore),peptide))
                    if num == "15":
                         fwrite.write(tab  + "WT") 
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
               fileData = functions.getFinalMer(merList)
               mer = fileData[0]
               score = fileData[1]
               peptide = fileData[2]
               fwrite.write(tab+ peptide + tab + mer + tab+ str(score))

               merList  = list()
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_TranscriptMap_SameSeq[transcript], num)
                    merList.append((num,float(bindingScore),peptide))
                    if num == "15":
                         fwrite.write(tab  + "Same_Seq")      
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
               fileData = functions.getFinalMer(merList)
               mer = fileData[0]
               score = fileData[1]
               peptide = fileData[2]
               fwrite.write(tab+ peptide + tab + mer + tab+ str(score))

               newHla = functions.getNewHLA(newOldHLAMap, hla)
               fwrite.write(tab + newHla+"\n")
                         
   
     fwrite.close()


#############################################
# End of Program
#############################################     


