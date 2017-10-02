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
patient_hla_file = ()
hlaSet = set()
patientHlaMAP = defaultdict(list)
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]

          syfpeithiStr = "False"
          IEDBStr = "False"
          netmhcpanStr = "False"
          if hla in syfpeithiList:
               syfpeithiStr = "True"

          key_.append(hla)     

          if hla in IEDBList:
               IEDBStr  = "True"
          if hla in netmhcpanList:
                netmhcpanStr = "True"
          updatedHLA = ""
          
          if IEDBStr == "False":
               hla =  functions.getClosestHLA(hla, IEDBList)
               value_.append(hla)
          else:
               value_.append(hla)

          if hla not in hlaSet:
               hlaSet.add(hla);
               patientHlaMAP[patientID].append(hla)

          if not os.path.exists(outputFilePath + patientID):
               os.makedirs(outputFilePath + patientID)


for patient in patientHlaMAP:
     hlas = patientHlaMAP[patient]
     for hla in hlas:
          os.makedirs(outputFilePath+ patientID +"/" + hla.replace(":","-") + "/")          
          for num in isomers:
               filename = "TCGA-" + patientID +"_Varscan_variants_filter.pass."+ str(num) +".peptide"
               functions.writeInputFile(filepath + "/peptides/", filename, patientID,str(num))
               
               proc = subprocess.Popen([sys.executable, 'epitope_new.py', hla, num , filename,  patientID, outputFolder, syfpeithiStr, IEDBStr, netmhcpanStr])
               procs.append(proc)
               
          



old_newHla_map = dict(zip(key_,value_))
for proc in procs:
     proc.wait()

for patient in patientHlaMAP:
     fwrite = open(outputFilePath+ patientID + "/"+patientID+ '.tsv','w')
     fwrite.write("PatientId"+tab+   "Allele"  + tab  + "Gene-transcript" + tab+ \
         "Mutant" +tab+    "Start_15mer"  +tab+  "Peptide_15mer" +tab+  "IEDB_Binding_15mer" +tab+ \
                           "Start_Position_17mer"  +tab+    "Peptide_17mer" +tab+ "IEDB_Binding_17mer"  +tab+ \
                           "Start_Position_19mer" +tab+    "Peptide_19mer" +tab+  "IEDB_Binding_19mer" +tab+  \
                           "Start_Position_21mer"  +tab+   "Peptide_21mer" +tab+ "IEDB_Binding_21mer" +tab+ \
                           "BestScoreMer_Mutant"  +tab+ \
          "WildType" +tab+ "Start_15mer"  +tab+  "Peptide_15mer"  +tab+ "IEDB_Binding_15mer"  +tab+ \
                         "Start_Position_17mer"   +tab+  "Peptide_17mer" +tab+ "IEDB_Binding_17mer" +tab+\
                          "Start_Position_19mer"   +tab+  "Peptide_19mer"+tab+  "IEDB_Binding_19mer"  +tab+ \
                          "Start_Position_21mer"    +tab+ "Peptide_21mer" +tab+ "IEDB_Binding_21mer" +tab+  "BestScoreMer_WildType"   +tab+ \
          "Same_Seq" +tab+ "Start_15mer"   +tab+ "Peptide_15mer" +tab+ "IEDB_Binding_15mer"  +tab+ \
                          "Start_Position_17mer"  +tab+   "Peptide_17mer" +tab+ "IEDB_Binding_17mer"  +tab+ \
                          "Start_Position_19mer"  +tab+   "Peptide_19mer" +tab+ "IEDB_Binding_19mer"  +tab+ \
                          "Start_Position_21mer"   +tab+  "Peptide_21mer" +tab+ "IEDB_Binding_21mer" +tab+ "BestScoreMer_SameSeq"  +tab + \
                           "isSameHLA\n") 
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
              
               #MTWTFile = functions.getFilePath(outputFilePath,"",patient) + "MTWT_File"  + peptideLenStr
               #transcript_SET, mutant_Map, wildType_Map = functions.initializeDataSets(MTWTFile)
               
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
               
               '''    
               #print hla + " ---> " + num
               #print peptideSet 
               IEDB_norm_map = functions.writeNormFile(filepath, outputFolder,"IEDB",peptideSet,IEDBMap,num)    
               bp_map = functions.writeTransformFile(filepath, outputFolder,IEDB_norm_map, peptideSet, num)
               ann_map = functions.getANNMap(outputFolder,num)
               print ann_map

               rFilePath = filepath + "/R/"

               commandR= "Rscript "+rFilePath+"/ANN_Immunogenicity.r " + outputFolder + " "  + rFilePath + " " + num
               returncode = subprocess.call(commandR,shell=True)

               while returncode != 0:
                    print "The ANN did not converge, running the ANN again for " + patient+ " - " + hla
                    print returncode
                    returncode = subprocess.call(commandR,shell=True)
               ann_prob = functions.getAnnProb(outputFolder, num)

               print ann_prob

               '''
               
          for t in transcript_SET_G:
               fwrite.write(patient+  tab + hla+ tab + t )
               
               merList  = list()
               mer =  namedtuple("mer",("num","score"))
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_transcriptMap_MT[t], num)
                    m = mer(num = num , score = bindingScore)                    
                    merList.append(m)
                    if num == "15":
                         fwrite.write(tab  + "MT") 
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))          
               fwrite.write(tab + functions.getFinalMer(merList))

               merList  = list()
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_transcriptMap_WT[t], num)
                    m = mer(num = num , score = bindingScore)                    
                    merList.append(m)
                    if num == "15":
                         fwrite.write(tab  + "WT") 
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
               fwrite.write(tab + functions.getFinalMer(merList))

               merList  = list()
               for num in isomers:
                    position, peptide, bindingScore = functions.getData(IEDB_TranscriptMap_SameSeq[t], num)
                    m = mer(num = num , score = bindingScore)                    
                    merList.append(m)
                    if num == "15":
                         fwrite.write(tab  + "Same_Seq")      
                    fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
               fwrite.write(tab + functions.getFinalMer(merList))

               newHla = functions.getNewHLA(old_newHla_map, hla)
               fwrite.write(tab + newHla+"\n")
                         
   
     fwrite.close()
'''
#############################################
# Create Sequence transcript map
#############################################
tab = "\t"`
epitope_score_map = dict()
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     fwrite = open(outputFilePath+ patientID + "/"+patientID+ '.tsv','w')
     fwrite.write("PatientId" + tab +	"Allele"+tab+"Position" +tab+	"Mutant"+tab+	"Gene-transcript"+tab+	"Peptide"+tab+	"IEDB_Binding_Scores"	+tab+"Start_Pos"+tab+	"End_Pos"+tab+	"Ann_Score"	+tab+"WT_Position"	+tab+ "WT"+tab+	"Gene-transcript"+tab+	"WT_peptide"+tab+	"WT_IEDB Binding_Scores"+tab+	"WT_Start_Pos" +tab+	"WT_End_Pos"+tab+	"DELTA-MT_WT" +tab+	"WT_Ann_Score" +tab+	"SameSeq_Position"	+tab+"Same_Seq"	+tab+"SameSeq_peptide"+tab+	"SameSeq_IEDB_Binding_Scores"+tab+	"SameSeq_Start_Pos" +tab+"SameSeq_End_Pos"	+tab+ "MT_WT(SameSeq)"+tab+	"SameSeq_Ann_Score\n")
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]
          IEDBStr = "False"
          netmhcpanStr = "False"
          isSameHla = "False"
          if hla in IEDBList:
               IEDBStr  = "True" 
          if hla in netmhcpanList:
                netmhcpanStr = "True"
          orignalHLA = ""
          morphedHLA = ""
          orignalHLA = hla
          if hla == old_newHla_map[hla]:
               isSameHla = "Yes"
          else:
               isSameHla = old_newHla_map[hla]
               isSameHla = isSameHla.replace(":","-")
               morphedHLA  = isSameHla
                   
          
          hla = hla.replace(":","-")
          filename = "Epitope_prob.txt"
          key_ =[]
          value_ =[]

#          print "inside if else " + hla + "  "+ patientID
          if isSameHla == "Yes":
               epitope_file_path  = outputFilePath + patientID + "/"+  hla + "/" + filename
          else:
               epitope_file_path  = outputFilePath + patientID + "/"+  morphedHLA + "/" + filename
#          print epitope_file_path      
          try:
               with open(epitope_file_path ) as f:
                    for line in f:
                         #fwrite.write(line)
                         if isSameHla == "Yes":
                              fwrite.write(line.strip() +tab+ "" + "\n")
                         else:
                              fwrite.write(line.strip() +tab+ hla + "\n")
          except IOError as e:
               print "I/O error({0}): {1}".format(e.errno, e.strerror)                            
     fwrite.close()




#############################################
# End of Program
#############################################     


'''
