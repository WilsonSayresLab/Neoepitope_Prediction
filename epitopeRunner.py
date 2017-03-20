from os import listdir
from os.path import isfile, join, dirname
import os
import sys
import subprocess
from datetime import datetime
from collections import defaultdict


def readSyphethi():
     syphethi = []
     with open("syphethi.txt" ) as f:
          for line in f:
               syphethi.append(line.strip())
     return syphethi

def readIEDB():
     IEDB = []
     with open("IEDB.txt" ) as f:
          for line in f:
               IEDB.append(line.strip())
     return IEDB

def readnetmhcpan():
     netmhcpan = []
     with open("netmhcpan.txt" ) as f:
          for line in f:
               netmhcpan.append(line.strip())
     return netmhcpan

def writeInputFile(filepath,file, patientID):
     seq=""
     header=""
     data=""
     fwrite = open(  filepath+patientID + '.txt','w')
     with open(filepath + file ) as f:
          for line in f:
               if ">" in line:
                    header = line.strip();
                    data=""
               else:
                    data = line.strip().replace("X","")
                    
               if len(data) > 20:
                    #print header + "\n" + data+"\n"
                    fwrite.write(header + "\n" + data+"\n")

     fwrite.close();

def getClosestHLA(hla,strList):
     minDiff = 1000000
     minDiffHLa = ""
     hlaLists =  [s for s in strList if hla[:6] in s]  
     for h in hlaLists:
          diff = abs(float(hla[6:].replace(":","")) - float(h[6:].replace(":",""))) 
          if(diff  <= minDiff):
               minDiff = diff
               minDiffHLa = h
     return minDiffHLa

def getAllSeq(seq, length):
     l = len(seq) + 1
     list_ = []
     for i in range(l-length):
          list_.append(seq[ : length])
          seq = seq[1:]
     return list_



## Get the current directory
filepath = os.path.dirname(os.path.realpath(__file__)) + "/"
hlapath = filepath+ "hla"
allfiles = [f for f in listdir(hlapath) if isfile(join(hlapath, f)) and f.endswith(".txt")]
allfiles = set(allfiles)
fwrite = open(  'errorHLA.txt','a')
length = 11
#print allfiles

##############################
##############################
syfpeithiList= set(readSyphethi())
IEDBList = set(readIEDB())
netmhcpanList = set(readnetmhcpan())

file_hlas_map =dict();
for file in allfiles:
     file_hlas_map[file] = list();
     with open(hlapath+ "/"+ file) as f:
          for line in f:
               hlas = line.strip().split("\t")
               file_hlas_map[file].extend(list(set(hlas[1:len(hlas)])))

outputFolder =  str(datetime.utcnow()).replace(" ","-").replace(":","-").replace(".","-") + "/"
outputFilePath = filepath + "output/" + outputFolder
#print file_hlas_map



#############################################
# Create Sequence transcript map
#############################################

seq_transcript_map = defaultdict(list)
transcript_SET = set()
mutant_Map = defaultdict(list)
wildType_Map = defaultdict(list)

for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     filename = "VarScan2-TCGA-" + patientID + "-TNBC.peptide"
     #filename = "VarScan2_TCGA-" + patientID + "-nonTNBC.peptide"
     writeInputFile(filepath + "/peptides/", filename, patientID);
     filename =  patientID + ".txt"
     transcript = ""
     trans_ = ""
     counter = 0
     if not os.path.exists(outputFilePath + patientID ):
          os.makedirs(outputFilePath + patientID) 
     fwriteMT = open(outputFilePath+ patientID + '/MTWT_File.txt','w')
     with open(filepath + "peptides/" + filename ) as f:
          for line in f:
               seq=""
               if ">" in line:
                    transcript = line.strip().replace(">","")
                    trans_ = transcript.replace("MT.","").replace("WT.","")
                    transcript_SET.add(trans_)
                    counter += 1
               else:
                    seq = line.strip()
                    counter += 1
               if counter % 2 ==0:
                    seq_transcript_map[seq].append(transcript)
                    if "MT." in transcript:
                         mutant_Map[trans_].append(getAllSeq(seq,length))
                         s = mutant_Map[trans_]
                         fwriteMT.write("MT\t" + str(trans_) + "\t" + str(s[0]) + "\n")
                    if "WT." in transcript:
                         wildType_Map[trans_].append(getAllSeq(seq,length))   
                         s = wildType_Map[trans_]
                         fwriteMT.write("WT\t" + str(trans_) + "\t" + str(s[0]) + "\n")       
     fwriteMT.close()


#############################################
# Core Processing
#############################################


procs = []
key_ = []
value_ = []

for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     #filename = "VarScan2_TCGA-" + patientID + "-nonTNBC.peptide"
     filename = "VarScan2-TCGA-" + patientID + "-TNBC.peptide"
    # writeInputFile(filepath + "/peptides/", filename, patientID);
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]
          if not os.path.exists(outputFilePath + patientID):
               os.makedirs(outputFilePath + patientID)
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
               hla =  getClosestHLA(hla, IEDBList)
               value_.append(hla)
          else:
               value_.append(hla)

          proc = subprocess.Popen([sys.executable, 'epitope_A1.py', hla, '11', filename,  patientID, outputFolder, syfpeithiStr, IEDBStr, netmhcpanStr])
          procs.append(proc)
          
fwrite.close()

old_newHla_map = dict(zip(key_,value_))
for proc in procs:
     proc.wait()


#############################################
# Create Sequence transcript map
#############################################
tab = "\t"
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
