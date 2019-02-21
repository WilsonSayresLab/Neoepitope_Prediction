import subprocess
import os
import sys
import operator
from subprocess import call
from bs4 import BeautifulSoup
from StringIO import StringIO
import pycurl
import re
from datetime import datetime
from collections import defaultdict
from collections import namedtuple

###################################################
# python epitope.py H-2-Db 9
###################################################

###################################################
class epitopes:
###################################################
###################################################

###################################################
     @staticmethod
     def getTop200(list_ ):
          counter=0
          result = []
          for val in list_:
                     if counter == 2000000:
                             break
                     result.append(val[0])
                     counter = counter+1
          return result
###################################################
     @staticmethod
     def getAnnProb(outputFilePath):
          key_ = []
          value_ = []
          counter = 0
          with open(outputFilePath+ "Ann_output.txt") as f:
               for line in f:
                    if counter == 0:
                         counter=counter+1
                         continue
                    else:
                         counter = counter +1
                         templine = line.split("\t")
                         key_.append(templine[0].strip().replace("\"",""))
                         value_.append(float(templine[1].strip().replace("\"","")))
          map_values = dict(zip(key_,value_))
          return map_values


     @staticmethod
     def getMapwithValues(fileName):
          key_ = []
          value_ = []
          counter = 0
          ALL_LISTS = []
          if("output_netmhcpan.txt" in fileName ):
            NetMHC = namedtuple("NetMHC", ("peptide","start","end","score"))
            NetMHC_TupleList = []
            NetMHC_Set= set()
          if("output_IEDB.txt" in fileName ):
            IEDB = namedtuple("IEDB", ("peptide","start","end","score"))
            IEDB_TupleList = []
            IEDB_Set = set()
          with open(fileName) as f:
                  for line in f:
                          if counter == 0:
                                  counter=counter+1
                                  continue
                          else:
                                  counter=counter+1
                                  tempLine = line.split('\t')
                          if("syfpeithi.txt" in  fileName):
                                       key_.append(tempLine[0].strip())
                                       value_.append(float(tempLine[1].strip()))
                          if("output_netmhcpan.txt" in fileName ):
                            peptideStr = tempLine[5].strip()
                            if  peptideStr not in NetMHC_Set:
                              NetMHC_Set.add(peptideStr)
                              key_.append(peptideStr)
                              value_.append(float(tempLine[6].strip()))
                              netmhc = NetMHC(peptide = peptideStr, start = tempLine[2].strip(), end = tempLine[3].strip(), score = tempLine[6].strip())
                              NetMHC_TupleList.append(netmhc)
                          if("output_IEDB.txt" in fileName):
                            peptideStr = tempLine[5].strip()
                            if  peptideStr not in IEDB_Set:
                              IEDB_Set.add(peptideStr)
                              key_.append(peptideStr)
                              value_.append(float(tempLine[6].strip()))
                              iedb = IEDB(peptide = peptideStr, start = tempLine[2].strip(), end = tempLine[3].strip(), score = tempLine[6].strip())
                              IEDB_TupleList.append(iedb)
          map_values = dict(zip(key_,value_))
          ALL_LISTS.append(map_values)
          if("output_IEDB.txt" in fileName ):
            ALL_LISTS.append(IEDB_TupleList)
          if("output_netmhcpan.txt" in fileName ):
            ALL_LISTS.append(NetMHC_TupleList)
          return ALL_LISTS


###################################################
     @staticmethod
     def getMapwithValuesIEDB(fileName):
          key_ = []
          value_ = []
          counter = 0
          ALL_LISTS = []
          if("output_netmhcpan" in fileName ):
            NetMHC = namedtuple("NetMHC", ("peptide","start","end","score"))
            NetMHC_TupleList = []
            NetMHC_Set= set()
          if("output_IEDB" in fileName ):
            print "in True"
            IEDB = namedtuple("IEDB", ("peptide","start","end","score"))
            IEDB_TupleList = []
            IEDB_Set = set()
          print "printing filename"
          print fileName
          with open(fileName) as f:
                  for line in f:
                          if counter == 0:
                                  counter=counter+1
                                  continue
                          else:
                                  counter=counter+1
                                  tempLine = line.split('\t')
                          if("syfpeithi.txt" in  fileName):
                                       key_.append(tempLine[0].strip())
                                       value_.append(float(tempLine[1].strip()))
                          if("output_netmhcpan" in fileName ):
                            peptideStr = tempLine[5].strip()
                            if  peptideStr not in NetMHC_Set:
                              NetMHC_Set.add(peptideStr)
                              key_.append(peptideStr)
                              value_.append(float(tempLine[6].strip()))
                              netmhc = NetMHC(peptide = peptideStr, start = tempLine[2].strip(), end = tempLine[3].strip(), score = tempLine[6].strip())
                              NetMHC_TupleList.append(netmhc)
                          if("output_IEDB" in fileName):
                            peptideStr = tempLine[5].strip()
                            if  peptideStr not in IEDB_Set :
                              IEDB_Set.add(peptideStr)
                              key_.append(peptideStr)
                              value_.append(float(tempLine[14].strip()))
                              iedb = IEDB(peptide = peptideStr, start = tempLine[2].strip(), end = tempLine[3].strip(), score = tempLine[14].strip())
                              IEDB_TupleList.append(iedb)
          map_values = dict(zip(key_,value_))
          ALL_LISTS.append(map_values)
          if("output_IEDB" in fileName ):
            ALL_LISTS.append(IEDB_TupleList)
          if("output_netmhcpan" in fileName ):
            ALL_LISTS.append(NetMHC_TupleList)
          return ALL_LISTS

###################################################
     @staticmethod
     def writeNormFile(rFilePath, filepath,fileName,final_set,final_map):
          result = []
          key_ = []
          value_ = []
          fwrite = open(filepath + fileName + '_norm.csv','w')
          fwrite.write("Sequence,"+ fileName+ ".bind\n")
          for val in final_set:
               fwrite.write(str(val) + ","+ str(final_map[val]) + "\n")
          fwrite.close()
          commandR = "Rscript "+ rFilePath +fileName + "_norm.R " + filepath
          returncode = subprocess.call(commandR,shell=True)
          counter= 0
          with open(filepath + fileName+"_normalized.txt") as f:
               for line in f:
                          if counter == 0:
                                  counter=counter+1
                                  continue
                          else:
                                  counter=counter+1
                                  tempLine = line.split('\t')
                                  key_.append(tempLine[1].strip().replace("'","").replace("\"",""))
                                  value_.append(float(tempLine[3].strip()))
          map_values = dict(zip(key_,value_))
          return map_values
###################################################
     @staticmethod
     def writeCombinedProb(final_output,bp_map,ann_map,ann_prob, outputFilePath, tuples, teanscriptMap_MT, teanscriptMap_WT, transcript_SET):
          fwrite = open(outputFilePath +'Epitope_prob.txt','w')

          for t in transcript_SET:
            value = teanscriptMap_MT[t]
            prob = bp_map[value] * (1- ann_prob[ann_map[value]])
            peptide_tuple = [item for item in tuples if item.peptide == value]
            fwrite.write("MT\t" + t + "\t" + value + "\t" + str(prob)+"\t" + peptide_tuple[0].start + "\t" + peptide_tuple[0].end  + "\t"+ peptide_tuple[0].score + "\n" )
            value_WT = teanscriptMap_WT[t]
            prob = bp_map[value_WT] * (1- ann_prob[ann_map[value_WT]])
            peptide_tuple = [item for item in tuples if item.peptide == value_WT]
            fwrite.write("WT\t" + t + "\t" + value_WT + "\t" + str(prob)+"\t" + peptide_tuple[0].start + "\t" + peptide_tuple[0].end  + "\t" + peptide_tuple[0].score + "\t" + (bp_map[value] - bp_map[value_WT]) + "\n" )
          fwrite.close()

###################################################
     @staticmethod
     def getANNMap(outputFilePath):
          result = []
          key_ = []
          value_ = []
          templine = []
          counter = 0
          fwrite = open(outputFilePath + "ANN_input_final.txt",'w')
          with open(outputFilePath + "Ann_input.txt") as f:
            for line in f:
                templine = line.strip().split('\t')
                if counter == 0:
                  counter= counter + 1
                  fwrite.write(templine[2].strip()+"\t"+templine[3].strip()+"\t"+ templine[4].strip()+"\t"+templine[5].strip()+"\t"+templine[6].strip()+"\t"+templine[7].strip()+"\t"+templine[8].strip()+"\t"+templine[9].strip()+"\t"+templine[10].strip()+ "\n")
                else:
                  counter = counter+1
                  key_.append(templine[1].strip().replace("\"",""))
                  value_.append(templine[0].strip().replace("\"",""))
                  fwrite.write(templine[2].strip()+"\t"+templine[3].strip()+"\t"+templine[4].strip()+"\t"+templine[5].strip()+"\t"+templine[6].strip()+"\t"+templine[7].strip()+"\t"+templine[8].strip()+"\t"+templine[9].strip()+"\t"+templine[10].strip()+"\n")
          fwrite.close()
          map_values = dict(zip(key_,value_))
          return map_values

###################################################
###################################################
     @staticmethod
     def parseSypethi(motif,amers,seq):
          motif = motif.replace("H-2","H2")
          url= "http://www.syfpeithi.de/bin/MHCServer.dll/EpitopePrediction?Motif="+str(motif)+"&amers="+str(amers)+"&SEQU="+str(seq)+"&DoIT=++Run++"
          url = str( url)
          storage = StringIO()
          c = pycurl.Curl()
          c.setopt(c.URL, url)
          c.setopt(c.WRITEFUNCTION, storage.write)
          c.perform()
          c.close()
          content = storage.getvalue()
          soup = BeautifulSoup(content,'html.parser')
          tables = soup.find_all("table")
          #temprows = unicode(rows, "utf-8",errors="ignore")
          if(len(tables) >1):
            rows = tables[1].find_all("tr")
            key_ =[]
            val_ =[]
            for row in rows:
                    tds= row.find_all("td")
                    counter = 0
                    for td in tds:
                            if counter == 0:
                                    counter = counter + 1
                                    continue
                            if counter == 1:
                                    counter = counter + 1
                                    tempvar = re.sub(r'[^a-zA-Z]', "",td.getText())
                                    if tempvar == "gototop":
                                            continue
                                    key_.append(tempvar.strip())
                                    continue
                            if counter == 2:
                                    counter = counter + 1
                                    tempnum = td.getText().strip()
                                    val_.append(float(tempnum))
                                    continue
            map_values = dict(zip(key_,val_))
            return map_values
          else:
            return -1


###################################################
     @staticmethod
     def getSeq(file):
          seqList=[]
          seq=""
          with open(file) as f:
               for line in f:
                          if ">" in line:
                                  continue
                          else:
                            seqList.append(str(line.strip()))
          seqSet = set(seqList)
          for s in seqSet:
            seq += str(s.strip())
          return seq


#############################################
# Create Sequence transcript map
#############################################
     @staticmethod
     def initializeDataSets(filepath,fileName):
      ALL_LISTS = []
      transcript_SET = set()
      mutant_Map = defaultdict(list)
      wildType_Map = defaultdict(list)
      transcript = ""
      counter = 0
      print "MT-WT "
      print fileName
      with open(filepath + fileName ) as f:
        for line in f:
          data = line.split("\t")
          if data[0] == "MT":
            transcript_SET.add(data[1].strip())
            mutant_Map[data[1]] = data[2].strip()
          else:
            transcript_SET.add(data[1].strip())
            wildType_Map[data[1]] = data[2].strip()
      ALL_LISTS.append(transcript_SET)
      ALL_LISTS.append(mutant_Map)
      ALL_LISTS.append(wildType_Map)
      return ALL_LISTS
#############################################
# Create Sequence transcript map
#############################################
     @staticmethod
     def getLowestScore(peptideList,peptide_score):
      key_ =[]
      val_ =[]
      peptideList = peptideList.replace("'","").replace("[","").replace("]","").replace(" ","").replace("'","")
      peptideList = peptideList.split(",")
      for p in peptideList:
        key_.append(p)
        val_.append(peptide_score[p])
      map_scores = dict(zip(key_,val_))
      map_scores = sorted(map_scores.iteritems(),key=operator.itemgetter(1))
      return map_scores[0]




###################################################
argument = sys.argv
## argument[1] = Allele
## argument[2] = size of peptide
## argument[3] = inputFileName
## argument[4] = patientId
## argument[5] = timestamp of folder
## argument[6] = syfpeithi allele boolen flag

ep = epitopes()

filepath = os.path.dirname(os.path.realpath(__file__))
## Some Constants
hla_allele = argument[1].replace(":","-")
outputFilePath = filepath + "/output/" + argument[5] + argument[4] + "/" + hla_allele +  "/"



seqLen = argument[2]
peptideLen = (int(seqLen) / 2) + 1
peptideLenStr = "."+ seqLen + ".txt"

inputFilePath = filepath + "/peptides/"
outputNetmhcpanFile = outputFilePath + "output_netmhcpan"  + peptideLenStr
outputIEDBFile= outputFilePath + "output_IEDB" + peptideLenStr
rFilePath = filepath + "/R/"
inputFile = inputFilePath + argument[4] + "."+ seqLen+ ".txt"
patientId = argument[4]
tab = "\t"


command = "python "+filepath+"/mhc_i/src/predict_binding.py netmhcpan " + argument[1] +" "+ argument[2] +" "+ inputFile + " > "  +   outputNetmhcpanFile
command1 = "python "+filepath+"/mhc_i/src/predict_binding.py IEDB_recommended " + argument[1] +" "+ str(peptideLen) +" "+ inputFile + " > " + outputIEDBFile



#returncode = subprocess.call(command,shell=True)

returncode = subprocess.call(command1,shell=True)
method = "IEDB"

if os.path.getsize(outputIEDBFile) ==0:
  method= "NETMHC"
  commandnetMhc = "python "+filepath+"/mhc_i/src/predict_binding.py netmhcpan " + argument[1] +" "+ str(peptideLen) +" "+ inputFile + " > " + outputIEDBFile
  returncode = subprocess.call(commandnetMhc,shell=True)

'''
###################################################
# get the Map of key and values
###################################################
print "begining of IEDB"
IEDB_top200= []
IEDB_map = defaultdict(list)
IEDB_tuples = defaultdict(list)
if method == "IEDB":
  IEDB_data = ep.getMapwithValuesIEDB(outputIEDBFile)
else:
  IEDB_data = ep.getMapwithValues(outputIEDBFile)
IEDB_map[seqLen] = IEDB_data[0]
IEDB_tuples[seqLen] =  IEDB_data[1]
#sorted_IEDB = sorted(IEDB_map.iteritems(),key=operator.itemgetter(1))
print "printing IEDB data"
print IEDB_data
print "printing IEDB MAP data"
print IEDB_map
print "printing IEDB TUPLES data"
print IEDB_tuples


###################################################
# Transcript - peptide - min score
###################################################
IEDB_transcriptMap_MT = defaultdict(list)
IEDB_transcriptMap_WT = defaultdict(list)

sameSeqTupleList= list()
for t in transcript_SET:
  IEDB_transcriptMap_MT[t][seqLen] = ep.getLowestScore(mutant_Map[t],IEDB_map[seqLen])
  IEDB_transcriptMap_WT[t][seqLen] = ep.getLowestScore(wildType_Map[t],IEDB_map[seqLen])
  print IEDB_transcriptMap_MT
  print IEDB_transcriptMap_MT
  value = IEDB_transcriptMap_MT[t]
  value = value[0]
  peptide_tuple_MT = [item for item in IEDB_tuples if item.peptide == value]
  mutantStart = peptide_tuple_MT[0].start
  peptide_tuple_sameSeq = [item for item in IEDB_tuples if str(item.start) == str(mutantStart)]
  peptideList = wildType_Map[t]
  peptideList = peptideList.replace("'","").replace("[","").replace("]","").replace(" ","").replace("'","")
  peptideList = peptideList.split(",")
  for tup in peptide_tuple_sameSeq:
    if tup.peptide in peptideList:
      sameSeqTupleList.append(tup.peptide)

  #IEDB_top200 = ep.getTop200(sorted_IEDB)
IEDB_top200 = list(value[0] for key , value in IEDB_transcriptMap_WT.iteritems())
IEDB_top200 = list(set(IEDB_top200 + ( list(value[0] for key , value in IEDB_transcriptMap_MT.iteritems()))))
IEDB_top200 = IEDB_top200 +  sameSeqTupleList

IEDB_norm_map = ep.writeNormFile(rFilePath, outputFilePath,"IEDB",IEDB_top200,IEDB_map)

print "END of IEDB"


###################################################
### Final Set with top 200
###################################################

final_set=[]
bp_key_ = []
bp_value_ = []

fwrite = open(outputFilePath+'transform_input.csv','w')
fwrite.write("Epitope,IEDB.Norm,BP Score\n")

final_set = set(IEDB_top200)
for val in final_set:
  val = val.replace("'","")
  bp_score = IEDB_norm_map[val];
  bp_key_.append(val)
  bp_value_.append(bp_score)
  fwrite.write(val + "," + str(IEDB_norm_map[val])  + ","  + str(bp_score) +"\n" )

fwrite.close()
bp_map = dict(zip(bp_key_,bp_value_))

###################################################
commandR= "Rscript "+rFilePath+"/Peptidematrix.R " + outputFilePath + " " + rFilePath
returncode = subprocess.call(commandR,shell=True)
ann_map = ep.getANNMap(outputFilePath)
###################################################

commandR= "Rscript "+rFilePath+"/ANN_Immunogenicity.r " + outputFilePath + " "  + rFilePath
returncode = subprocess.call(commandR,shell=True)

while returncode != 0:
     print "The ANN did not converge, running the ANN again for " + argument[4]+ " - " + argument[1]
     print returncode
     commandR= "Rscript "+rFilePath+"/ANN_Immunogenicity.r " + outputFilePath + " " + rFilePath
     returncode = subprocess.call(commandR,shell=True)
ann_prob = ep.getAnnProb(outputFilePath)


fwrite = open(outputFilePath +'Epitope_prob.txt','w')
for t in transcript_SET:
  value = IEDB_transcriptMap_MT[t]
  value = value[0]
  prob = bp_map[value] * (1- ann_prob[ann_map[value]])
  peptide_tuple_MT = [item for item in IEDB_tuples if item.peptide == value]
  mutantStart = peptide_tuple_MT[0].start
  fwrite.write(patientId+  tab + hla_allele+ tab +  peptide_tuple_MT[0].start + tab + "MT" + tab+  t + tab+ value + tab + peptide_tuple_MT[0].score + tab  + peptide_tuple_MT[0].start + tab + peptide_tuple_MT[0].end  + tab  +str(prob) + tab  )
  value_WT = IEDB_transcriptMap_WT[t]
  value_WT = value_WT[0]
  prob = bp_map[value_WT] * (1- ann_prob[ann_map[value_WT]])
  peptide_tuple_WT = [item for item in IEDB_tuples if item.peptide == value_WT]
  fwrite.write(peptide_tuple_WT[0].start + tab + "WT" + tab  + t + tab + value_WT + tab + peptide_tuple_WT[0].score + tab + peptide_tuple_WT[0].start  + tab + peptide_tuple_WT[0].end + tab + str(float(peptide_tuple_MT[0].score) - float(peptide_tuple_WT[0].score)) + tab+ str(prob)+ tab )
  peptide_tuple_sameSeq = [item for item in IEDB_tuples if str(item.start) == str(mutantStart)]
  peptideList = wildType_Map[t]
  peptideList = peptideList.replace("'","").replace("[","").replace("]","").replace(" ","").replace("'","")
  peptideList = peptideList.split(",")
  for tup in peptide_tuple_sameSeq:
    if tup.peptide in peptideList:
      sameSeqTuple = tup
  prob = bp_map[sameSeqTuple.peptide] * (1- ann_prob[ann_map[sameSeqTuple.peptide]])
  fwrite.write(sameSeqTuple.start + tab + "SameSeq" + tab  + sameSeqTuple.peptide + tab + sameSeqTuple.score + tab + sameSeqTuple.start  + tab + sameSeqTuple.end + tab + str(float(peptide_tuple_MT[0].score) - float(sameSeqTuple.score)) + tab+ str(prob)+  "\n" )

fwrite.close()


###################################################
# End of Program
###################################################
'''
