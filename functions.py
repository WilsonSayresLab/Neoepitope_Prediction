import subprocess
import os
import sys
import operator
import ast
from subprocess import call
from bs4 import BeautifulSoup
from StringIO import StringIO
import pycurl
import re
from datetime import datetime
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter


###################################################
# python epitope.py H-2-Db 9
###################################################


###################################################

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

def getAnnProb(outputFilePath,num):
    key_ = []
    value_ = []
    counter = 0
    with open(outputFilePath+ "Ann_output."+num+".txt") as f:
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
                        val=0
                        if "netmhcpan" in tempLine[6]: 
                          val = float(tempLine[14].strip())
                          value_.append(val)
                        else:
                          val = float(tempLine[8].strip())
                          value_.append(val)
                        iedb = IEDB(peptide = peptideStr, start = tempLine[2].strip(), end = tempLine[3].strip(), score = val)
                        IEDB_TupleList.append(iedb)
    map_values = dict(zip(key_,value_))
    ALL_LISTS.append(map_values)
    if("output_IEDB" in fileName ):
      ALL_LISTS.append(IEDB_TupleList)
    if("output_netmhcpan" in fileName ):
      ALL_LISTS.append(NetMHC_TupleList)
    return ALL_LISTS
    
###################################################

def writeNormFile(rFilePath, filepath,fileName,final_set,final_map,num):
    result = []
    key_ = []
    value_ = []          
    fwrite = open(filepath + fileName + '_norm.'+num+'.csv','w')
    fwrite.write("Sequence,"+ fileName+ ".bind\n")
    for val in final_set:
         fwrite.write(str(val) + ","+ str(final_map[val]) + "\n")
    fwrite.close()
    rFilePath = rFilePath + "/R/"
    commandR = "Rscript "+ rFilePath +fileName + "_norm.R " + filepath + " "+ num
    returncode = subprocess.call(commandR,shell=True)
    counter= 0
    
    with open(filepath + fileName+"_normalized."+num+".txt") as f:
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

def getANNMap(outputFilePath, num):
    result = []
    key_ = []
    value_ = []
    templine = []
    counter = 0
    fwrite = open(outputFilePath + "ANN_input_final."+num+".txt",'w')
    with open(outputFilePath + "Ann_input."+num+".txt") as f:
      for line in f:
          templine = line.strip().split('\t')
          if counter == 0:
            counter= counter + 1
            if (int(num) >14):
              fwrite.write(templine[2].strip()+"\t"+templine[3].strip()+"\t"+ templine[4].strip()+"\t"+templine[5].strip()+"\t"+templine[6].strip()+"\t"+templine[7].strip()+"\t"+templine[8].strip()+"\t"+templine[9].strip())
            if (int(num) >16):   
              fwrite.write("\t"+templine[10].strip())
            if (int(num) >18):   
              fwrite.write("\t"+templine[11].strip())
            if (int(num) >20):   
              fwrite.write("\t"+templine[12].strip())
            fwrite.write("\n")
          else:
            counter = counter+1
            key_.append(templine[1].strip().replace("\"",""))
            value_.append(templine[0].strip().replace("\"",""))

            if (int(num) >14): 
              fwrite.write(templine[2].strip()+"\t"+templine[3].strip()+"\t"+templine[4].strip()+"\t"+templine[5].strip()+"\t"+templine[6].strip()+"\t"+templine[7].strip()+"\t"+templine[8].strip()+"\t"+templine[9].strip())
            if (int(num) >16):   
              fwrite.write("\t"+templine[10].strip())
            if (int(num) >18):   
              fwrite.write("\t"+templine[11].strip())
            if (int(num) >20):   
              fwrite.write("\t"+templine[12].strip())
            fwrite.write("\n")
    fwrite.close()
    map_values = dict(zip(key_,value_))
    return map_values

###################################################
###################################################

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

def initializeDataSets(filepath):
  ALL_LISTS = []
  transcript_SET = set()
  mutant_Map = defaultdict(list)
  wildType_Map = defaultdict(list)
  transcript = ""
  counter = 0
  with open(filepath ) as f:
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

def getLowestScore(peptideList, peptideTuples):
  peptideTuples = sorted(peptideTuples, key=attrgetter('score'))
  
  for peptideTuple in peptideTuples:
    if peptideTuple.peptide in peptideList:
      return peptideTuple
      break

#############################################
# Create Sequence transcript map
#############################################  

def getFilePath(outputFilePath, hla,patientID):
  ## Some Constants
  hla_allele = hla.replace(":","-")
  return outputFilePath +patientID + "/" + hla_allele +  "/"

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

def writeInputFile(filepath,file, patientID, num):
     seq=""
     header=""
     data=""
     filename = filepath+patientID + "." + num+ '.txt'
     peptideFileName = filepath +num+"mer/" + file 
     fwrite = open(  filename,'w')
     with open(peptideFileName) as f:
          for line in f:
               if ">" in line:
                    header = line.strip();
                    data=""
               else:
                    data = line.strip().replace("X","")
               n = int(num) * 2
               if len(data) > (int(num) / 2):
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

def getMutantWildTypeData(filepath, patientId, num):
#############################################
# Create Sequence transcript map
#############################################

  ALL_LISTS = []
  seq_transcript_map = defaultdict(list)
  transcript_SET = set()

  filename =  patientId +"."+str(num)+ ".txt"
  transcript = ""
  trans_ = ""
  counter = 0
  mutant_Map = defaultdict(list)
  wildType_Map = defaultdict(list)
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
           length = (int(num)/2) + 1 
           seq_transcript_map[seq].append(transcript)
           if "MT." in transcript:
                mutant_Map[trans_].append(getAllSeq(seq,length))
                #s = mutant_Map[trans_]
                #fwriteMT.write("MT\t" + str(trans_) + "\t" + str(s[0]) + "\n")
           if "WT." in transcript:
                wildType_Map[trans_].append(getAllSeq(seq,length))
                #s = wildType_Map[trans_]
                #fwriteMT.write("WT\t" + str(trans_) + "\t" + str(s[0]) + "\n")       
  ALL_LISTS.append(transcript_SET)
  ALL_LISTS.append(mutant_Map)
  ALL_LISTS.append(wildType_Map)
  return ALL_LISTS

def getSameSeqScore(mutantLowestScoreTuple,wildTypePeptideList, peptideTuples):

  peptide_tuple_sameSeq = [item for item in peptideTuples if str(item.start) == str(mutantLowestScoreTuple.start)]

  for tup in peptide_tuple_sameSeq:
    if tup.peptide in wildTypePeptideList:
      return tup

def getPeptides(IEDB_transcriptMap, num):
  peptideSet = set()
  for item in IEDB_transcriptMap:
    if item.mer == num:
      peptideSet.add( item.data.peptide)
  return peptideSet

def writeTransformFile(filepath, outputFilePath, norm_map, peptideSet, num):

  bp_key_ = []
  bp_value_ = []
  fwrite = open(outputFilePath+'transform_input.'+num+'.csv','w')
  fwrite.write("Epitope,IEDB.Norm,BP Score\n")
  final_set = set(peptideSet)
  for val in final_set:
    val = val.replace("'","")
    bp_score = norm_map[val];
    bp_key_.append(val)
    bp_value_.append(bp_score)
    fwrite.write(val + "," + str(norm_map[val])  + ","  + str(bp_score) +"\n" )

  fwrite.close()
  bp_map = dict(zip(bp_key_,bp_value_))
  rFilePath = filepath + "/R/"
  commandR= "Rscript "+rFilePath+"/Peptidematrix.R " + outputFilePath + " " + rFilePath + " " + num
  returncode = subprocess.call(commandR,shell=True)
  return bp_map
    

def getData(dataMapArr, num):
  dataTuples = namedtuple("dataTuple", ("mer","data"))
  dataTuples = dataMapArr
  t = list()
  for s in dataTuples:
    if s.mer == num:
      dataField = namedtuple("IEDB", ("peptide","start","end","score"))
      dataField = s.data
      t.append(str(dataField.start))
      t.append(str(dataField.peptide )) 
      t.append(str(dataField.score))
  return tuple(t)

def getNewHLA(old_newHla_map, hla):
  if hla == old_newHla_map[hla]:
    isSameHla = "Yes"
    return ""
  else:
    isSameHla = old_newHla_map[hla]
    isSameHla = isSameHla.replace(":","-")
    morphedHLA  = isSameHla
    return morphedHLA
    
def getFinalMer(merList): 
  mer = sorted(merList, key=attrgetter('score'))
  return mer[0].num


