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
                     if counter == 200:
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

###################################################
     @staticmethod
     def getMapwithValues(fileName):
          key_ = []
          value_ = []
          counter = 0
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
                                       key_.append(tempLine[5].strip())
                                       value_.append(float(tempLine[6].strip()))
                          if("output_IEDB.txt" in fileName):
                                       key_.append(tempLine[5].strip())
                                       value_.append(float(tempLine[8].strip()))
          map_values = dict(zip(key_,value_))
          return map_values
          
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
     def writeCombinedProb(final_output,bp_map,ann_map,ann_prob, outputFilePath):
          fwrite = open(outputFilePath +'Epitope_prob.txt','w')
          for value in final_output:
               prob = bp_map[value] * (1- ann_prob[ann_map[value]])
               fwrite.write(value + "\t" + str(prob)+"\n")
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
          seq=""
          with open(file) as f:
               for line in f:
                          if ">" in line:
                                  continue
                          else:
                            seq += str(line.strip())
          return seq
          

###################################################
argument = sys.argv
## argument[1] = Allele
## argument[2] = size of peptide
## argument[3] = inputFileName
## argument[4] = patientId
## argument[5] = timestamp of folder
ep = epitopes()

filepath = os.path.dirname(os.path.realpath(__file__))
## Some Constants

outputFilePath = filepath + "/output/" + argument[5] + argument[4] + "/" + argument[1].replace(":","-") +  "/"

if not os.path.exists(outputFilePath):
	os.makedirs(outputFilePath)

inputFilePath = filepath + "/peptides/"
outputNetmhcpanFile = outputFilePath + "output_netmhcpan.txt"
outputIEDBFile = outputFilePath + "output_IEDB.txt"
rFilePath = filepath + "/R/"
inputFile = inputFilePath + argument[3]
#rOutputFilePath = rFilePath + argument[4] + "/"
#if not os.path.exists(rOutputFilePath):
#  os.makedirs(ro)

command = "python "+filepath+"/mhc_i/src/predict_binding.py netmhcpan " + argument[1] +" "+ argument[2] +" "+ inputFile + " > "  +   outputNetmhcpanFile
command1 = "python "+filepath+"/mhc_i/src/predict_binding.py IEDB_recommended " + argument[1] +" "+ argument[2] +" "+ inputFile + " > " + outputIEDBFile

returncode = subprocess.call(command,shell=True)
returncode = subprocess.call(command1,shell=True)

###################################################
# get the Map of key and values
###################################################

netmhcpan_top200=[]
netmhcpan_map = ep.getMapwithValues(outputNetmhcpanFile)
sorted_netmhcpan = sorted(netmhcpan_map.iteritems(),key=operator.itemgetter(1))
netmhcpan_top200 = ep.getTop200(sorted_netmhcpan)
netmhcpan_norm_map = ep.writeNormFile(rFilePath ,outputFilePath,"NetMHC",netmhcpan_top200,netmhcpan_map)

###################################################
# get the Map of key and values
###################################################

IEDB_top200= []
IEDB_map = ep.getMapwithValues(outputIEDBFile)
sorted_IEDB = sorted(IEDB_map.iteritems(),key=operator.itemgetter(1))
IEDB_top200 = ep.getTop200(sorted_IEDB)
IEDB_norm_map = ep.writeNormFile(rFilePath, outputFilePath,"IEDB",IEDB_top200,IEDB_map)

###################################################
# get the Map of key and values
###################################################
syfpeithi_top200 = []
#syfpeithi_map = ep.getMapwithValues("syfpeithi.txt")
seq = ep.getSeq(inputFile)
syfpeithi_map = ep.parseSypethi(argument[1],argument[2],seq)
if(syfpeithi_map != -1):

  sorted_syfpeithi = sorted(syfpeithi_map.iteritems(),key=operator.itemgetter(1),reverse=True)

  syfpeithi_top200 = ep.getTop200(sorted_syfpeithi)
  syfpeithi_norm_map = ep.writeNormFile(rFilePath, outputFilePath ,"Syfpethi",syfpeithi_top200,syfpeithi_map)

###################################################
### Final Set with top 200
###################################################

final_set=[]
bp_key_ = []
bp_value_ = []
if(syfpeithi_map != -1):
  final_set = set(IEDB_top200).intersection(set(netmhcpan_top200).intersection(syfpeithi_top200))
else:
  final_set = set(IEDB_top200).intersection(set(netmhcpan_top200))
fwrite = open(outputFilePath+'transform_input.csv','w')
if(syfpeithi_map != -1):
  fwrite.write("Epitope,IEDB.Norm,NetMHC.Norm,Syfpethi.Norm,BP Score\n")
else:
  fwrite.write("Epitope,IEDB.Norm,NetMHC.Norm,BP Score\n")

for val in final_set:
     val = val.replace("'","")
     if(syfpeithi_map != -1):
      bp_score = (IEDB_norm_map[val]  + netmhcpan_norm_map[val] + syfpeithi_norm_map[val] ) / 3 ;
      bp_key_.append(val)
      bp_value_.append(bp_score)
      fwrite.write(val + "," + str(IEDB_norm_map[val]) + "," + str(netmhcpan_norm_map[val]) + "," + str(syfpeithi_norm_map[val]) + "," + str(bp_score) +"\n" )
     else:
      bp_score = (IEDB_norm_map[val]  + netmhcpan_norm_map[val]  ) / 2 ;
      bp_key_.append(val)
      bp_value_.append(bp_score)
      fwrite.write(val + "," + str(IEDB_norm_map[val]) + "," + str(netmhcpan_norm_map[val]) + ","  + str(bp_score) +"\n" )

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
     print "The ANN did not converge, running the ANN again for" +argument[4]+ " - " + argument[1]
     print returncode
     commandR= "Rscript "+rFilePath+"/ANN_Immunogenicity.r " + outputFilePath + " " + rFilePath
     returncode = subprocess.call(commandR,shell=True)
ann_prob = ep.getAnnProb(outputFilePath)
ep.writeCombinedProb(final_set,bp_map,ann_map,ann_prob,outputFilePath)
commandSort = "sort -nk2 "+outputFilePath+"Epitope_prob.txt > "+outputFilePath+"Epitope_FinalCombined_Prob.txt"
returncode = subprocess.call(commandSort,shell=True)

###################################################
# End of Program
###################################################
