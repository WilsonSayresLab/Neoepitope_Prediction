from os import listdir
from os.path import isfile, join, dirname
import os
import sys
import subprocess
from datetime import datetime
from collections import defaultdict

## Get the current directory
filepath = os.path.dirname(os.path.realpath(__file__)) + "/"
hlapath = filepath+ "hla"
allfiles = [f for f in listdir(hlapath) if isfile(join(hlapath, f)) and f.endswith(".txt")]
allfiles = set(allfiles)


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

procs = []
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     filename = "VarScan2-TCGA-" + patientID + "-TNBC.peptide"
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]
          proc = subprocess.Popen([sys.executable, 'epitope.py', hla, '9', filename, outputFolder + patientID])
          procs.append(proc)
for proc in procs:
     proc.wait()

#############################################
# Create Sequence transcript map
#############################################
seq_transcript_map = defaultdict(list)
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     filename = "VarScan2-TCGA-" + patientID + "-TNBC.peptide"
     transcript = ""
     counter = 0
     with open(filepath + "peptides/" + filename ) as f:
          for line in f:
               seq=""
               if ">" in line:
                    transcript = line.strip().replace(">","")
                    counter += 1
                    print transcript
               else:
                    seq = line.strip()
                    counter += 1
               if counter % 2 ==0:
                    seq_transcript_map[seq].append(transcript)

#############################################
# Create Sequence transcript map
#############################################
epitope_score_map = dict()
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     fwrite = open(outputFilePath+ patientID + "/"+patientID+ '.csv','w')
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]
          hla = hla.replace(":","-")
          filename = "Epitope_FinalCombined_Prob.txt"
          key_ =[]
          value_ =[]
          epitope_file_path  = outputFilePath + patientID + "/"+  hla + "/" + filename
          with open(epitope_file_path ) as f:
               for line in f:
                    vals = line.split("\t")
                    key_ = vals[0]
                    value_ = vals[1].strip()
                    for seq_key in seq_transcript_map.keys():
                         if vals[0] in seq_key:
                              transcript = seq_transcript_map.get(seq_key)
                              for t in transcript:
                                   fwrite.write( patientID + "\t" + t + "\t" +  key_ + "\t" + hla + "\t"  + "\t" + value_ + "\n")
     fwrite.close()
#############################################
# End of Program
#############################################     