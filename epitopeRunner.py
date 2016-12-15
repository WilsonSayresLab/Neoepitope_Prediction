from os import listdir
from os.path import isfile, join, dirname
import os
import sys
import subprocess

## Get the current directory
filepath = os.path.dirname(os.path.realpath(__file__))
hlapath = filepath+ "/hla"
allfiles = [f for f in listdir(hlapath) if isfile(join(hlapath, f))]
allfiles = set(allfiles)

file_hlas_map =dict();
for file in allfiles:
     file_hlas_map[file] = list();
     with open(hlapath+ "/"+ file) as f:
          for line in f:
               hlas = line.strip().split("\t")
               file_hlas_map[file].extend(list(set(hlas[1:len(hlas)])))

#print file_hlas_map
procs = []
for key, value in file_hlas_map.items():
     filenameArr = key.split("-")
     patientID = filenameArr[1] +"-"+ filenameArr[2];
     filename = "VarScan2-TCGA-" + patientID + "-TNBC.peptide"
     for val in value:
          hlas = val.upper().split("_")
          hla = hlas[0] + "-" + hlas[1] + "*" + hlas[2] + ":" + hlas[3]
          proc = subprocess.Popen([sys.executable, 'epitope.py', hla, '9', filename, patientID])
          procs.append(proc)
for proc in procs:
     proc.wait()
