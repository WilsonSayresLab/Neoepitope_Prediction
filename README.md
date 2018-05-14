# Neoepitope Prediction

Here we will lay out our pipeline for inferring neoepitopes in cancer. 
Steps
1) Create a direcdory “mkdir NeoepitopePrediction”
2) Download the IEDB tool from http://tools.iedb.org/mhci/download/
3) untar the folder "tar -xvf IEDB_MHC_I-2.15.1.tar.gz”
4) configure the tool cd mhc_i;./configure
cd ..
5) clone the github source code - git clone https://github.com/WilsonSayresLab/Neoepitope_Prediction.git
6) Add all the peptide files to the peptide folder in the downloaded git folder. 
7) Add about 5-7 hla to the hla folder 
8) Run the epitopeRunner.py by - python epitopeRunner.py. 
