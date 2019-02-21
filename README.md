# Neoepitope Prediction


Here we will lay out our pipeline for inferring neoepitopes in cancer.

_**Steps**_

- Create a direcdory “mkdir NeoepitopePrediction”
- Download the IEDB tool from http://tools.iedb.org/mhci/download/
- untar the folder 
		`tar -xvf IEDB_MHC_I-2.15.1.tar.gz`
- Configure the tool 
	`cd mhc_i; ./configure`
- Change directory to exit IEDB folder `cd ..`
- Clone the github source code - `git clone https://github.com/WilsonSayresLab/Neoepitope_Prediction.git`
- Add all the peptide files to the respective peptide folder 15mers, 17mers, 19mers and 21mers in the downloaded git folder. 
- Add hlas to the hla folder 
- Run the epitopeHunter.py by - `python epitopeHunter.py`.
