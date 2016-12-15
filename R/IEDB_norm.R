#data= read.table("ConsensusGag (1-200) H2Db revised.csv",header=T,sep=",")
#data= read.table("Combined H2Db Glyco(1-200) weighted.csv",header=T,sep=",")
#data= read.table("SMMpmbec-Comimb(prot)-IEDB-NetMHCpan-Syfpethi(bind)1-200 weighted.csv",header=T,sep=",")
#data= read.table("E7_B8(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E7_B35(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E7_B7(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E7_A11(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E7_A24 (1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E6_B8(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E6_B7(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E6_A11(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E6_B35(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("E6_A24(1-200) combined pre-immuno.csv",header=T,sep=",")
#data= read.table("FluA-Neuramin Db combined pre-immuno hydro.csv",header=T,sep=",")
#--data= read.table("LCMV-NP(1-200) combined hydro.csv",header=T,sep=",")
args <- commandArgs(trailingOnly = TRUE)

filename = paste(args[1], "IEDB_norm.csv", sep="")
data= read.table(filename,header=T,sep=",")

# IEDB binding normalization (lower the better)
amin = min(data$IEDB.bind)
amax = max(data$IEDB.bind)
Norm.IEDB.bind = ((data$IEDB.bind-amin)/(amax-amin))

# IEDB proteasome cleavage normalization (higher the better) - using absolute values
#dmax = max(abs(data$IEDB.prot))
#dmin = min(abs(data$IEDB.prot))
#Norm.IEDB.prot = (abs(data$IEDB.prot)-dmin)/abs((dmax-dmin))

# ANN proteasome cleavage normalization (higher the better) - using absolute values
#emax = max(abs(data$ANN.prot))
#emin = min(abs(data$ANN.prot))
#Norm.ANN.prot = (abs(data$ANN.prot)-emin)/abs((emax-emin))

# Comimb proteasome cleavage normalization (higher the better)
#fmax = min(data$Comimb.prot)
#fmin = max(data$Comimb.prot)
#Norm.Comimb.prot = ((data$Comimb.prot-fmin)/(fmax-fmin))

# SMMPMBC proteasome cleavage normalization (higher the better) - using absolute values
#gmax = max(abs(data$Smmpmbc.prot))
#gmin = min(abs(data$Smmpmbc.prot))
#Norm.Smmpmbc.prot = (abs(data$Smmpmbc.prot)-gmin)/abs((gmax-gmin))

data$IEDB.Norm = Norm.IEDB.bind
#data$IEDB.prot.Norm = Norm.IEDB.prot
#data$ANN.prot.Norm = Norm.ANN.prot
#data$Comimb.prot.Norm = Norm.Comimb.prot
#data$Smmpmbc.prot.Norm = Norm.Smmpmbc.prot
filename = paste(args[1], "IEDB_normalized.txt", sep="")
write.table(data,file =filename, sep = "\t", col.names = NA, qmethod = "double")

