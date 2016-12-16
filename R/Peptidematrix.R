args <- commandArgs(trailingOnly = TRUE)
filename = paste(args[1], "transform_input.csv", sep="")
#data= read.table("H2Db CTL-test.csv",header=T,sep=",")
#data= read.table("H2Kb eluted-test.csv",header=T,sep=",")
#data= read.table("human-mouse self CTL.csv",header=T,sep=",")
#data = read.table("All pathogens H2Db CTL.csv",header=T,sep=",")
#data= read.table("htlvTAX.pep.csv",header=T,sep=",")
#data = read.table("H2Db full.csv",header=T,sep=",")
#data = read.table("H2Db training.csv",header=T,sep=",")
#ata = read.table("H2Db testing.csv",header=T,sep=",")
#data = read.table("Glyco-combined pre-immuno.csv",header=T,sep=",")
#data = read.table("A2 full.csv",header=T,sep=",")
#data = read.table("A2 training.csv",header=T,sep=",")
#data = read.table("A2 testing.csv",header=T,sep=",")
#data = read.table("H2Db CTL pos 9mers.csv",header=T,sep=",")
#data = read.table("ConsensusGag (1-200) H2Db combined weighted.csv",header=T,sep=",")
#data = read.table("LCMV-NP (1-200) H2Db combined weighted-1.csv",header=T,sep=",")
#data = read.table("T antigen (1-200) combined hydro.csv",header=T,sep=",")
#data = read.table("FLuA-NP_H2Db 9mers combined.csv",header=T,sep=",")
#data = read.table("97CN54Gag (1-200) H2Db combined weighted.csv",header=T,sep=",")
#data = read.table("ZM96Gag (1-200) combined weighted.csv",header=T,sep=",")

data = read.table(filename,header=T,sep=",")

filename = paste(args[2], "HydrophobicityTable.csv", sep="")

HydrophoTable = read.table(filename,header=T,sep=",")
#HydrophoTable = read.table("Weighted Hydro table.csv.csv",header=T,sep=",")

Epitope=data$Epitope



a=vector()

n=nrow(data)

for (i in 1:n){
  
  chars_positives=strsplit(gsub("([[:alnum:]]{1})", "\\1 ",Epitope[i]), " ")[[1]]
  
  a=c(a,length(chars_positives))
  
}

aa=sort(unique(a))
matrixx=list()

for (k in 1:length(aa)){
  
  index_peptides=which(aa[k]==a)
  peptides=Epitope[index_peptides]
  
  cc=vector()
  
  for (i in 1:length(peptides)){
    
    chars=strsplit(gsub("([[:alnum:]]{1})", "\\1 ", peptides[i]), " ")[[1]]
    
    for (j in 1:length(chars)){
      
      index_amino=which(chars[j]==HydrophoTable$Amino.acid)
      chars[j]=HydrophoTable$Hydrophobicity[index_amino]
      
    }
    
    cc=c(cc,chars)
    cc=as.numeric(cc)
    
  }
  
  matrix=matrix(cc,nrow=aa[k],ncol=length(peptides))
  matrix=t(matrix)
  
  matrixx[[k]]=data.frame(Peptide=data$Epitope[index_peptides],matrix)
  #matrixx[[k]]=data.frame(Name=data$Name[index_peptides],Peptide=data$Epitope[index_peptides],matrix)
  
}
filename = paste(args[1], "Ann_input.txt", sep="")
write.table(matrixx[[1]],file =filename, sep = "\t", col.names = NA, qmethod = "double")
