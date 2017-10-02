rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

filename = paste(args[2], "09.14 H2Db training p5N hydro mtrx.csv", sep="")
# Training and testing data sets
#data_training = read.table("H2Db training mtrx2.csv",header=T,sep=",")
#data_testing = read.table("H2Db testing mtrx2.csv",header=T,sep=",")
#data_training = read.table("H2Db training-prediction-self mtrx.csv",header=T,sep=",")
#data_training = read.table("H2Db only NPpred-training mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db training p5N hydro mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db training p5N hydro mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db training p5N hydro mtrx2.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db training p5N hydro mtrx2.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db 2.3 training-DFCI(pos.incl) mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N training-DFCI mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N training-DFCI nop5N-NB mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N 2.3 training-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full 2.3 training-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N training-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N-DFCI NB-nop5N pre-valid-train mtrx.csv",header=T,sep=",")
##data_training = read.table("09.14 H2Db full training-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
#data_training = read.table("data_training.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full training-DFCI NB-nop5N no Glyco mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db p5N training-DFCI NB-nop5N no FluNP mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full training-DFCI-full NB no glyco mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full training-DFCI NB-nop5-N-others no glyco mtrx.csv",header=T,sep=",")
#data_training = read.table("data_training_noglyco.csv",header=T,sep=",")
#data_training = read.table("data_training_H2Db p5N(almost) no NB no glyco.csv",header=T,sep=",")
#data_training = read.table("data_training_H2Db p5N(almost) no NB no fluNP.csv",header=T,sep=",")                                                     
#data_training = read.table("data_training_H2Db full training-no NB no fluNP.csv",header=T,sep=",")
#data_training = read.table("H2Db p5N no NB no glyco mtrx.csv",header=T,sep=",")
#data_training = read.table("03.10 H2Db p5N no NB no fluNP mtrx.csv",header=T,sep=",")
#data_training = read.table("data_training.csv",header=T,sep=",")
#data_training = read.table("datatrainingtantigen.csv",header=T,sep=",")
#data_training = read.table("data_train.csv",header=T,sep=",")
data_training = read.table(filename,header=T,sep=",")
#data_training = read.table("03.12 H2Db p5N training-p5N NB mtrx.csv",header=T,sep=",")
#data_training = read.table("03.12 H2Db p5N training-p5N NB no glyco mtrx.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full training-DFCI NB-nop5-N-others no Flu-NP mtrx.csv",header=T,sep=",")
#data_training = read.table("data_training_glyco.csv",header=T,sep=",")
#data_training = read.table("09.14 H2Db full training-DFCI NB-nop5N no FluNP mtrx.csv",header=T,sep=",")
#2data_training = read.table("09.14 H2Db training full hydro mtrx.csv",header=T,sep=",")
#3data_training = read.table("H2Db full matrix.csv",header=T,sep=",")   .
#data_testing = read.table("H2Db only NPpred-testing mtrx.csv",header=T,sep=",")
#data_testing = read.table("H2Db full-NPpred-testing mtrx.csv",header=T,sep=",")
#data_testing = read.table("H2Db testing-prediction-self mtrx.csv",header=T,sep=",")
#data_testing = read.table("Glyco hydro.csv",header=T,sep=",")                                   
##data_testing = read.table("data_testing_glyco.csv",header=T,sep=",")
#data_testing = read.table("Tantigen hydro.csv",header=T,sep=",")
#data_testing = read.table("NP hydro.csv",header=T,sep=",")
#data_testing = read.table("Flu_NP hydro.csv",header=T,sep=",")
#data_testing = read.table("ConsensusGag hydro.csv",header=T,sep=",")
#data_testing = read.table("ConsensusGag hydro2.csv",header=T,sep=",")
####data_testing = read.table("FluA-Neuramin Db pre-immuno hydro mtrx.csv",header=T,sep=",")
#data_testing = read.table("09.14 H2Db 1.3 testing-DFCI(pos.incl) mtrx.csv",header=T,sep=",")
#data_testing = read.table("09.14 H2Db p5N 1.3 testing-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
##data_testing = read.table("09.14 H2Db full 1.3 testing-DFCI NB-nop5N mtrx.csv",header=T,sep=",")
#data_testing = read.table("97CN54Gag hydro.csv",header=T,sep=",")
#data_testing = read.table("ZM96Gag hydro.csv",header=T,sep=",")
#data_testing = read.table("ConsensusGag (1-200) H2Db pre-immuno Apr01.csv",header=T,sep=",")
#data_testing = read.table("09.14 H2Db testing p5N hydro mtrx.csv",header=T,sep=",")
#data_testing = read.table("ConsensusGag hydro.csv",header=T,sep=",")
#data_testing = read.table("Glyco_weighted.csv",header=T,sep=",")
#data_training = read.table("A2 training matrix.csv",header=T,sep=",")
#data_testing = read.table("A2 testing matrix.csv",header=T,sep=",")
#data_training = read.table("H2Db training weighted.csv",header=T,sep=",")
#data_testing = read.table("H2Db testing weighted.csv",header=T,sep=",")
#data_training = read.table("H2Db full matrix2.csv",header=T,sep=",")
#data_training = read.table("H2Db training mtrx 23578.csv",header=T,sep=",")
#data_testing = read.table("H2Db testing mtrx 23578.csv",header=T,sep=",")
filename = paste(args[1], "ANN_input_final.", sep="")
filename = paste(filename,args[3],sep="")
filename = paste(filename,".txt",sep="")
data_testing = read.table(filename,header=T,sep="\t")    
                                                
library(neuralnet)
#library(ROCR)
#library(pROC)
#library(ggplot)

#install.packages("neuralnet")

#vec1=vector()

#install.packages("ROCR")
#install.packages("pROC")

vec1=list()
for (j in 1:60) {            # 60 realizations of the predictions 

nnet_train = data_training
nnet_train = cbind(nnet_train, data_training$Immunogenecity == 'Positive')
nnet_train = cbind(nnet_train, data_training$Immunogenecity == 'Negative')
names(nnet_train)[11] <- 'Positive'
names(nnet_train)[12] <- 'Negative'
#nn <- neuralnet(Positive+Negative ~ X2 + X3 + X5 + X7 + X8, data=nnet_train, hidden=c(6),linear.output=FALSE,rep=1)
#nn <- neuralnet(Positive+vecNegative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data=nnet_train, hidden=c(9), linear.output=FALSE, rep=1,startweights=w)
#nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data=nnet_train, hidden=c(9), linear.output=FALSE, rep=1,startweights=w5)
if (args[3] == "21") {
	nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 ,data=nnet_train,hidden=c(3),linear.output=FALSE,rep=1)        # 3 units in the hidden layer	
}
if (args[3] == "19") {
	nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 ,data=nnet_train,hidden=c(3),linear.output=FALSE,rep=1)        # 3 units in the hidden layer	
}
if (args[3] == "17") {
	nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 ,data=nnet_train,hidden=c(3),linear.output=FALSE,rep=1)        # 3 units in the hidden layer	
}
if (args[3] == "15") {
	nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 ,data=nnet_train,hidden=c(3),linear.output=FALSE,rep=1)        # 3 units in the hidden layer	
}
#nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data=nnet_train, hidden=c(4), linear.output=FALSE,rep=1,startweights=w32)
#nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data=nnet_train, hidden=c(9), linear.output=FALSE, rep=1,startweights=w52)
#nn <- neuralnet(Positive+Negative ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data=nnet_train, hidden=c(5),linear.output=FALSE,learningrate=0.001,threshold=0.01,algorithm="backprop")
#nn$data
#nn$covariate
#nn$response
#nn$net.result
#nn$result.matrix
#out <- cbind(nn$covariate, nn$net.result[[1]])

#plot(nn)

#mypredict <- compute(nn, data_testing[-10])$net.result
#mypredict1 <- compute(nn, data_training[-10])$net.result
#mypredict2 <- compute(nn, data_testing[-10])$net.result
mypredict2 <- compute(nn, data_testing)$net.result
#out <- compute(nn, data_testing[-10])
#out$net.result

# Put multiple binary output to categorical output
maxidx <- function(arr) {
return(which(arr == max(arr)))
}

#idx1 <- apply(mypredict1, c(1), maxidx)
idx2 <- apply(mypredict2, c(1), maxidx)

#predictions1 <- c('Positive', 'Negative')[idx1]
predictions2 <- c('Positive', 'Negative')[idx2]

vec1[[j]]=mypredict2[,1]

}

vec2=vector()

for (i in 1:length(vec1[[1]])){

# average of output predictions (probabilities) from the different realizations
avg=(vec1[[1]][i]+vec1[[2]][i]+vec1[[3]][i]+vec1[[4]][i]+vec1[[5]][i]+vec1[[6]][i]+vec1[[7]][i]+vec1[[8]][i]+vec1[[9]][i]+vec1[[10]][i]+vec1[[11]][i]+vec1[[12]][i]+vec1[[13]][i]+vec1[[14]][i]+vec1[[15]][i]+vec1[[16]][i]+vec1[[17]][i]+vec1[[18]][i]+vec1[[19]][i]+vec1[[20]][i]+vec1[[21]][i]+vec1[[22]][i]+vec1[[23]][i]+vec1[[24]][i]+vec1[[25]][i]+vec1[[26]][i]+vec1[[27]][i]+vec1[[28]][i]+vec1[[29]][i]+vec1[[30]][i]+vec1[[31]][i]+vec1[[32]][i]+vec1[[33]][i]+vec1[[34]][i]+vec1[[35]][i]+vec1[[36]][i]+vec1[[37]][i]+vec1[[38]][i]+vec1[[39]][i]+vec1[[40]][i]+vec1[[41]][i]+vec1[[42]][i]+vec1[[43]][i]+vec1[[44]][i]+vec1[[45]][i]+vec1[[46]][i]+vec1[[47]][i]+vec1[[48]][i]+vec1[[49]][i]+vec1[[50]][i]+vec1[[51]][i]+vec1[[52]][i]+vec1[[53]][i]+vec1[[54]][i]+vec1[[55]][i]+vec1[[56]][i]+vec1[[57]][i]+vec1[[58]][i]+vec1[[59]][i]+vec1[[60]][i])/60
vec2=c(vec2,avg)

}

#table(predictions1, data_training$Immunogenecity)
#table(predictions2, data_testing$Immunogenecity)

#sort(mypredict2[,1],decreasing=TRUE)           

filename = paste(args[1], "Ann_output.", sep="")
filename = paste(filename,args[3],sep="")
filename = paste(filename,".txt",sep="")
write.table(vec2,file = filename, sep = "\t", col.names = NA, qmethod = "double") 
