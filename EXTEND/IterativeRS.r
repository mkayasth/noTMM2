IterativeRS <- function(data){


Datasize = nrow(data)

genes = c('C18orf54','CCT6P1','CEP72','POLE2','MCM4','HELLS','TERC','CCT6AP1','C21orf81','LOC81691','AC004381.6','ANKRD20A11P','EHMT2','LIN9','PAXIP1','TERT')

idx=match(genes,data[,1])

DataRanks = data.frame(genes,rank(data[,2])[idx])

colnames(DataRanks)=c('Genes',colnames(data)[2])

for(j in 3:ncol(data)){

sub = data.frame(rank(data[,j])[idx])   

colnames(sub)=colnames(data)[j]

DataRanks = cbind(DataRanks,sub)

}

row.names(DataRanks)=DataRanks$Genes

DataRanks = DataRanks[,-1]

DataRanks2 = DataRanks[complete.cases(DataRanks),]

TERTrank = match("TERT",rownames(DataRanks2))

TERCrank = match("TERC",rownames(DataRanks2))

if(!is.na(TERCrank) & !is.na(TERTrank)){

ComponentAndMarkerFunction(DataRanks2,TERTrank,TERCrank,Datasize)  ## calling ComponentCalculation Function1

}else if(is.na(TERCrank) & !is.na(TERTrank)){

ComponentOneAndMarkerFunction(DataRanks2,TERTrank,Datasize)  ## calling ComponentCalculation Function2

}else if(is.na(TERTrank) & !is.na(TERCrank)){

ComponentTwoAndMarkerFunction(DataRanks2,TERCrank,Datasize)  ## calling ComponentCalculation Function3

}else if(is.na(TERTrank) & is.na(TERCrank)){

MarkerFunction(DataRanks2,Datasize)  ## calling ComponentCalculation Function4

}

}#### IRS function ends here
