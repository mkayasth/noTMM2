ComponentAndMarkerFunction <- function(DataRanks2,TERTrank,TERCrank,Datasize){


TM_Facts = DataRanks2[c(TERTrank,TERCrank),]
TM_Facts = t(TM_Facts)


TM_FactsMaxAll = data.frame(rownames(TM_Facts)[1],max(TM_Facts[1,1],TM_Facts[1,2]))
colnames(TM_FactsMaxAll)=c('SampleID','Max')
for(i in 2:nrow(TM_Facts)){
TM_FactsMax = data.frame(rownames(TM_Facts)[i],max(TM_Facts[i,1],TM_Facts[i,2]))
colnames(TM_FactsMax)=c('SampleID','Max')
TM_FactsMaxAll = rbind(TM_FactsMaxAll,TM_FactsMax)
}

DataRanks3 = DataRanks2[-c(TERCrank,TERTrank),]
DataRanks3 = t(DataRanks3)
DataRanks3 = as.data.frame(DataRanks3)
SampleID = rownames(DataRanks3)
DataRanks3 = data.frame(SampleID,DataRanks3)
row.names(DataRanks3) = NULL
DataRanks3 = merge(TM_FactsMaxAll,DataRanks3,by='SampleID')

row.names(DataRanks3) = DataRanks3$SampleID
DataRanks3 = DataRanks3[,-1]
DataRanks3 = t(DataRanks3)

DataRanks3 = 10000*DataRanks3/nrow(DataRanks3)
CumSum = data.frame(colSums(DataRanks3, na.rm = FALSE, dims = 1))    ### Cumsum of each sample calculated based on the ranks
colnames(CumSum)[1]='RawRankSum'
SampleID = rownames(CumSum)
CumSum1 = data.frame(SampleID,CumSum)
row.names(CumSum1)=NULL

DataRanks4 = DataRanks3
DataRanks4 = DataRanks4/nrow(DataRanks4)

CumSum = data.frame(colSums(DataRanks4, na.rm = FALSE, dims = 1))    ### Cumsum of each sample calculated based on the ranks
colnames(CumSum)[1]='NormRankSum'
SampleID = rownames(CumSum)
CumSum2 = data.frame(SampleID,CumSum)
row.names(CumSum2)=NULL
CumSum2$NormRankSum = (CumSum2$NormRankSum-min(CumSum2$NormRankSum))/(max(CumSum2$NormRankSum)-min(CumSum2$NormRankSum))
CumSumAll = merge(CumSum1,CumSum2,by='SampleID')



Maxrank = match("Max",rownames(DataRanks3))
Correl3 = as.numeric(cor.test(as.numeric(CumSum1$RawRankSum),as.numeric(DataRanks3[Maxrank,]),method='spearman')$estimate)   ### correlation of Maxrank and the Cumsum
#print(Correl3)
PenalizingFactor = 1-Correl3    ### Maxrank PenalizingFactor
DataRanks3[Maxrank,] = DataRanks3[Maxrank,]/PenalizingFactor

CumSum = data.frame(colSums (DataRanks3, na.rm = FALSE, dims = 1))  
SampleID = rownames(CumSum)

CumSum3 = data.frame(SampleID,CumSum)
row.names(CumSum3)=NULL
colnames(CumSum3)[2]='RawEXTENDScores'
#CumSumAll = merge(CumSumAll,CumSum3,by='SampleID')

DataRanks5 = DataRanks3
DataRanks5 = DataRanks5/nrow(DataRanks5)

CumSum = data.frame(colSums(DataRanks5, na.rm = FALSE, dims = 1))    ### Cumsum of each sample calculated based on the ranks
colnames(CumSum)[1]='NormEXTENDScores'
SampleID = rownames(CumSum)
CumSum4 = data.frame(SampleID,CumSum)
row.names(CumSum4)=NULL
CumSum4$NormEXTENDScores = (CumSum4$NormEXTENDScores-min(CumSum4$NormEXTENDScores))/(max(CumSum4$NormEXTENDScores)-min(CumSum4$NormEXTENDScores))
CumSumFinal = merge(CumSum3,CumSum4,by='SampleID')


write.table(CumSumFinal,'TelomeraseScores.txt',sep='\t',quote=FALSE,row.names=FALSE)
print("Finished Computation of Telomerase Activity")


}### ComponentAndMarkerFunction ends here
