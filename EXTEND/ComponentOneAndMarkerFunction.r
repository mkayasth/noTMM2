ComponentOneAndMarkerFunction <- function(DataRanks2,TERTrank,Datasize){

DataRanks3 = DataRanks2
DataRanks3 = 10000*DataRanks3/Datasize
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
TERTrank = match("TERT",rownames(DataRanks3))
Correl3 = as.numeric(cor.test(as.numeric(CumSum1$RawRankSum),as.numeric(DataRanks3[TERTrank,]),method='spearman')$estimate)   ### correlation of TERTrank and the Cumsum
PenalizingFactor = 1-Correl3    
DataRanks3[TERTrank,] = DataRanks3[TERTrank,]/PenalizingFactor

CumSum = data.frame(colSums (DataRanks3, na.rm = FALSE, dims = 1))  
SampleID = rownames(CumSum)

CumSum3 = data.frame(SampleID,CumSum)
row.names(CumSum3)=NULL
colnames(CumSum3)[2]='RawEXTENDScores'

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


} ### ComponentOneAndMarkerFunction ENDS HERE
