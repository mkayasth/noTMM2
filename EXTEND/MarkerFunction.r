MarkerFunction <- function(DataRanks2,Datasize){

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
colnames(CumSumAll)[2:3]=c('RawEXTENDScores','NormEXTENDScores')


write.table(CumSumAll,'TelomeraseScores.txt',sep='\t',quote=FALSE,row.names=FALSE)
print("Finished Computation of Telomerase Activity")

}### MarkerFunction ENDS HERE

