InputData <- function(data){

if(is.matrix(data)){

Name = as.character(rownames(data))

data = data.frame(Name,data)

row.names(data)=NULL

IterativeRS(data) ## calling IterativeRS Function

}else{

print("Wrong Input: Please Enter the data in matrix form where Row Names correspond to genes and column Names correspond to Samples")

}
} ### Input function ends here
