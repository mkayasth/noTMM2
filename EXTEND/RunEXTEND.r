#' EXTEND(EXpression based Telomerase ENzyme Detection)
#'
#' This function takes the Gene expression data (in matrix format) as an input, where the Row Names correspond to Gene Names (as Gene Symbols) and the column Names correspond to sample Names.
#'
#' RunEXTEND function processes the input matrix for calculation of Telomerase activity scores via internal modules.
#' @param infiles: Input = Gene expression matrix.Use RunEXTEND(data) command to execute the EXTEND function.
#' @return Output would be generated in the form of a text file containing raw and normalized scores for Telomerase activity prediction.
#' @export
RunEXTEND = function(data){

InputData(data)

}### EXTEND ENDS HERE
