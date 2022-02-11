
#' @title anctract
#' @description Generating tract from call
#' @param inputdir input file directory
#' @param outputdir output file directory
#' @param chrlength length of the chrmosome
#' @param exclude length to be excluded from both end
#'
# @return
#' @export
#'
# @examples
anctract<-function(inputdir,outputdir,chrlength,exclude=0){
  anctract_c(chrlen = chrlength,
             excludelen = exclude,
             input = inputdir,
             output = outputdir)
}