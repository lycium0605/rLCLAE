
#' anccall
#'
#' @param delta The cutoff for ancestral allele frequency difference
#' @param window_size The size of the sliding window
#' @param inputdir The file to the ancestral likelihood file
#' @param outputdir The file to store the ancestry call
#'
#' @return nothing
#' @export
#'
# @examples
anccall<-function(delta=0.2, window_size=50000, inputdir,outputdir){

  #Calculating snp number in the input
  system(paste("echo $(wc -l",inputdir,")"), intern = TRUE)->lines
  nline<-unlist(strsplit(lines,split = ' '))
  smax=as.numeric(nline[1])

  #Call ancestry
  anccall_c(deltaf=delta,window=window_size,SMAX=smax,anclikdir=inputdir,output=outputdir)
}

test_anccall<-function(){
  testdir="/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/anclik_R"
  genolik2="/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test.call_R"
  anccall(inputdir=testdir,outputdir=genolik2)
}
