
#' @title  preprocess
#' @description  A function to preprocess the vcf file
#'
#' @param inputdir dir to orginial vcf file
#' @param outputdir dir to store cleaned_up vcf file
#'
# @return
#' @export
#'

preprocess<-function(inputdir,outputdir){
  clean="sed 's/|/\\//g' | sed '/^#/d' | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\\./999/g'"
  command_=paste(sep='','cat ',inputdir, ' | ',clean,'>',outputdir)
  system(command=command_,intern = FALSE)
  inputcheck(outputdir)
}

#' @title inputcheck
#' @description A function to check the format of the input cleaned-up vcf file
#' @param inputdir dir to clean vcf
#'
#' @return A string indicating the result of the input check
#' @export
#'

inputcheck<-function(inputdir){
  #Check characters
  com = '| cut -f 2- | grep [^0-9,[:space:]/] > /dev/null&& echo "Unexpected character, please double check your input" || echo "Pass character test"'
  command_ = paste(sep=' ', 'cat',inputdir,com)
  system(command = command_,intern = TRUE) -> charcheck
  #return(charcheck)
  #message(charcheck)
  if(charcheck=="Pass character test"){
    message(charcheck)
    #Check individual number
    input=file(inputdir,'r')
    line=unlist(strsplit(readLines(input,n=1),split = '\t'))
    indnum=length(line)-2
    close(input)
    if(indnum>0){
    message(paste(sep=' ',"Finding",indnum,"individuals from this vcf file."))
    }else{
      stop(paste(sep=' ',"Finding",indnum,"individuals from this vcf file."))
    }
    #Check chromosomes
    x<-get_chrlist(inputdir)
    if(length(x)<=0){
      stop("Not finding chromosomes in the processed input")
    }else{
      chr<-x[seq(2,length(x),2)]
      rep<-as.numeric(x[seq(1,length(x),2)])
      message(paste(sep = ' ',"Finding",length(chr),"unique chromosomes in this vcf file, Including:"))
      for(i in 1:length(chr)){
        message(paste(rep[i],"snps in",chr[i]))
    }
  }
  }else{
    stop(charcheck)
  }
  }
