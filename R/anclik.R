#void ancfreq_c(int sum, int type, int testid,
#std::string genolik, std::string ancfreq, std::string output)

#' anclik
#'
#' @param genodir genotype likelihood file
#' @param ancfreqdir ancestral allele frequency file
#' @param outputdir ancestral state estimation
#' @param type 'hap' or 'dip'
#' @param test can be 'all' or a numeric vector containing the numeric index (e.g. c(1,2,3) will represent the first, second and third individuals in the vcf file) of individuals to be estimated. A suffix "_n" will be added to each output file.
#'
#' @return nothing
#' @export
#'
# @examples
anclik<-function(genodir,ancfreqdir,outputdir,type='dip',test='all'){
  flag=datacheck(ancfreqdir,'[^0-9.[:space:]]')
  if(!file.exists(genodir)){
    stop("The genotype likelihood file you provided does not exist. Please check.")
    flag=1
  }
  if(!file.exists(ancfreqdir)){
    stop("The ancestral frequency file you provided does not exist. Please check.")
    flag=1
  }
  input=file(genodir,'r')
  line=unlist(strsplit(readLines(input,n=1),split = ' '))
  on.exit(close(input))
  if(type=='dip'){
    typenum=2
    indnum=(length(line)-2)/3
    if(floor(indnum)!=indnum || indnum < 1){
      stop(paste("Unexpectedly finding",indnum,"individuals in the file,
                  please double check your input format"))
      flag = 1
    }
  }
  else if(type=='hap'){
    typenum=1
    indnum=(length(line)-2)/2
    if(floor(indnum)!=indnum || indnum < 1){
      stop(paste("Unexpectedly finding",indnum,"individuals in the file,
                  please double check your input format"))
      flag = 1
    }
  }
  else{
    stop("Wrong data type, should be either hap or dip.")
    flag = 1
  }

  if(flag==0){
    if(length(test==1)&&test=='all'){
      ind_seq<-1:indnum
      for(i in 1:indnum){
        outputdir_=paste(outputdir,i,sep = '_')
        anclik_c(sum=indnum,type=typenum,testid=i,
                 genolik=genodir,ancfreq=ancfreqdir,output=outputdir_)
      }
    }
    else{
      ind_seq<-test
      for(i in test){
        i_=as.numeric(i)
        if(i_<=indnum){
          outputdir_=paste(outputdir,i_,sep = '_')
          anclik_c(sum=indnum,type=typenum,testid=i_,
                   genolik=genodir,ancfreq=ancfreqdir,output=outputdir_)
        }
      }
    }

  }
return(ind_seq)
}

