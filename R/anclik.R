#void ancfreq_c(int sum, int type, int testid,
#std::string genolik, std::string ancfreq, std::string output)

#' anclik
#'
#' @param genodir genotype likelihood file
#' @param ancfreqdir ancestral allele frequency file
#' @param outputdir ancestral state estimation
#' @param type 'hap' or 'dip'
#' @param test the individual number to be estimated
#'
#' @return nothing
#' @export
#'
# @examples
anclik<-function(genodir,ancfreqdir,outputdir,type='dip',test){
  flag=0
  if(!file.exists(genodir)){
    print("The genotype likelihood file you provided does not exist. Please check.")
    flag=1
  }
  if(!file.exists(ancfreqdir)){
    print("The ancestral frequency file you provided does not exist. Please check.")
    flag=1
  }
  input=file(genodir,'r')
  line=unlist(strsplit(readLines(input,n=1),split = ' '))
  if(type=='dip'){
    typenum=2
    indnum=(length(line)-2)/3
    if(floor(indnum)!=indnum || indnum < 1){
      print(paste("Unexpectedly finding",indnum,"individuals in the file,
                  please double check your input format"))
      flag = 1
    }
  }
  else if(type=='hap'){
    typenum=1
    indnum=(length(line)-2)/2
    if(floor(indnum)!=indnum || indnum < 1){
      print(paste("Unexpectedly finding",indnum,"individuals in the file,
                  please double check your input format"))
      flag = 1
    }
  }
  else{
    print("Wrong data type, should be either hap or dip.")
    flag = 1
  }

  if(flag==0){
    anclik_c(sum=indnum,type=typenum,testid=test,
             genolik=genodir,ancfreq=ancfreqdir,output=outputdir)
  }

}

#' test_anclik
#'
#' @return nothing
#' @export
#'
# @examples
test_anclik<-function(){
  data="/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/"
  input1=paste(data,"test.genolik_R",sep = '')
  input2=paste(data,"ancfreq_dip",sep = '')
  output_=paste(data,"anclik_R",sep = '')
  anclik(genodir = input1,ancfreqdir = input2,outputdir =output_,
         type = 'dip',test = 3 )
}
