#' @title ancfreq
#' @description A function to calculate ancestral allele frequency based on two reference populations.
#' @param inputdir_dip The dir to a file containing diploid genotype likelihood data.
#' @param pop1_dip The dir to a file containing information of reference population 1 for diploid data.
#' @param pop2_dip The dir to a file containing information of reference population 2 for diploid data.
#' @param inputdir_hap The dir to a file containing haploid genotype likelihood data.
#' @param pop1_hap The dir to a file containing information of reference population 1 for haploid data.
#' @param pop2_hap The dir to a file containing information of reference population 2 for haploid data.
#' @param outputdir The dir to the output file. A postfix ('_dip','_hap','_merge') will be added.
#' @param mergetype A parameter for merging type ("intersect","union","full_dip","full_hap")
#' @return Nothing
#' @export
#'
# @examples
ancfreq<-function(outputdir,inputdir_dip="missing",pop1_dip="missing",pop2_dip="missing",
                  inputdir_hap="missing",pop1_hap="missing",pop2_hap="missing",
                  mergetype="union"){

  #Check input
  flag=0
  if(inputdir_dip!='missing'){
    flag=datacheck(inputdir_dip,character='[^0-9.[:space:]-]',field = '2- -d \" \"')
    if(pop1_dip=='missing'){
      print("Please provide input for diploid reference population 1.")
      flag=1
    }
    else if(pop2_dip=='missing'){
      print("Please provide input for diploid reference population 2.")
      flag=1
    }
    else{
    input=file(inputdir_dip,'r')
    line=unlist(strsplit(readLines(input,n=1),split = ' '))
    #print(line)
    indnum_dip=(length(line)-2)/3
    if(floor(indnum_dip)!=indnum_dip || indnum_dip < 1){
      print(paste("Unexpectedly finding",indnum_dip,"individuals in the file,
                  please double check your input format"))
      flag = 1
    }
    else{
      print(paste("Finding",indnum_dip,"individuals in the file"))
      close(input)
      outdip=paste(outputdir,'_dip',sep = '')
    }
    }
  }
  if(flag==0){
    if(inputdir_hap!='missing'){
      flag=datacheck(inputdir_hap,character='[^0-9.[:space:]-]',field = '2- -d \" \"')
      if(pop1_hap=='missing'){
        print("Please provide input for haploid reference population 1.")
        flag=1
      }
      else if(pop2_hap=='missing'){
        print("Please provide input for haploid reference population 2.")
        flag=1
      }
      else{
        input=file(inputdir_hap,'r')
        line=unlist(strsplit(readLines(input,n=1),split = ' '))
        indnum_hap=(length(line)-2)/2
        close(input)
        outhap=paste(outputdir,'_hap',sep = '')
      }
    }

    #Do ancestral allele frequency calculation
    if(inputdir_dip!="missing"&&inputdir_hap!="missing"){
      print('Generating ancestral allele frequency for diploid, haploid data and merge them.')
      typenum=-1
      #"intersect","union","full_dip","full_hap"
      if(mergetype=="intersect"){
        typenum=0
        print("Taking only the intersect of snps during merging.")
      }
      else if(mergetype=="full_dip"){
        typenum=2
        print("Keeping all snps in diploid data during merging.")
      }
      else if(mergetype=="full_hap"){
        typenum=1
        print("Keeping all snps in haploid data during merging.")
      }
      else if(mergetype=="union"){
        typenum=3
        print("Keeping the union of snps during merging.")
      }
      else{
        print(paste("Merge type",mergetype,"not found, please choose from intersect,union,full_dip and full_hap"))
      }
      if(typenum>=0){
        outmerge=paste(outputdir,'_merge',sep = '')
        #int n, int type, std::string pop1, std::string pop2, std::string input, std::string output
        ancfreq_c(n=indnum_dip,type=2,pop1=pop1_dip,pop2=pop2_dip,input=inputdir_dip,output=outdip)
        ancfreq_c(n=indnum_hap,type=1,pop1=pop1_hap,pop2=pop2_hap,input=inputdir_hap,output=outhap)
        #std::string hapfreq, std::string dipfreq, std::string outputdir, int type
        ancfreq_merge(hapfreq=outhap,dipfreq=outdip,outputdir=outmerge,type=typenum)
      }
    }
    else if(inputdir_dip!='missing'){
      print('Generating ancestral allele frequency for diploid data only.')
      ancfreq_c(n=indnum_dip,type=2,pop1=pop1_dip,pop2=pop2_dip,input=inputdir_dip,output=outdip)
    }
    else if(inputdir_hap!='missing'){
      print('Generating ancestral allele frequency for haploid data only.')
      ancfreq_c(n=indnum_hap,type=1,pop1=pop1_hap,pop2=pop2_hap,input=inputdir_hap,output=outhap)
    }
    else{
      print("The input can't be missing for both types.")
    }
  }
}

#' @title mergefreq
#' @description A function to merge existing ancestral frequency files
#' @param hap The ancestral frequency calculated for haploid data
#' @param dip The ancestral frequency calculated for diploid data
#' @param merge The output directory for merged ancestral frequency
#' @param type A parameter for merging type ("intersect","union","full_dip","full_hap"). 'union' for default.
#'
#' @return A merged ancfreq file
#' @export
#'
# @examples
mergefreq<-function(hap,dip,merge,type='union'){
  dict=c("union"=3,"full_dip"=2,"full_hap"=1,"intersect"=0)
  typenum=as.numeric(dict[type])
  #print(typenum)
  ancfreq_merge(hapfreq=hap,dipfreq=dip,outputdir=merge,type=typenum)
}




