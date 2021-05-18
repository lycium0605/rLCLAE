#' ancfreq
#'
#' @param inputdir_dip The dir to a file containing diploid genotype likelihood data.
#' @param pop1_dip The dir to a file containing information of reference population 1 for diploid data.
#' @param pop2_dip The dir to a file containing information of reference population 2 for diploid data.
#' @param inputdir_hap The dir to a file containing haploid genotype likelihood data.
#' @param pop1_hap The dir to a file containing information of reference population 1 for haploid data.
#' @param pop2_hap The dir to a file containing information of reference population 2 for haploid data.
#' @param outputdir The dir to the output file. A postfix ('_dip','_hap','_merge') will be added.
#'
# @return Nothing
# @export
#'
# @examples
ancfreq<-function(inputdir_dip="missing",pop1_dip="missing",pop2_dip="missing",
                  inputdir_hap="missing",pop1_hap="missing",pop2_hap="missing",
                  outputdir){

  #Check input
  if(inputdir_dip!='missing'){
    if(pop1_dip=='missing'){
      print("Please provide input for diploid reference population 1.")
      #quit()
    }
    if(pop2_dip=='missing'){
      print("Please provide input for diploid reference population 2.")
      #quit()
    }
    input=file(inputdir_dip,'r')
    line=unlist(strsplit(readLines(input,n=1),split = ' '))
    #print(line)
    indnum_dip=(length(line)-2)/3
    print(paste("Finding",indnum_dip,"individuals in the file"))
    close(input)
    outdip=paste(outputdir,'_dip',sep = '')
  }
  if(inputdir_hap!='missing'){
    if(pop1_hap=='missing'){
      print("Please provide input for haploid reference population 1.")
      #quit()
    }
    if(pop2_hap=='missing'){
      print("Please provide input for haploid reference population 2.")
      #quit()
    }
    input=file(inputdir_hap,'r')
    line=unlist(strsplit(readLines(input,n=1),split = ' '))
    indnum_hap=(length(line)-2)/2
    close(input)
    outhap=paste(outputdir,'_hap',sep = '')
  }

  #Do ancestral allele frequency calculation
  if(inputdir_dip!="missing"&&inputdir_hap!="missing"){
    print('Generating ancestral allele frequency for diploid, haploid data and merge them.')
    outmerge=paste(outputdir,'_merge')
    #int n, int type, std::string pop1, std::string pop2, std::string input, std::string output
    #ancfreq_c(n=indnum_dip,type=2,pop1=pop1_dip,pop2=pop2_dip,input=inputdir_dip,output=outdip)
    #ancfreq_c(n=indnum_hap,type=1,pop1=pop1_hap,pop2=pop2_hap,input=inputdir_hap,output=outhap)
  }
  else if(inputdir_dip!='missing'){
    print('Generating ancestral allele frequency for diploid data only.')
    ancfreq_c(n=indnum_dip,type=2,pop1=pop1_dip,pop2=pop2_dip,input=inputdir_dip,output=outdip)
  }
  else if(inputdir_hap!='missing'){
    print('Generating ancestral allele frequency for haploid data only.')
    #ancfreq_c(n=indnum_hap,type=1,pop1=pop1_hap,pop2=pop2_hap,input=inputdir_hap,output=outhap)
  }
  else{
    print("The input can't be missing for both types.")
    #quit()
  }
}

test_ancfreq<-function(){
  data="/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/"
  input=paste(data,"test.genolik_R",sep = '')
  output=paste(data,"ancfreq",sep = '')
  ref1=paste(data,"fullref_anubis.h",sep = '')
  ref2=paste(data,"fullref_yellow.h",sep = '')
  ancfreq(inputdir_dip = input,pop1_dip = ref2,pop2_dip = ref1,outputdir = output)
  }

