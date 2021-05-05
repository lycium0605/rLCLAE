#Step 0. filter your vcf with the linux command that looks like
#cat recode.vcf | grep ^20 | cut -f 1,2,10- | sed -e 's/:/ /g' -e 's/\./999/g' | sed -e 's/\// /g' > test_clean.vcf
#filtbaboon1b.c - Takes filtered vcf files (e.g., in1)
#and converts them into a list of genotype likelihood (PL) files.

#Inputfile Outputfile individual_num

#Output format: Each line is a SNP, consisting of the base position
#followed by the 3 genotype likelihoods (or -1. for missing data) for each
#of n different individuals.  n=47 for the example data.


#inputdir='/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test_clean.vcf'
#outputdir='/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test.genolik_R'
#options(scipen = 200)

#Roxygen commands
#' Calculating the relative genotype likelihood from a vcf file.
#'
#' @param inputdir The input directory of the cleaned vcf file.
#' @param outputdir The output directory of the genolik file.
#'
#' @return A file in which each line represents a snp. The line looks like chr, pos, -1. -1. -1. (missing data) or gl1, gl2, gl3 transformed from the phred score.


genolik<-function(inputdir,outputdir){
  #library(Rcpp)
  #Rcpp::sourceCpp('/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rcode/Aut/genolik.cpp')
  #Rcpp::sourceCpp('./R/genolik.cpp')
  input=file(inputdir,'r')
  output=file(outputdir,'w')
  line=unlist(strsplit(readLines(input,n=1),split = '\t'))
  #print(line)
  print(paste('Individuals number:',length(line)-2))
  while(length(line)!=0){
    sname=line[1]
    pos=line[2]
    newline<-c(sname,pos)
    cat(newline,file = output,sep = ' ')
    cat(' ',file = output)
    for (ind in 3:length(line)){
      indinfo=line[ind]
      indinfo<-unlist(strsplit(indinfo,split = ' '))
      if(indinfo[1]=='999/999'){
        nl1='-1.'
        nl2='-1.'
        nl3='-1.'
      }
      else{
      gl<-indinfo[length(indinfo)]
      gl<-unlist(strsplit(gl,split = ',')) # chr of 3, gl1, gl2,gl3
      gl<-as.numeric(gl)
      # Error report
      if(length(gl)!=3){
        print(paste('Error at',sname,pos,'! Each individual should have 3 PL score, not',length(gl)))
        print(line[ind])
        break
      }
      #Detecting missing data
      #l1=10^(-gl[1]/10)
      #l2=10^(-gl[2]/10)
      #l3=10^(-gl[3]/10)
      l1=glpow_R(gl[1])
      l2=glpow_R(gl[2])
      l3=glpow_R(gl[3])


      sum=l1+l2+l3
      nl1=l1/sum
      nl2=l2/sum
      nl3=l3/sum
      nl1=sprintf("%.6f",nl1)
      nl2=sprintf("%.6f",nl2)
      nl3=sprintf("%.6f",nl3)
      #nl1=format(round(l1/sum,digits=6),nsmall=6)
      #nl2=format(round(l2/sum,digits=6),nsmall=6)
      #nl3=format(round(l3/sum,digits=6),nsmall=6)
      }
      newline<-paste(nl1,nl2,nl3)
      cat(newline,file = output,sep = ' ')
      cat(' ',file = output)
    }
    cat('\n',file = output)
    #writeLines(newline,con = output,sep = '\n')
    line=unlist(strsplit(readLines(input,n=1),split = '\t'))
  }
  close(input)
  close(output)
}






