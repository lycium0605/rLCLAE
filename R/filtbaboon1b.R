#Step 0. filter your vcf with the linux command that looks like
#cat recode.vcf | grep ^20 | cut -f 1,2,10- | sed -e 's/:/ /g' -e 's/\./999/g' | sed -e 's/\// /g' > test_clean.vcf
#filtbaboon1b.c - Takes filtered vcf files (e.g., in1)
#and converts them into a list of genotype likelihood (PL) files.

#Inputfile Outputfile individual_num

#Output format: Each line is a SNP, consisting of the base position
#followed by the 3 genotype likelihoods (or -1. for missing data) for each
#of n different individuals.  n=47 for the example data.

## usethis namespace: start
#' @useDynLib rLCLAE, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#inputdir='/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test_clean.vcf'
#outputdir='/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test.genolik_R'
#options(scipen = 200)

#Roxygen commands
#' Calculating the relative genotype likelihood from a vcf file.
#'
#' @param inputdir The input directory of the cleaned vcf file.
#' @param outputdir The output directory of the genolik file.
#' @param type haploid (i.e. male-X) or diploid (Autosome), dipoid as default
#' @return A file in which each line represents a snp. The line looks like chr, pos, -1. -1. -1. (missing data) or gl1, gl2, gl3 transformed from the phred score.
#' @export

genolik<-function(inputdir,outputdir,type='dip'){
  root=outputdir
  skip_line=0

  input=file(inputdir,'r')
  line=unlist(strsplit(readLines(input,n=1),split = '\t'))
  indnum=length(line)-2
  close(input)

  x<-get_chrlist(inputdir)
  chr<-x[seq(2,length(x),2)]
  rep<-as.numeric(x[seq(1,length(x),2)])

  print(paste("Finding",length(chr),"unique chromosomes in the input vcf."))

  for(i in 1:length(chr)){
    print(paste("Calculating genotype likelihood for chromosome",chr[i]))
    root_chr=paste(root,"_",chr[i],sep = "")
    #print(root_chr)
    read_line=as.integer(rep[i])
    skip_line=as.integer(skip_line)
    #print(data.class(read_line))
    if(type == 'dip'){
      filt1_dip(indnum,skip_line,read_line,inputdir,root_chr)
    }
    else if(type=='hap'){
      filt1_hap(indnum,skip_line,read_line,inputdir,root_chr)
    }
    else{
      print("Please specify a type, dip or hap")
    }
    skip_line=skip_line+read_line
    #print(skip_line)
  }
}


#' get_chrlist
#'
#' @param vcf the input vcf file
#'
#' @return a character vector rep,chr
#' @export
#'
# @examples
get_chrlist<-function(vcf="/Users/lycium/Desktop/Jennylab/rpackage_LCLAE/rawtestdata/test_clean.vcf"){
  c=paste("echo $(cat ",vcf, " | cut -f 1 | uniq -c)",sep = '')
  system(c,intern = TRUE) -> out
  out<-strsplit(unlist(out),split = ' ')
  unlist(out)->x
  return (x)
}




