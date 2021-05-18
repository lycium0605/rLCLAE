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
#' @param type haplotype (i.e. male-X) or diploid (Autosome), dipoid as default
#' @return A file in which each line represents a snp. The line looks like chr, pos, -1. -1. -1. (missing data) or gl1, gl2, gl3 transformed from the phred score.


genolik<-function(inputdir,outputdir,type='dip'){
  input=file(inputdir,'r')
  line=unlist(strsplit(readLines(input,n=1),split = '\t'))
  indnum=length(line)-2
  close(input)
  if(type == 'dip'){
    filt1_dip(indnum,inputdir,outputdir)
  }
  else if(type=='hap'){
    filt1_hap(indnum,inputdir,outputdir)
  }
  else{
    print("Please specify a type, dip or hap")
  }

}






