
#' anccall_check_int
#'
#' @param inputdir The file to the ancestral likelihood file
#' @param outputdir The file to store the ancestry call
#' @param chrom_name The name of the chromosome
#' @param indiv_name The name of the individual
#' @param mode_min The minimum percentage for the mode to make a call, default set to 0.5
#' @param n_min the minimum number of SNPs within the window to make a call, default set to zero
#' @param delta The cutoff for ancestral allele frequency difference
#' @param window_size The size of the sliding window
#' @return nothing
#' @export
#'
# @examples
anccall_check_int<-function(inputdir,outputdir,
                  chrom_name,indiv_name,
                  mode_min=0.5,n_min=0,
                  delta=0.2, window_size=50000){
  flag=0
  if(!file.exists(inputdir)){
    print("The ancestral likelihood file you provided does not exist. Please check.")
    flag=1
    }
  if(flag==0){
    #Calculating snp number in the input
    system(paste("echo $(wc -l",inputdir,")"), intern = TRUE)->lines
    nline<-unlist(strsplit(lines,split = ' '))
    smax=as.numeric(nline[1])
    print(paste("Set SMAX to",smax))
    #得到一个包含snp,call,chrom,indiv,n,mode,n_mode,perc的表格maj_rule
    #Call ancestry
    # anccall_c_test(deltaf=delta,window=window_size,
    #           SMAX=smax,
    #           anclikdir=inputdir,output=outputdir,
    #           chrom=chrom_name,indiv=indiv_name,
    #           mode=mode_min,n=n_min)
    # anccall_c_nozero_int(deltaf=delta,window=window_size,
    #           SMAX=smax,
    #           anclikdir=inputdir,
    #           int1=paste0(outputdir,"_nozero_int1"),
    #           int2=paste0(outputdir,"_nozero_int2"),
    #           output=paste0(outputdir,"_nozero"),
    #           chrom=chrom_name,indiv=indiv_name,
    #           mode=mode_min,n=n_min)
    anccall_c_int(deltaf=delta,window=window_size,
                         SMAX=smax,
                         anclikdir=inputdir,
                         int1=paste0(outputdir,"_int1"),
                         int2=paste0(outputdir,"_int2"),
                         output=outputdir,
                         chrom=chrom_name,indiv=indiv_name,
                         mode=mode_min,n=n_min)
  }
  }
