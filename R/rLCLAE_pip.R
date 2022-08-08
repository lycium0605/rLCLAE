#' rLCLAE
#' @description The summarized pipeline for rLCLAE
#' @param input vcf input
#' @param tmp intermediate file storage
#' @param output final output (not yet used)
#' @param chrlist a data frame for chr name & length, from fai
#' @param refpop1_dip numeric vector for index of ref ind for pop1, need expand
#' @param refpop2_dip numeric vector for index of ref ind for pop2, need expand
#' @param test_ind numeric vector for index of test individual, default all
#' @param inputtype hap or dip
#' @param freqtype hap, dip or merge
#' @param window window size
#' @param mode_min mode min
#' @param n_min n min in window
#' @param delta freq diff lower cutoff
#' @param exclude exclude length
#' @param inputdir_hap missing, do not touch, not finished for ancfreq
#' @param refpop1_hap missing, do not touch, not finished for ancfreq
#' @param refpop2_hap missing, do not touch, not finished for ancfreq
#' @param mergetype union, do not touch, not finished for ancfreq
#' @param zero_value a small zero value added to the zeros, if zero then throw out sites with any zero
#'
#' @export
#'
#'
rLCLAE<-function(
  input, # input vcf directory
  tmp, # a folder for temporary files
  output, #output directory
  chrlist, #a data frame with two columns, the first is chrname, second is length in bp
  refpop1_dip,
  refpop2_dip,
  test_ind="all",
  inputtype = "dip", #input type
  freqtype = "dip", #ancfreq type, could be hap or merge as well
  window=50000,
  mode_min=0.3,
  n_min=10,
  delta=0.1,
  exclude=1000,
  #ancfreq input
  inputdir_hap = "missing",
  refpop1_hap = "missing",
  refpop2_hap = "missing",
  mergetype = "union",
  zero_value = 10^(-7)
){

  # Check if the directory end with /
  if(substr(tmp,nchar(tmp),nchar(tmp))=='/'){
    tmp<-gsub('.{1}$','',tmp)
  }

# 01. Preprocess input vcf ---------------------------------------
  message(paste("Pre-processing file from",input))
  preout<-paste(tmp,"01_preprocessed",sep = '/')
  preprocess(input,preout)


# 02.  Calculate genotype likelihood --------------------------------------
  message("Starting calculation of genotype likelihood...")
  genoout<-paste(tmp,"02_genolik",sep = '/')
  genoout_suf<-genolik(inputdir=preout,
          outputdir = genoout,
          type = inputtype) # Get a list of suffix (chromosome names)
  message("Genotype likelihood stored in:")
  print(paste0(genoout,"_",genoout_suf))


# Split by chromosome if multiple were present ----------------------------
  for(chr in genoout_suf){

## 03. Calculate ancestral allele frequency -------------------------------------
    ancfreq_in<-paste0(genoout,"_",chr)
    ancfreq_out<-paste0(tmp,"/03_ancfreq_",chr)
    ancfreq(outputdir = ancfreq_out,
            inputdir_dip = ancfreq_in,
            pop1_dip = refpop1_dip,
            pop2_dip = refpop2_dip,
            inputdir_hap = inputdir_hap,
            pop1_hap = refpop1_hap,
            pop2_hap = refpop2_hap,
            mergetype = mergetype)

# 04. Calculate ancestral likelihood --------------------------------------
    anclik_in<-ancfreq_in
    ancfreq_in_type<-paste(ancfreq_out,freqtype,sep="_")
    anclik_out<-paste0(tmp,"/04_anclik_",chr)
    ind_suffix<-anclik(genodir = anclik_in,
           ancfreqdir = ancfreq_in_type,
           outputdir = anclik_out,
           type = freqtype,
           test = test_ind)
# For each individual -----------------------------------------------------
    # Get chromosome length
    chrlen<-as.numeric(chrlist[which(chrlist[,1]==chr),2])
    #print(chrlen)
    for(ind in ind_suffix){

### 05. Generating ancestral call -------------------------------------------
      anccall_in<-paste0(tmp,"/04_anclik_",chr,"_",ind)
      anccall_out<-paste0(tmp,"/05_anccall_",chr,"_",ind)
      anccall(inputdir=anccall_in,
              outputdir=anccall_out,
              chrom_name=chr,
              indiv_name=as.character(ind),
              mode_min=mode_min,
              n_min=n_min,
              chrlength = chrlen,
              window_size = window,
              delta = delta,
              ploidy = inputtype,
              zero_value = zero_value)
      # fixed parameters: round1=1, round2=1, round3=2, auto_optimize=F, zero_value=10^-7

# 06. Generating ancetry tract ----------------------------------------------
      gettract_in<-anccall_out
      gettract_out<-paste0(tmp,"/06_anctract_",chr)
      gettract(
        outputdir = gettract_out,
        datadir = gettract_in,
        indiv_name = ind,
        chr = chr,
        chrlength = as.character(chrlen),
        value = window,
        exclude = as.numeric(exclude),
        min_n = as.numeric(n_min),
        mode_n = as.numeric(mode_min)
      )
    }


  }

}
