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
  exclude=1000
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
            pop2_dip = refpop2_dip)

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
      anccall(anccall_in,anccall_out,chr,as.character(ind),
              mode_min=mode_min,n_min=n_min,chrlength = chrlen,
              window_size = window,
              delta = delta)

# 06. Generating ancetry tract ----------------------------------------------
      gettract_in<-anccall_out
      gettract_out<-paste0(tmp,"/06_anctract_",chr)
      gettract(
        outputdir = gettract_out,
        datadir = gettract_in,
        indiv_name = as.character(ind),
        chr = chr,
        chrlength = as.character(chrlen),
        value = window,
        exclude = exclude,
        min_n = n_min,
        mode_n = mode_min
      )

    }


  }

}
