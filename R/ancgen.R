#' ancgen
#' @description A general wrapped-up function for steps after ancestral frequency estimation.
#' @param genodir The directory to genolik file.
#' @param ancfreqdir The directory to ancestral frequency file.
#' @param outputdir The folder to save output files.
#' @param type 'dip' for diploid data (default) and 'hap' for haploid data.
#' @param test The individuals to be tested. 'all' for all individuals (default), or a numeric vector containing the individual number to be estimated
#' @param delta The cutoff for ancestral allele frequency difference
#' @param window The size of the sliding window for round 1 and 2, and half of the sliding window for round 3.
#' @param chrlength A list of the length to the chromosome.
#' @param mode_n The minimum percentage for the mode to make a call, default set to 0.5.
#' @param min_n The minimum number of SNPs within the window to make a call, default set to zero.
#' @param exclude The length of the start/end of each chromosome to ignore.
#'
#' @return
#' @export
ancgen<-function(genodir,ancfreqdir,outputdir,
                 type = "dip", test = "all", #anclik parameter
                 delta, window, #anccall parameter
                 chrlength, mode_n = 0.5, min_n = 20, exclude = 50000 # Get tract parameter
){
  a<-1;
}
