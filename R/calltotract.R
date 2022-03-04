#library(data.table)

#' getmode
#'
#' @param v a input
#'
#' @return nothing
# @export
#'
# @examples
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

utils::globalVariables(c(".", "chrom","loc_Dummy", "loc_Minus100", "loc_Plus100",
                         "nxt", "nxt_chrom", "nxt_state",
                         "snp"))

#' @title gettract
#' @description Generating tract from snp level ancestry call.
#' @param datadir input snp level ancestry call file
#' @param outputdir output tract file
#' @param indiv_name name or index to specify the individual
#' @param chr the chromosome it is on
#' @param chrlength the length of the chromosome
#' @param value window size
#' @param mode_n the minimum percentage for a call to be made using majority rule
#' @param min_n the minimum number of SNPs within the window to make a call
#' @param exclude how much of the start/end of each chromosome to ignore
#'
#' @return An output file including tract and majority rule
#' @export
#' @import Rcpp
#' @import data.table
#' @importFrom graphics hist
#' @importFrom stats median
#' @importFrom utils read.delim write.table
# @examples
gettract <- function(datadir,outputdir,
                     indiv_name='50',
                     chr='chrX',chrlength='142711496',
                     value=35000,
                     mode_n=0.5,
                     min_n=20,
                     exclude=50000){

# This function could only run one chromosome at a time
filedir<-datadir
maj_dir<-paste(outputdir,indiv_name,'MajorityRule.txt',sep = '.')
tractdir<-paste(outputdir,indiv_name,'tracts.txt',sep = '.')

# Prepare a data column with snp, call, chrom, indiv_name
fread(filedir,showProgress = T) -> indv_ #chr17 pos call
data.table(indv_)->indv
#indv$chrom <- 'chrX' #add a line full of chrX
indv$chrom<-chr
colnames(indv) <- c('snp', 'call', 'chrom')
indv$indiv <- as.character(indiv_name)
chroms<-cbind(chr,chrlength)
j=1

indv -> tmp
subset(tmp, !(tmp$call == -1)) -> tmp
tmp -> indv
if (nrow(tmp) > 0) {
  # index temporary file of chromosome calls and use setkeys to prepare the data for matching sites within $VALUE of each SNP
  # run match to get sets of sites within `value` of each SNP
  # get the mode of all calls, the total number of calls, and the number of calls which were the mode.
  tmp[,loc_Dummy := snp]; tmp[,.(snp, call)] -> tmp2
  tmp2[,loc_Plus100 := snp + value]; tmp2[,loc_Minus100 := snp - value]
  setkey(tmp,snp,loc_Dummy); setkey(tmp2,loc_Minus100, loc_Plus100)
  #print(nrow(tmp))
  print(paste("Now doing matches for ", chroms[j,1], "....", sep=""))
  Matches <- foverlaps(tmp[,.(snp, loc_Dummy)], tmp2[,.(loc_Minus100,loc_Plus100,call)])
  Matches[,.(n = .N, mode = getmode(call), n_mode=sum(call==getmode(call))), by = .(snp)] -> i1
  rm(Matches); gc(); rm(tmp2)
  print("done with matches!")
  i1$perc <- i1$n_mode/i1$n
  # remove sites where there is not a consensus call by majority rule (at least `mode_n` percent of calls with the same state) or enough nearby ancestry informative sites (`min_n` within `value`)
  subset(tmp, i1$perc >= mode_n & i1$n >= min_n) -> tmp; subset(i1, i1$perc >= mode_n & i1$n >= min_n) -> i1

  # Some format adjustment
  tmp$snp<-as.numeric(as.character(tmp$snp))
  chroms[j,2]<-as.numeric(as.character(chroms[j,2]))

  if (nrow(i1) > 0) {
    # merge majority rule calls for this chromosome to the growing file, removing the first and last X kb (set by "exclude").
    #if (j==1) {cbind(tmp,i1)[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)] -> maj_rule} else {cbind(tmp,i1) -> te; rbind(maj_rule, te[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)]) -> maj_rule}
    cbind(tmp,i1)[as.numeric(tmp$snp) >= exclude & as.numeric(tmp$snp) <= as.numeric(chroms[j,2])-exclude,-c(5:6)] -> maj_rule

    # progress to calling ancestry tracts
    cbind(tmp, i1) -> modes; rm(tmp); rm(i1); if (j > 1) {rm(te)}

    modes$nxt <- c(as.numeric(modes$snp[-1]),max(modes$snp))
    modes$nxt_chrom <- c(as.character(modes$chrom[-1]),modes$chrom[nrow(modes)])
    modes$nxt_state <- c(modes$mode[-1],modes$call[nrow(modes)])

    # keep only sites where this site is different than the next one.
    # when this occurs, the first row is now the first break point
    subset(modes, !(modes$mode == modes$nxt_state & modes$chrom == modes$nxt_chrom)) -> blocks
    blocks[,.(chrom, nxt_chrom, snp, nxt, mode, nxt_state)] -> blocks

    # nxt_state is the one that defines that block
    # add the first block from position 1 to the first break pt
    first <- as.data.frame(t(matrix(c(blocks$chrom[1], blocks$chrom[1], 1, blocks$snp[1], modes$mode[1], blocks$mode[1]))))
    colnames(first) <- colnames(blocks); rbind(first, blocks) -> blocks; rm(first)

    last <- as.data.frame(matrix(c(as.character(chroms[j,1]), as.character(chroms[j,1]), max(modes$snp), chroms[j,2], modes$nxt_state[nrow(modes)], modes$nxt_state[nrow(modes)]), ncol=6)); as.numeric(as.character(last$V3)) -> last$V3; as.numeric(as.character(last$V4)) -> last$V4; as.numeric(as.character(last$V5)) -> last$V5
    ;as.numeric(as.character(last$V6)) -> last$V6;
    colnames(last) <- colnames(blocks); rbind(blocks, last) -> blocks; rm(last)

    blocks$chrom <- blocks$nxt_chrom <- as.character(chroms[j,1]) # If you want the numbers, this should be `j` instead of `chroms$V1[j]`

    # write and read back in to reset formatting
    #write.table(blocks, "./tmp2.INDIV.txt", row.names=F, col.names=T, sep="\t")
    #read.delim("./tmp2.INDIV.txt") -> blocks

    blocks$brk <- (as.numeric(blocks$nxt) + as.numeric(blocks$snp))/2
    blocks$brk[1] <- 1  # start at the first bp of the chromosome
    modes$snp[1] -> blocks$brk[1] # assume the a single ancestry tract up until the first called SNP

    chroms[j,2]  -> blocks$brk[nrow(blocks)] # extend the last tract until the end of the chromosome
    blocks$nxt_brk <- c(as.numeric(blocks$brk[-1]),'end')
    blocks[-nrow(blocks),] -> blocks  # last line no longer matters, we've already extended from the last SNP to the end of the chromosome

    blocks$length <- as.numeric(blocks$nxt_brk) - as.numeric(blocks$brk)
    #blocks$length[1] <- NA; blocks$length[nrow(blocks)-1] <- NA

    # remove uncertainty because we really don't need this now that I'm not doing as much methods testing.
    #uncertainty of the break points to each side
    blocks$u_prev <- (as.numeric(blocks$nxt) - as.numeric(blocks$snp))/2
    blocks$u_prev[1] <- NA
    blocks$u_next <- c(as.numeric(blocks$u_prev[-1]),NA)
    ## blocks u_prev is the number of bases at the start of that tract which were inferred (i.e. before the first AIM with that ancestry call)
    ## blocks u_next is the same thing for the end of that tract

    #Blocks starts at the first SNP and goes to the last one.

    t <- blocks[, c(1, 7, 8, 6, 9:11)]
    colnames(t) <- c("chrom", "start", "end", "state", "length", "prev_inferred", "after_inferred")
    t$name <- indiv_name
    #t$length2 <- as.numeric(t$end) - as.numeric(t$start)
    t$end <- as.numeric(as.character(t$end))
    t$start <- as.numeric(as.character(t$start))

    # fix the extremes of `t`, cropping the first and last `exclude`bp from the chromosome
    t <- subset(t, t$start <= (as.numeric(chroms[j,2])-exclude) & t$end >= exclude)
    t$start[t$start < exclude] <- exclude
    t$end[t$end > (as.numeric(chroms[j,2])-exclude)] <- (as.numeric(chroms[j,2])-exclude)
    t$length[c(1,nrow(t))] <- NA

    # merge with list of tracts per individual
    if (j==1) {t -> tracts} else {rbind(tracts,t) -> tracts}
    if(nrow(tracts) == 1 & length(unique(indv$call)==1)){
      tracts$state<-indv$call[1]
      tracts$length<-tracts$end - tracts$start
      print("Only one tract, using the unique value.")
    }
    rm(t); rm(modes)
    rm(blocks)
  }
}

write.table(maj_rule, maj_dir, row.names=F, col.names=T, quote=F, sep="\t")
write.table(tracts, tractdir, row.names=F, col.names=T, quote=F, sep="\t")
}

