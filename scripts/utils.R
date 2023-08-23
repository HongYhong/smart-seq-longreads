readBam <- function(bamfile) {
  require(parallel)
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM"),
                                              what = c("qname","flag", "rname", 
                                                       "pos")))
  ops <- GenomicAlignments::CIGAR_OPS
  wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = ops)
  keep.ops <- GenomicAlignments::explodeCigarOps(cigar(bam), ops = ops)
  explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                         unlist(keep.ops)), wdths))
  for (opts in setdiff(GenomicAlignments::CIGAR_OPS, c("=","P","X"))) {
    mcols(bam)[[paste0("nbr", opts)]] <- 
      unlist(mclapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE),mc.cores = 5))
  }

  mcols(bam)$readLength <- rowSums(as.matrix(mcols(bam)[, c("nbrS", "nbrH", "nbrM", "nbrI")]))
  bam
}

makeReadDf <- function(bam) {
  tmp <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::rename(read = qname,
                  nbrJunctions = njunc) %>%
    dplyr::select(-cigar) %>%
    dplyr::mutate(alignedLength = nbrM + nbrI) ## equivalent to readLength-nbrS-nbrH
  return(tmp)
}

get_baseaccuracy = function(df){
  (df$nbrM + df$nbrI + df$nbrD - df$NM)/(df$nbrM + df$nbrI + df$nbrD)
}

nblocks = function(se){
  return(lengths(rowRanges(se)))
}
