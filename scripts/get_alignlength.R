library(GenomicAlignments)
library(parallel)
library(stringr)
source("scripts/utils.R")

bamfiles = list.files("alignments_hj_human",pattern = "*.bam$",full.names = T)

bamdfs = mclapply(bamfiles,function(bamfile){
  bam = readBam(bamfile)
  bamdf = makeReadDf(bam)
  tibble(read = bamdf$qname,readLength = bamdf$readLength,alignedLength = bamdf$alignedLength)
},mc.cores =10)

for (i in 1:length(bamfiles)){
  outfile = str_replace(basename(bamfiles[[i]]),".sorted.bam",".lengthstat.tsv")
  write.table(bamdfs[[i]],file = paste0("data/",outfile),sep = "\t",quote = F,row.names =F,col.names =T)
}
