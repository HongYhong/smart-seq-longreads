library(GenomicAlignments)
library(parallel)
library(stringr)
source("scripts/utils.R")

# bamfiles = list.files("alignments_hj_human",pattern = "*.bam$",full.names = T)
# bamfiles = list.files("alignment_hj_transcriptome",pattern = "*.bam$",full.names = T)
bamfiles = list.files("SGNex_bam",pattern = "*.bam$",full.names = T)

bamdfs = mclapply(bamfiles,function(bamfile){
  bam = readBam(bamfile)
  bamdf = makeReadDf(bam)
},mc.cores =10)

for (i in 1:length(bamfiles)){
  outfile = str_replace(basename(bamfiles[[i]]),".bam",".stat.tsv")
  write.table(bamdfs[[i]],file = paste0("data/",outfile),sep = "\t",quote = F,row.names =F,col.names =T)
}
