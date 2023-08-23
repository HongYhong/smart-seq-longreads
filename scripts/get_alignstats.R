library(GenomicAlignments)
library(parallel)
library(stringr)
source("scripts/utils.R")

# bamfiles = list.files("alignments_hj_human",pattern = "*.bam$",full.names = T)
# bamfiles = list.files("alignment_hj_transcriptome",pattern = "*.bam$",full.names = T)
# bamfiles = list.files("SGNex_bam",pattern = "*.bam$",full.names = T)
bamfiles = c("SGNex_transcriptome/SGNex_K562_cDNA_replicate1_run3.sorted.bam")

bamdfs = lapply(bamfiles,function(bamfile){
  bam = readBam(bamfile)
  bamdf = makeReadDf(bam)
})

for (i in 1:length(bamfiles)){
  outfile = str_replace(basename(bamfiles[[i]]),".bam",".stat.transcriptome.tsv")
  write.table(bamdfs[[i]],file = paste0("control_data/",outfile),sep = "\t",quote = F,row.names =F,col.names =T)
}
