library(rtracklayer)
library(vroom)
library(plyranges)
library(dplyr)
source("scripts/utils.R")

Y1.genome = vroom("data/Y1.stat.tsv")
Y2.genome = vroom("data/Y2.stat.tsv")
Y3.genome = vroom("data/Y3.stat.tsv")

gtf = import("genome/Homo_sapiens.GRCh38.109.chr.gtf")
gtf.gene = gtf %>% filter(type == "gene")
Y1.bed = import("alignments_hj_human/Y1.sorted.bed")
Y2.bed = import("alignments_hj_human/Y2.sorted.bed")
Y3.bed = import("alignments_hj_human/Y3.sorted.bed")

Y1.annotated = tibble(id = Y1.bed$name[unique(queryHits(findOverlaps(Y1.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y2.annotated = tibble(id = Y2.bed$name[unique(queryHits(findOverlaps(Y2.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y3.annotated = tibble(id = Y3.bed$name[unique(queryHits(findOverlaps(Y3.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y1.alignedlengthgroupbygene = tibble(id = Y1.genome$read,alignedLength = Y1.genome$alignedLength) %>%
left_join(Y1.annotated,by = c("id"))
Y2.alignedlengthgroupbygene = tibble(id = Y2.genome$read,alignedLength = Y2.genome$alignedLength) %>%
left_join(Y2.annotated,by = c("id"))
Y3.alignedlengthgroupbygene = tibble(id = Y3.genome$read,alignedLength = Y3.genome$alignedLength) %>%
left_join(Y3.annotated,by = c("id"))

write.table(Y1.alignedlengthgroupbygene,file = "data/Y1.alignedlengthbygene.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y2.alignedlengthgroupbygene,file = "data/Y2.alignedlengthbygene.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y3.alignedlengthgroupbygene,file = "data/Y3.alignedlengthbygene.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
