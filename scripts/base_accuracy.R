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

Y1.baseaccuracy = tibble(id = Y1.genome$read,baseAccuracy = 
(Y1.genome$nbrM+Y1.genome$nbrI+Y1.genome$nbrD-Y1.genome$NM)/
(Y1.genome$nbrM+Y1.genome$nbrI+Y1.genome$nbrD)) %>% 
left_join(Y1.annotated,by = c("id"))

Y2.baseaccuracy = tibble(id = Y2.genome$read,baseAccuracy = 
(Y2.genome$nbrM+Y2.genome$nbrI+Y2.genome$nbrD-Y2.genome$NM)/
(Y2.genome$nbrM+Y2.genome$nbrI+Y2.genome$nbrD)) %>% 
left_join(Y2.annotated,by = c("id"))

Y3.baseaccuracy = tibble(id = Y3.genome$read,baseAccuracy = 
(Y3.genome$nbrM+Y3.genome$nbrI+Y3.genome$nbrD-Y3.genome$NM)/
(Y3.genome$nbrM+Y3.genome$nbrI+Y3.genome$nbrD)) %>% 
left_join(Y3.annotated,by = c("id"))

write.table(Y1.baseaccuracy,file = "data/Y1.baseaccuracy.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y2.baseaccuracy,file = "data/Y2.baseaccuracy.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y3.baseaccuracy,file = "data/Y3.baseaccuracy.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
