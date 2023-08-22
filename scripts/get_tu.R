library(rtracklayer)
library(dplyr)
library(plyranges)

Y1.bed = import("alignments_hj_human/Y1.sorted.bed")
Y2.bed = import("alignments_hj_human/Y2.sorted.bed")
Y3.bed = import("alignments_hj_human/Y3.sorted.bed")

gtf = import("genome/Homo_sapiens.GRCh38.109.chr.gtf")
gtf.gene = gtf %>% filter(type == "gene")

Y1.annotated = tibble(id = Y1.bed$name[unique(queryHits(findOverlaps(Y1.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y2.annotated = tibble(id = Y2.bed$name[unique(queryHits(findOverlaps(Y2.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y3.annotated = tibble(id = Y3.bed$name[unique(queryHits(findOverlaps(Y3.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

Y1.tu = reduce(Y1.bed,ignore.strand = T)
Y2.tu = reduce(Y2.bed,ignore.strand = T)
Y3.tu = reduce(Y3.bed,ignore.strand = T)

Y1.ovlps = findOverlaps(Y1.bed,Y1.tu,ignore.strand = T,type = "within")
Y2.ovlps = findOverlaps(Y2.bed,Y2.tu,ignore.strand = T,type = "within")
Y3.ovlps = findOverlaps(Y3.bed,Y3.tu,ignore.strand = T,type = "within")

Y1.cluster = tibble(read = Y1.bed$name[queryHits(Y1.ovlps)],cluster = subjectHits(Y1.ovlps)) %>% 
group_by(cluster) %>%
mutate(cluster_size = length(unique(read))) %>%
left_join(Y1.annotated,by = c("read" = "id"))

Y2.cluster = tibble(read = Y2.bed$name[queryHits(Y2.ovlps)],cluster = subjectHits(Y2.ovlps)) %>%
group_by(cluster) %>%
mutate(cluster_size = length(unique(read)))%>%
left_join(Y2.annotated,by = c("read" = "id"))

Y3.cluster = tibble(read = Y3.bed$name[queryHits(Y3.ovlps)],cluster = subjectHits(Y3.ovlps)) %>%
group_by(cluster) %>%
mutate(cluster_size = length(unique(read)))%>%
left_join(Y3.annotated,by = c("read" = "id"))

write.table(Y1.cluster,file = "data/Y1.cluster.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y2.cluster,file = "data/Y2.cluster.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(Y3.cluster,file = "data/Y3.cluster.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)




