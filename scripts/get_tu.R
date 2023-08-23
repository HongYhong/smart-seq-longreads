library(rtracklayer)
library(dplyr)
library(plyranges)

# Y1.bed = import("alignments_hj_human/Y1.sorted.bed")
# Y2.bed = import("alignments_hj_human/Y2.sorted.bed")
# Y3.bed = import("alignments_hj_human/Y3.sorted.bed")

rep1.bed = import("SGNex_bam/SGNex_K562_cDNA_replicate1_run3.bed")
rep2.bed = import("SGNex_bam/SGNex_K562_cDNA_replicate3_run4.bed")

gtf = import("genome/Homo_sapiens.GRCh38.109.chr.gtf")
gtf.gene = gtf %>% filter(type == "gene")

# Y1.annotated = tibble(id = Y1.bed$name[unique(queryHits(findOverlaps(Y1.bed,gtf.gene,ignore.strand = T)))],
# type = "annotated")

# Y2.annotated = tibble(id = Y2.bed$name[unique(queryHits(findOverlaps(Y2.bed,gtf.gene,ignore.strand = T)))],
# type = "annotated")

# Y3.annotated = tibble(id = Y3.bed$name[unique(queryHits(findOverlaps(Y3.bed,gtf.gene,ignore.strand = T)))],
# type = "annotated")

rep1.annotated = tibble(id = rep1.bed$name[unique(queryHits(findOverlaps(rep1.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")
rep2.annotated = tibble(id = rep2.bed$name[unique(queryHits(findOverlaps(rep2.bed,gtf.gene,ignore.strand = T)))],
type = "annotated")

# Y1.tu = reduce(Y1.bed,ignore.strand = T)
# Y2.tu = reduce(Y2.bed,ignore.strand = T)
# Y3.tu = reduce(Y3.bed,ignore.strand = T)

rep1.tu = reduce(rep1.bed,ignore.strand = T)
rep2.tu = reduce(rep2.bed,ignore.strand = T)

# Y1.ovlps = findOverlaps(Y1.bed,Y1.tu,ignore.strand = T,type = "within")
# Y2.ovlps = findOverlaps(Y2.bed,Y2.tu,ignore.strand = T,type = "within")
# Y3.ovlps = findOverlaps(Y3.bed,Y3.tu,ignore.strand = T,type = "within")

rep1.ovlps = findOverlaps(rep1.bed,rep1.tu,ignore.strand = T,type = "within")
rep2.ovlps = findOverlaps(rep2.bed,rep2.tu,ignore.strand = T,type = "within")

# Y1.cluster = tibble(read = Y1.bed$name[queryHits(Y1.ovlps)],cluster = subjectHits(Y1.ovlps)) %>% 
# group_by(cluster) %>%
# mutate(cluster_size = length(unique(read))) %>%
# left_join(Y1.annotated,by = c("read" = "id"))

# Y2.cluster = tibble(read = Y2.bed$name[queryHits(Y2.ovlps)],cluster = subjectHits(Y2.ovlps)) %>%
# group_by(cluster) %>%
# mutate(cluster_size = length(unique(read)))%>%
# left_join(Y2.annotated,by = c("read" = "id"))

# Y3.cluster = tibble(read = Y3.bed$name[queryHits(Y3.ovlps)],cluster = subjectHits(Y3.ovlps)) %>%
# group_by(cluster) %>%
# mutate(cluster_size = length(unique(read)))%>%
# left_join(Y3.annotated,by = c("read" = "id"))

rep1.cluster = tibble(read = rep1.bed$name[queryHits(rep1.ovlps)],cluster = subjectHits(rep1.ovlps)) %>% 
group_by(cluster) %>%
mutate(cluster_size = length(unique(read))) %>%
left_join(rep1.annotated,by = c("read" = "id"))

rep2.cluster = tibble(read = rep2.bed$name[queryHits(rep2.ovlps)],cluster = subjectHits(rep2.ovlps)) %>% 
group_by(cluster) %>%
mutate(cluster_size = length(unique(read))) %>%
left_join(rep2.annotated,by = c("read" = "id"))

# write.table(Y1.cluster,file = "data/Y1.cluster.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y2.cluster,file = "data/Y2.cluster.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y3.cluster,file = "data/Y3.cluster.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)

write.table(rep1.cluster,file = "control_data/rep1.cluster.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(rep2.cluster,file = "control_data/rep2.cluster.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)




