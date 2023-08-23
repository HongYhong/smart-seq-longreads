library(rtracklayer)
library(vroom)
library(plyranges)
library(dplyr)

gtf = import("genome/Homo_sapiens.GRCh38.109.chr.gtf")
gtf.gene = gtf %>% filter(type == "gene")
# Y1.bed = import("alignments_hj_human/Y1.sorted.bed")
# Y2.bed = import("alignments_hj_human/Y2.sorted.bed")
# Y3.bed = import("alignments_hj_human/Y3.sorted.bed")

rep1.bed = import("SGNex_bam/SGNex_K562_cDNA_replicate1_run3.bed")
rep2.bed = import("SGNex_bam/SGNex_K562_cDNA_replicate3_run4.bed")

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

# Y1.junction.width = tibble(id = Y1.bed$name,n_junction = lengths(blocks(Y1.bed))-1) %>% 
# left_join(Y1.annotated,by = c("id"))
# Y2.junction.width = tibble(id = Y2.bed$name,n_junction = lengths(blocks(Y2.bed))-1)%>% 
# left_join(Y2.annotated,by = c("id"))
# Y3.junction.width = tibble(id = Y3.bed$name,n_junction = lengths(blocks(Y3.bed))-1)%>% 
# left_join(Y3.annotated,by = c("id"))

rep1.junction.width = tibble(id = rep1.bed$name,n_junction = lengths(blocks(rep1.bed))-1) %>% 
left_join(rep1.annotated,by = c("id"))
rep2.junction.width = tibble(id = rep2.bed$name,n_junction = lengths(blocks(rep2.bed))-1) %>% 
left_join(rep2.annotated,by = c("id"))

# write.table(Y1.junction.width,file = "data/Y1.junction.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y2.junction.width,file = "data/Y2.junction.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y3.junction.width,file = "data/Y3.junction.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)

write.table(rep1.junction.width,file = "control_data/rep1.junction.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(rep2.junction.width,file = "control_data/rep2.junction.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)



