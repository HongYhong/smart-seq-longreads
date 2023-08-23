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

# Y1.introns = psetdiff(Y1.bed,rtracklayer::blocks(Y1.bed))
# Y2.introns = psetdiff(Y2.bed,rtracklayer::blocks(Y2.bed))
# Y3.introns = psetdiff(Y3.bed,rtracklayer::blocks(Y3.bed))

# rep1.introns = psetdiff(rep1.bed,rtracklayer::blocks(rep1.bed))
# rep2.introns = psetdiff(rep2.bed,rtracklayer::blocks(rep2.bed))

# Y1.exons = rtracklayer::blocks(Y1.bed)
# Y2.exons = rtracklayer::blocks(Y2.bed)
# Y3.exons = rtracklayer::blocks(Y3.bed)

rep1.exons = rtracklayer::blocks(rep1.bed)
rep2.exons = rtracklayer::blocks(rep2.bed)

# Y1.junction = tibble(id = rep(Y1.bed$name,lengths(width(Y1.introns))),
# junction_width = unlist(width(Y1.introns))) %>% 
# left_join(Y1.annotated,by = c("id"))
# Y2.junction = tibble(id = rep(Y2.bed$name,lengths(width(Y2.introns))),
# junction_width = unlist(width(Y2.introns)))%>% 
# left_join(Y2.annotated,by = c("id"))
# Y3.junction = tibble(id = rep(Y3.bed$name,lengths(width(Y3.introns))),
# junction_width = unlist(width(Y3.introns)))%>% 
# left_join(Y3.annotated,by = c("id"))

rep1.junction = tibble(id = rep(rep1.bed$name,lengths(width(rep1.exons))),
junction_width = unlist(width(rep1.exons))) %>% 
left_join(rep1.annotated,by = c("id"))

rep2.junction = tibble(id = rep(rep2.bed$name,lengths(width(rep2.exons))),
junction_width = unlist(width(rep2.exons))) %>% 
left_join(rep2.annotated,by = c("id"))

# write.table(Y1.junction,file = "data/Y1.junctionwidth.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y2.junction,file = "data/Y2.junctionwidth.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)
# write.table(Y3.junction,file = "data/Y3.junctionwidth.groupbygene.tsv",sep = "\t",
# quote = F,row.names =F,col.names =T)

write.table(rep1.junction,file = "control_data/rep1.exonwidth.groupbygene.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
write.table(rep2.junction,file = "control_data/rep2.exonwidth.groupbygene.tsv",sep = "\t",
quote = F,row.names =F,col.names =T)
