library(qs)
library(bambu)
library(rtracklayer)
library(plyranges)
library(dplyr)
source("scripts/utils.R")

gtf = import("genome/Homo_sapiens.GRCh38.109.chr.gtf")
gtf.gene = gtf %>% filter(type == "gene")

TE = import("genome/GRCh38_Ensembl_rmsk_TE.gtf.gz")
TE.tx = asBED(split(TE,TE$transcript_id))

# se1 = qread("data/bambu/se.NDR_1.qs")
# se1.gene = transcriptToGeneExpression(se1)
# se1.gene.filtered = se1.gene[rowSums(assay(se1.gene))>=1,]
# se1.gene.filtered.gr = asBED(rowRanges(se1.gene.filtered))

# se1.geneovlps = tibble(id = se1.gene.filtered.gr$GENEID[unique(queryHits(findOverlaps(se1.gene.filtered.gr,
#                   gtf.gene)))],type = "Intra-gene")

# se1.teovlps = tibble(id = se1.gene.filtered.gr$GENEID[unique(queryHits(findOverlaps(se1.gene.filtered.gr,
#                   TE.tx)))],TE = "TE")

# se1.geneovlps = tibble(id = se1.gene.filtered.gr$GENEID,genomic_len = width(se1.gene.filtered.gr),rna_len = 
# sum(width(rowRanges(se1.gene.filtered))),count = rowSums(assay(se1.gene.filtered)),
# njunction = nblocks(se1.gene.filtered)-1,nexon = nblocks(se1.gene.filtered)
# ) %>% 
# dplyr::left_join(se1.geneovlps,by = c("id")) %>% dplyr::left_join(se1.teovlps,by = c("id"))
# se1.geneovlps$type[is.na(se1.geneovlps$type)] = "Inter-gene"
# se1.geneovlps$TE[is.na(se1.geneovlps$TE)] = "Non-TE"

# se1.geneovlps$newclass = ifelse(grepl("Bambu",se1.geneovlps$id),"novel","known")

# write.table(se1.geneovlps,"data/bambu/transcript_ovlps.tsv",row.names = F,col.names = T,quote = F,sep = "\t")

se2 = qread("control_data/bambu/se.NDR_1.qs")
se2.gene = transcriptToGeneExpression(se2)
se2.gene.filtered = se2.gene[rowSums(assay(se2.gene))>=1,]
se2.gene.filtered = se2.gene.filtered[which(!any(strand(rowRanges(se2.gene.filtered)) == "*")),]
se2.gene.filtered.gr = asBED(rowRanges(se2.gene.filtered))

se2.geneovlps = tibble(id = se2.gene.filtered.gr$GENEID[unique(queryHits(findOverlaps(se2.gene.filtered.gr,
                  gtf.gene)))],type = "Intra-gene")

se2.teovlps = tibble(id = se2.gene.filtered.gr$GENEID[unique(queryHits(findOverlaps(se2.gene.filtered.gr,
                  TE.tx)))],TE = "TE")

se2.geneovlps = tibble(id = se2.gene.filtered.gr$GENEID,genomic_len = width(se2.gene.filtered.gr),rna_len = 
sum(width(rowRanges(se2.gene.filtered))),count = rowSums(assay(se2.gene.filtered)),
njunction = nblocks(se2.gene.filtered)-1,nexon = nblocks(se2.gene.filtered)
) %>%
left_join(se2.geneovlps,by = c("id")) %>% left_join(se2.teovlps,by = c("id"))
se2.geneovlps$type[is.na(se2.geneovlps$type)] = "Inter-gene"
se2.geneovlps$TE[is.na(se2.geneovlps$TE)] = "Non-TE"

se2.geneovlps$newclass = ifelse(grepl("Bambu",se2.geneovlps$id),"novel","known")

write.table(se2.geneovlps,"control_data/bambu/transcript_ovlps.tsv",row.names = F,col.names = T,quote = F,sep = "\t")

