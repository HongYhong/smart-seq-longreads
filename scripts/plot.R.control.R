library(ggplot2)
library(vroom)
library(LSD)
library(ggplotify)
library(patchwork)
library(dplyr)
library(qs)
source("scripts/utils.R")
source("scripts/plot_utils.R")

mapped = read.table("control_data/mapped.tsv",header = T)
mapped_rate.p = stats_mapping_rate(mapped)

qsave(mapped_rate.p,file = "control_data/figures/mapped_rate.p.qs")

rep1.stat = vroom("control_data/rep1.stat.tsv",num_threads = 20)

rep1.stat.tx = vroom("control_data/rep1.transcriptome.tsv",num_threads = 20)

rep1.stat.p = as.ggplot(~heatscatter(rep1.stat$readLength,rep1.stat$alignedLength,log = "xy",
xlab = "Read length",ylab = "Aligned length",main = "rep1 genome"))

rep1.stat.tx.p = as.ggplot(~heatscatter(rep1.stat.tx$readLength,rep1.stat.tx$alignedLength,log = "xy",
xlab = "Read length",ylab = "Aligned length",main = "rep1 transcriptome"))

alignlen.p = (rep1.stat.p)/(rep1.stat.tx.p)
qsave(alignlen.p,file = "control_data/figures/alignlen.p.qs")

rep1.junction = vroom("control_data/rep1.junction.tsv")
rep1.junction$type[is.na(rep1.junction$type)] = "inter-gene"
rep1.junction$sample = "rep1"
rep1.junction$type[which(rep1.junction$type == "annotated")] = "intra-gene"
junction = bind_rows(rep1.junction)

rep1.n_junction.p = rep1.junction %>% unique() %>% filter(!is.na(n_junction)) %>% group_by(type,n_junction) %>% 
summarise(count = n())  %>% ggplot(aes(x = n_junction,y = count,fill = type)) + 
geom_bar(stat="identity",position = position_dodge2(preserve = "single")) + theme_bw() + 
scale_y_continuous(trans = "log1p",breaks = 4^(1:13)) + ggtitle("rep1 Junction number statistics") + 
scale_fill_manual(values = c("#00DFA2","#FF0060")) + scale_color_manual(values = c("#00DFA2","#FF0060"))

njunction.p = rep1.n_junction.p

qsave(njunction.p,file = "control_data/figures/njunction.p.qs")

rep1.junction.width = vroom("control_data/rep1.junctionwidth.groupbygene.tsv")

rep1.junction.width$type[is.na(rep1.junction.width$type)] = "inter-gene"
rep1.junction.width$sample = "rep1"

rep1.junction.width$type[which(rep1.junction.width$type == "annotated")] = "intra-gene"

junction.width = rep1.junction.width

dodge <- position_dodge(width = 0.8)
junction.width.p = junction.width %>% unique() %>% ggplot(aes(x = factor(sample),y = junction_width,fill = factor(type)))+
geom_violin(position=dodge,linewidth=1)+geom_boxplot(position=dodge,width=0.1,outlier.shape = NA,linewidth=1)+ 
theme_bw() + 
scale_y_continuous(trans = "log1p",breaks = 4^(1:13)) + ggtitle("Junction width")+
scale_fill_manual(values = c("#00DFA2","#FF0060"))

qsave(junction.width.p,file = "control_data/figures/junction.width.p.qs")

rep1.junction.width = vroom("control_data/rep1.exonwidth.groupbygene.tsv")

rep1.junction.width$type[is.na(rep1.junction.width$type)] = "inter-gene"
rep1.junction.width$sample = "rep1"

rep1.junction.width$type[which(rep1.junction.width$type == "annotated")] = "intra-gene"

junction.width = rep1.junction.width

dodge <- position_dodge(width = 0.8)
exon.width.p = junction.width %>% unique() %>% ggplot(aes(x = factor(sample),y = junction_width,fill = factor(type)))+
geom_violin(position=dodge,linewidth=1)+geom_boxplot(position=dodge,width=0.1,outlier.shape = NA,linewidth=1)+ 
theme_bw() + 
scale_y_continuous(trans = "log1p",breaks = 4^(1:13)) + ggtitle("Exon width")+
scale_fill_manual(values = c("#00DFA2","#FF0060"))

rep1.baseaccuracy = vroom("control_data/rep1.baseaccuracy.tsv")

rep1.baseaccuracy$type[is.na(rep1.baseaccuracy$type)] = "inter-gene"
rep1.baseaccuracy$sample = "rep1"

rep1.baseaccuracy.p = rep1.baseaccuracy %>% unique()  %>% ggplot(aes(x = sample,y = baseAccuracy,fill = type))+
geom_violin(position=dodge,linewidth=1)+geom_boxplot(position=dodge,width=0.1,outlier.shape = NA,linewidth=1)+ 
theme_bw() + 
scale_y_continuous(trans = "log1p",breaks = 4^(1:13)) + ggtitle("Junction width")+
scale_fill_manual(values = c("#00DFA2","#FF0060"))

qsave(rep1.baseaccuracy.p,file = "control_data/figures/rep1.baseaccuracy.p.qs")

rep1.tu = vroom("control_data/rep1.cluster.tsv")

rep1.tu$type[is.na(rep1.tu$type)] = "inter-gene"
rep1.tu$sample = "rep1"
rep1.tu$type[which(rep1.tu$type == "annotated")] = "intra-gene"

rep1.tu.p = rep1.tu %>% unique()  %>% ggplot(aes(x = cluster_size,color = type))+
geom_density(linewidth = 1) + theme_classic()+scale_x_continuous(trans = "log1p",breaks = 4^(1:8)) + 
ggtitle("Cluster size")+
scale_color_manual(values = c("#00DFA2","#FF0060"))

qsave(rep1.tu.p,file = "control_data/figures/rep1.tu.p.qs")

rep1.algnedbygene = vroom("control_data/rep1.alignedlengthbygene.tsv")

rep1.algnedbygene$type[is.na(rep1.algnedbygene$type)] = "inter-gene"
rep1.algnedbygene$sample = "rep1"
rep1.algnedbygene$type[which(rep1.algnedbygene$type == "annotated")] = "intra-gene"

rep1.algnedbygene.p = rep1.algnedbygene %>% unique()  %>% ggplot(aes(x = alignedLength,color = type))+
geom_density(alpha=0.5,size = 1) + theme_classic()+scale_x_continuous(trans = "log1p",breaks = 4^(1:8)) + 
ggtitle("Aligned length") + scale_color_manual(values = c("#00DFA2","#FF0060"))

qsave(rep1.algnedbygene.p,file = "control_data/figures/rep1.algnedbygene.p.qs")
