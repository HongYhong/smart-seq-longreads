stats_mapping_rate = function(mapped){
  require(dplyr)
  mapped = mapped %>% group_by(sample) %>% mutate(prop = paste0(round(value/max(value),2) * 100,"%"))
  p = ggplot(mapped,aes(x = sample,fill = type,y = value)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw()+ 
  scale_fill_manual(values = c("#468B97","#EF6262","#F3AA60"))+
  geom_text(aes(label = prop),position = position_dodge(0.9)
            ,vjust = -1,hjust = 0.5)
  return(p)
}
