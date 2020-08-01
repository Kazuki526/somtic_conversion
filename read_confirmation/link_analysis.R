library(tidyverse)
library(pipeR)
library(ggsignif)
loadNamespace('cowplot')
setwd('/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/')


write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}

all_tbl = read_tsv("hetero_germ_linked_candidate.tsv")
purity_cutoff=0.7
tvaf_cutoff=0.8

confirmed=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% (?.%>>%count(sample_id,chr,start,ref,alt)%>>%count()) %>>%
  filter(`alt-ref` >0 & `alt-alt` >0)%>>%(?.%>>%count(sample_id,chr,start,ref,alt)%>>%count())
  #dplyr::select(sample_id,chr,start,ref,alt,purity,germ_variant,distance,linked_read,`ref-ref`,`ref-alt`,`alt-ref`,`alt-alt`,tVAF)%>>%
  #arrange(sample_id,chr,start)%>>%
  #filter(!(`alt-ref` >0 & `alt-alt` >0))%>>%View

all_tbl%>>%inner_join(confirmed%>>%count(sample_id,chr,start,ref,alt))%>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter((`alt-ref` + `alt-alt`>5) & !(`alt-ref` >0 & `alt-alt` >0))%>>%View

hist=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
#  mutate(linked_convert=ifelse(`alt-ref` >4 & `alt-alt` >4, "confirmed","not confirmed"))%>>%
#  dplyr::select(-variant_type,-mutect_dcv_posi,-mutect_mut_num,-linked_read,-screening)%>>%
#  group_by(sample_id,chr,start)%>>%filter(any(linked_convert=="confirmed"))%>>%
#  write_df("~/Dropbox/work/somatic_gene_conversion/confirmed_gene_conversion.tsv")
  ggplot()+geom_histogram(aes(x=distance,fill=linked_convert),color="black",binwidth = 10) +
  theme_bw()+xlab("Distance to candidate mutation")+ylab("Count")
  

ratio=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(distance_class=ifelse(distance<=100,"1 - 100",ifelse(distance<=200,"101 - 200","301 -")))%>>%
  #mutate(distance_class=factor(distance_class,levels = c("1 - 100","101 - 200","301 -")))%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
  (?.%>>%count(distance_class,linked_convert))%>>%
  ggplot(aes(x=distance_class,fill=linked_convert))+geom_bar(color="black",position = "fill",width = 0.95)+
#  geom_signif(stat = "identity",aes(x=2,xend=3,y=1.1,yend=1.1,annotation="**"),textsize = 5)+
#  geom_signif(data=data.frame(x=c(1.05,2.05),xend=c(1.95,2.95), y=c(1.1,1.1),yend=c(1.1,1.1),ann=c("NS","**")),
#              stat="identity",aes(x=x,xend=xend, y=y,yend=yend,annotation=ann))+
#  scale_y_continuous(breaks = c(0,0.5,1),limits = c(0,1.2),expand = c(0,0))+
  xlab("Distance")+ylab("Ratio")+
  theme_bw()+
  theme(axis.ticks =element_blank(),legend.title = element_blank(),
        axis.title = element_text(size=24),legend.text = element_text(size=15),
        legend.position = "top",legend.justification = "center",
        axis.text =element_text(size=15),panel.grid.major.x = element_blank())
ratio
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_distance.pdf",ratio,height = 6,width = 4)
#test 1-100 vs 101-200 : p-value = 0.1037948
fisher.test(matrix(c(549,3419,224,1212),nrow = 2))
#test 101-200 vs 301- : p-value = 0.009443576
fisher.test(matrix(c(224,1212,107,411),nrow = 2))

cowplot::plot_grid(hist+theme(legend.position ="none" ),NULL,ratio,nrow = 1,rel_widths = c(10,1,10))
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_distance.pdf",height = 3,width = 10)

all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
  group_by(sample_id,chr,start,ref,alt,tVAF)%>>%
  summarise(linked_convert=ifelse(any(linked_convert=="confirmed"),"confirmed","not confirmed"))%>>%
  ungroup()%>>%(?.%>>%count(linked_convert))%>>%
  mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  ggplot()+geom_histogram(aes(x=tVAF,fill=linked_convert),color="black",position = "fill",bins=25)+
  scale_y_continuous(breaks = c(0,0.5,1),limits = c(0,1),expand = c(0,0))+
  scale_x_continuous(breaks = c(0.75,0.8,0.85,0.9,0.95,1),expand=c(0,0))+
  ylab("Ratio")+theme_bw()+
  theme(axis.ticks =element_blank(),legend.title = element_blank(),
        axis.text.x =element_text(size=10),panel.grid.major.x = element_blank())
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_tVAF.pdf",height = 3,width = 6)
