library(tidyverse)
library(pipeR)
library(gridExtra)
setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/")
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

patient_list = read_tsv("patient_list.tsv")
all_maf = read_tsv("all_pass_with_ascat.maf.gz")
purity = read_tsv("/Volumes/areca42TB2/gdc/purity/Dvir_purity_data.tsv")
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role=`Role in Cancer`)%>>%dplyr::select(gene,role)

patient_list = patient_list %>>%
  left_join(purity,by=c("sample_id","cancer_type","CPE"))%>>%
  mutate(ALL_MAX = pmax(ASCAT,ABSOLUTE,HE_staining,ESTIMATE,LUMP,na.rm = T))

puri_tvaf_lm = all_maf%>>%
  filter(allele_num!=0)%>>%
  filter_at(c("ASCAT","HE_staining","CPE"),all_vars(!is.na(.)))%>>%
  group_by(patient_id,sample_id)%>>%filter(n()>10)%>>%ungroup()%>>%
  tidyr::gather(key=purity_class,value=purity,ASCAT,HE_staining,CPE)%>>%
  mutate(tVAF=t_alt/(t_depth*(allele_num*purity/(allele_num*purity+2*(1-purity)))))%>>%mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  group_by(purity_class,purity,patient_id,sample_id)%>>%summarise(`mean tVAF`=mean(tVAF))%>>%ungroup()%>>%
  group_by(purity_class)%>>%summarise(cor=cor(purity,`mean tVAF`))
all_maf%>>%filter(allele_num!=0)%>>%
  filter_at(c("ASCAT","HE_staining","CPE"),all_vars(!is.na(.)))%>>%
  group_by(patient_id,sample_id)%>>%filter(n()>10)%>>%ungroup()%>>%
  tidyr::gather(key=purity_class,value=purity,ASCAT,HE_staining,CPE)%>>%
  mutate(tVAF=t_alt/(t_depth*(allele_num*purity/(allele_num*purity+2*(1-purity)))))%>>%mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  group_by(purity_class,purity,patient_id,sample_id)%>>%summarise(`mean tVAF`=mean(tVAF))%>>%ungroup()%>>%
  (?.%>>%count(purity_class))%>>%
  ggplot(aes(x=purity,y=`mean tVAF`))+
  geom_point()+
  facet_wrap(~ purity_class)+
  geom_text(data = puri_tvaf_lm,aes(x=0.3,y=0.2,label=paste0("r = ",signif(cor,4))))+
  stat_smooth(method = "lm",colour="red")+
  theme_bw()
ggsave("~/Dropbox/work/somatic_gene_conversion/draft/figs/purity_tvaf_cor.pdf")
puri05_tvaf_lm = all_maf%>>%
  filter(allele_num!=0)%>>%
  filter_at(c("ASCAT","HE_staining","CPE"),all_vars(!is.na(.)))%>>%
  group_by(patient_id,sample_id)%>>%filter(n()>10)%>>%ungroup()%>>%
  tidyr::gather(key=purity_class,value=purity,ASCAT,HE_staining,CPE)%>>%
  mutate(tVAF=t_alt/(t_depth*(allele_num*purity/(allele_num*purity+2*(1-purity)))))%>>%mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  group_by(purity_class,purity,patient_id,sample_id)%>>%summarise(`mean tVAF`=mean(tVAF))%>>%ungroup()%>>%
  filter(purity>0.5)%>>%
  group_by(purity_class)%>>%summarise(cor=cor(purity,`mean tVAF`))
all_maf%>>%filter(allele_num!=0)%>>%
  filter_at(c("ASCAT","HE_staining","CPE"),all_vars(!is.na(.)))%>>%
  group_by(patient_id,sample_id)%>>%filter(n()>10)%>>%ungroup()%>>%
  tidyr::gather(key=purity_class,value=purity,ASCAT,HE_staining,CPE)%>>%
  mutate(tVAF=t_alt/(t_depth*(allele_num*purity/(allele_num*purity+2*(1-purity)))))%>>%mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  group_by(purity_class,purity,patient_id,sample_id)%>>%summarise(`mean tVAF`=mean(tVAF))%>>%ungroup()%>>%
  filter(purity>0.5)%>>%  (?.%>>%count(purity_class))%>>%
  ggplot(aes(x=purity,y=`mean tVAF`))+
  geom_point()+
  facet_wrap(~ purity_class)+
  geom_text(data = puri05_tvaf_lm,aes(x=0.55,y=0.2,label=paste0("r = ",signif(cor,4))),size=5)+
  stat_smooth(method = "lm",colour="red")+
  theme_bw()
ggsave("~/Dropbox/work/somatic_gene_conversion/draft/figs/purity05_tvaf_cor.pdf",width = 18,height = 8)

all_maf%>>%filter(allele_num!=0)%>>%
  mutate(purity=ifelse(!is.na(CPE),CPE,ifelse(!is.na(HE_staining),HE_staining,ASCAT)))%>>%
  filter(purity>0.5)%>>%
  mutate(tVAF=t_alt/(t_depth*(allele_num*purity/(allele_num*purity+2*(1-purity)))))%>>%mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  group_by(purity,patient_id,sample_id)%>>%summarise(`mean tVAF`=mean(tVAF))%>>%ungroup()%>>%(?.)%>>%
  (?.%>>%{cor(.$purity,.$`mean tVAF`)})%>>%
  ggplot(aes(x=purity,y=`mean tVAF`))+
  geom_point()+
  geom_text(data=~.%>>%summarise(cor=cor(purity,`mean tVAF`)),
    aes(x=0.55,y=0.15,label=paste0("r = ",signif(cor,4))),size=10)+
  stat_smooth(method = "lm",colour="red")+
  theme_bw()
ggsave("~/Dropbox/work/somatic_gene_conversion/draft/figs/allv_purity05_tvaf_cor.pdf")








if(0){
all_maf%>>%
  filter_at(c("ASCAT","HE_staining","CPE","MAX"),all_vars(!is.na(.)))%>>%
  tidyr::gather(purity_class,purity,ASCAT,HE_staining,CPE,MAX)%>>%
  mutate(vaf=t_alt/(t_depth*purity))%>>%
  filter(vaf>0.75)%>>%
  count(sample_id,purity_class)%>>%dplyr::rename(tumor_sample_id=sample_id)%>>%
  right_join(patient_list%>>%
               filter_at(c("ASCAT","HE_staining","CPE","MAX"),all_vars(!is.na(.)&.!=0))%>>%
               tidyr::gather(purity_class,purity,ASCAT,HE_staining,CPE,MAX))%>>%
  mutate(n=ifelse(is.na(n),0,n))%>>%mutate(ratio=n/mutation_num)%>>%
  filter(purity>0.5)%>>% #change this value
  ggplot(aes(x=purity,y=n/mutation_num))+
  geom_point()+
  facet_wrap(~ purity_class)+
  stat_smooth(method = "lm",colour="red")+
  theme_bw()
}