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

patient_list = read_tsv("patient_list.tsv")%>>%
all_maf = read_tsv("all_pass_with_ascat.maf.gz")
purity = read_tsv("/Volumes/areca42TB2/gdc/purity/Dvir_purity_data.tsv")
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role=`Role in Cancer`)%>>%dplyr::select(gene,role)

patient_list = patient_list %>>%
  left_join(purity,by=c("sample_id","cancer_type","CPE"))%>>%
  mutate(ALL_MAX = pmax(ASCAT,ABSOLUTE,HE_staining,ESTIMATE,LUMP,na.rm = T))


all_maf%>>%
#  filter(ascat_minor==0)%>>%
  tidyr::gather(purity_class,purity,ASCAT,HE_staining,CPE,MAX)%>>%
  mutate(vaf=t_alt/(t_depth*purity))%>>%mutate(vaf=ifelse(vaf>1,1,vaf))%>>%
  ggplot()+
  geom_histogram(aes(x=vaf),binwidth = 0.01)+
  facet_wrap(~ purity_class)

#if using ALL_MAX
all_maf %>>%
  mutate(sample_id_short = str_extract(sample_id,"TCGA-..-....-..."))%>>%
  left_join(purity,by=c("sample_id_short"="sample_id","cancer_type","CPE"))%>>%
  mutate(ALL_MAX = pmax(ASCAT,ABSOLUTE,HE_staining,ESTIMATE,LUMP,na.rm = T)) %>>%
  mutate(vaf=t_alt/(t_depth*ALL_MAX))%>>%
  filter(vaf>0.75)%>>%
  count(sample_id,ALL_MAX)%>>%
  right_join(patient_list%>>%mutate(sample_id=tumor_sample_id))%>>%
  mutate(n=ifelse(is.na(n),0,n))%>>%
  filter(ALL_MAX>0.5)%>>% #change this value
  ggplot(aes(x=ALL_MAX,y=n/mutation_num))+
  geom_point()+
  stat_smooth(method = "lm",colour="red")


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
