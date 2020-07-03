library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/patients_mutect/")
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
print_tbl.df = function(x,..){print(as.data.frame(x))}

patient_list = read_tsv("../patient_list.tsv")
sample_list = patient_list%>>%
  group_by(patient_id)%>>%
  filter(mutation_num==max(mutation_num))%>>%
  filter(tumor_sample_id==first(tumor_sample_id))%>>%ungroup

dcv_status=function(file,purity_class="CPE"){
  read_tsv(file,col_types = "ccccicccciiidddd")%>>%
    filter(allele_num!=0, t_depth>=10, n_depth>=10) %>>%
    {if(purity_class=="CPE"){.%>>%mutate(CPE=ifelse(is.na(CPE),ifelse(is.na(HE_staining),ASCAT,HE_staining),CPE))}else{.}}%>>%
    mutate(purity=get(purity_class))%>>%
    filter(!is.na(purity))%>>%
    mutate(dcv =t_depth/n_depth/((1-purity)+allele_num*purity/2))%>>%
    dplyr::select(sample_id,purity,dcv)%>>%
    group_by(sample_id)%>>%
    inner_join(sample_list%>>%dplyr::select(tumor_sample_id)%>>%
                dplyr::rename(sample_id=tumor_sample_id))%>>%
    group_by(sample_id)%>>%
    mutate(n=n(),median=median(dcv),mean=mean(dcv),sd=sd(dcv),rank=row_number(dcv))%>>%
    filter(rank<n*0.95)%>>%
    summarise(n=first(n),median=first(median),mean=first(mean),sd=first(sd),purity=first(purity),
              median95=median(dcv),mean95=mean(dcv),sd95=sd(dcv))
}
if(0){
dcv_status("patient_ACC.maf.gz")
dcv_table = tibble(file=list.files(".",pattern = "^patient_[A-Z]+.maf.gz",full.names = T))%>>%
  mutate(tbl=purrr::map(file,~dcv_status(.)))%>>%
  unnest()%>>%dplyr::select(-file)
write_df(dcv_table,"dcv_table_CPE.tsv")
}
dcv_table=read_tsv("dcv_table_CPE.tsv")
all=dcv_table %>>%filter(purity>0.5)%>>%(?.)%>>%
  mutate(`sd/median`=sd95/median95)%>>%
  ggplot()+geom_histogram(aes(x=`sd/median`),bins=100)+theme_bw()
und2=dcv_table %>>%filter(purity>0.5)%>>%
  mutate(`sd/median`=sd95/median95)%>>%
  filter(`sd/median`<2)%>>%
  ggplot()+geom_histogram(aes(x=`sd/median`),bins=100)+theme_bw()
plot=cowplot::plot_grid(all,und2,nrow = 2)
plot
ggsave("~/Dropbox/innan/somatic_gene_conversion_MS/supply/DCV_SD_mean_distribution.pdf",height = 10,width = 10)

sample_list %>>%dplyr::select(-ASCAT,-HE_staining,-CPE,-MAX)%>>%
  left_join(dcv_table%>>%dplyr::rename(tumor_sample_id=sample_id,mutect_n=n,dcv_median95=median95,dcv_sd95=sd95)%>>%
              dplyr::select(tumor_sample_id, purity, mutect_n, dcv_median95, dcv_sd95))%>>%
  mutate(screening=ifelse(purity>0.5,ifelse(dcv_median95/2 > dcv_sd95,NA,"dcv_out"),"purity_out")) %>>%(?.%>>%count(screening))%>>%
  summarise(mean(mutation_num),sd(mutation_num),sum(mutation_num))%>>%print_tbl.df()
  #filter(is.na(screening))%>>%summarise(mean(mutation_num),sd(mutation_num),sum(mutation_num))%>>%print_tbl.df()
  #summarise(mean(mutect_n),sd(mutect_n),min(mutect_n))%>>%print_tbl.df()
  #write_df("../sample_list.tsv")
####################################################################################################
# DCV plot
pick_dcv = function(file,purity_class="CPE"){
  read_tsv(file,col_types = "ccccicccciiidddd")%>>%
    filter(allele_num!=0, t_depth>=10, n_depth>=10) %>>%
    {if(purity_class=="CPE"){.%>>%mutate(CPE=ifelse(is.na(CPE),ifelse(is.na(HE_staining),ASCAT,HE_staining),CPE))}else{.}}%>>%
    mutate(purity=get(purity_class))%>>%
    filter(!is.na(purity),purity>0.5)%>>%
    mutate(dcv =t_depth/n_depth/((1-purity)+allele_num*purity/2))%>>%
    dplyr::select(sample_id,purity,dcv)%>>%
    group_by(sample_id)%>>%mutate(dcv_correct=dcv/median(dcv))%>>%ungroup()
}
all_dcv = tibble(file=list.files(".",pattern = "^patient_[A-Z]+.maf.gz",full.names = T))%>>%
  mutate(tbl=purrr::map(file,~pick_dcv(.)))%>>%
  unnest()%>>%dplyr::select(-file)

all_dcv_plot=all_dcv%>>%mutate(dcv_correct=ifelse(dcv_correct>5,5,dcv_correct))%>>%
  ggplot()+geom_histogram(aes(x=dcv_correct),bins = 100)+
  theme_bw()+ylab("Count")+xlab("")

sample1="TCGA-CR-6472-01A-11D-1870-08"
sample2="TCGA-CV-7104-01A-11D-2012-08"
sample3="TCGA-A2-A259-01A-11D-A16D-09"
sample4="TCGA-AR-A2LL-01A-11D-A17W-09"
sample5="TCGA-D8-A27F-01A-11D-A16D-09"

sample_plot=function(sampleN){
  all_dcv%>>%filter(sample_id == sampleN)%>>%
    ggplot()+geom_histogram(aes(x=dcv_correct),binwidth = 0.05)+
    theme_bw()+ylab("Count")+xlab("")+
    ggtitle(paste0(sampleN))
}
sample_plot("TCGA-CR-6472-01A-11D-1870-08")
sample_plot("TCGA-CV-7104-01A-11D-2012-08")

by_samples=cowplot::plot_grid(sample_plot(sample1),sample_plot(sample3),sample_plot(sample2),sample_plot(sample4),
                              nrow = 2,labels = c("B","D","C","E"))
by_samples
.plot=cowplot::plot_grid(all_dcv_plot,by_samples,nrow=2,rel_heights = c(1.5,2),labels=c("A",NULL))
ggsave("~/Dropbox/innan/somatic_gene_conversion_MS/supply/DCV_distributions.pdf",.plot,height = 6,width = 10)
dcv_table%>>%filter(median95<2*sd95,purity>0.5)%>>%View()
