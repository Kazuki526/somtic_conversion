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
print_tbl.df = function(x,..){print(as.data.frame(x))}

#purity_class = "MAX"
purity_class = "CPE"
#patient_list = read_tsv("patient_list.tsv")
sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_",purity_class,".maf.gz"))
all_maf%>>%count(patient_id)
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)%>>%
  filter(str_detect(role,"TSG"))

###################################### allele count == 2 ######################################
if(0){# without binom
  purity_cutoff=0.7
  tvaf_cutoff=0.8
all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  (?.%>>%count(patient_id)%>>%count())%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>tvaf_cutoff)%>>%View
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%
  mutate(freq=n/sum(n))
}


#cutoff by binom
LOH_filtering = function(purity_cutoff,tvaf_cutoff){
  sample_list = sample_list%>>%
    mutate(screening=ifelse(purity>purity_cutoff,ifelse(dcv_median95/2 > dcv_sd95,NA,"dcv_out"),"purity_out"))
  all_maf %>>%
    filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
    inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
    mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
    mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
    #(?.%>>%count(patient_id)%>>%count())%>>%
    #(?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
    group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
    filter(allele_num<=2,allele_num>0)%>>% 
    #(?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
    group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
    mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
    mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
    mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
    mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
    filter((allele_num==2 & p>=0.95)|(allele_num==1),tVAF>tvaf_cutoff)%>>%
    count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%
    mutate(freq=n/sum(n))%>>%
    dplyr::select(genotype,variant_type,n,freq)
}

for(p in c(0,0.5,0.6,0.7,0.8,0.9)){
  print(paste0("purity == ",p))
  outbl=LOH_filtering(p,0.75)%>>%dplyr::rename(`tVAF>0.75`= n, `tVAF>0.75 freq`=freq)
  for(tvaf in c(0.8,0.85,0.9,0.95)){
    outbl=left_join(outbl,LOH_filtering(p,tvaf)%>>%dplyr::rename(!!paste0("tVAF>",tvaf):=n,!!paste0("tVAF>",tvaf," freq"):=freq))
  }
  outbl
  write_df(outbl,paste0("~/Dropbox/work/somatic_gene_conversion/cutoff_table/purity",p*100,".tsv"))
}




##### purityと　genotype (A, AA, AB)の関係 #########
all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(screening!="dcv_out"|is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  filter(allele_num<=2,allele_num>0)%>>%
  mutate(purity_class=((purity-0.000001) %/%0.1)/10)%>>%
  mutate(purity_class=ifelse(purity_class<0.5,0,purity_class))%>>%
  mutate(purity_class=ifelse(purity_class==0,"0-0.5",ifelse(purity_class==1,"0.9-1",paste0(purity_class,"-",purity_class+0.1))))%>>%
  count(purity_class,genotype)%>>%
  ggplot()+
  geom_bar(aes(x=purity_class,y=n,fill=genotype),color="black",stat="identity",position = "fill")+
  theme_bw()
all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(screening!="dcv_out"|is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  filter(allele_num<=2,allele_num>0)%>>%
  mutate(purity_class=((purity-0.000001) %/%0.1)/10)%>>%
  mutate(purity_class=ifelse(purity_class<0.5,0,purity_class))%>>%
  mutate(purity_class=ifelse(purity_class==0,"0-0.5",ifelse(purity_class==1,"0.9-1",paste0(purity_class,"-",purity_class+0.1))))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  ggplot()+
  geom_histogram(aes(x=tVAF),binwidth = 0.01)+
  facet_grid(genotype ~ purity_class,scales = "free")+
  theme_bw()





