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
ac2_maf = all_maf %>>% filter(allele_num<=2,allele_num>0)%>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))
all_maf %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>%
  count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10)%>>%print_tbl.df()


#### excluding variants on short indel or wrong copy number estimation #####
ac2_maf %>>%group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()


tVAFscore=0.75
### extruct truncal LOH mutation (tVAF > )
ac2_maf %>>%filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  ggplot()+geom_histogram(aes(x=tVAF),binwidth = 0.01)+
  facet_wrap(~genotype,scales = "free")+scale_x_continuous(limits = c(0,2))

ac2_maf %>>%filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()


#####################################################################################################
#################################### focus the gene mutation ########################################
#####################################################################################################
all_maf %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  filter(impact=="MODERATE"|impact=="HIGH")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP","nonsynonymous","inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  (?.%>>%count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()






ac2_maf_driver = all_maf %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  inner_join(driver_gene) %>>%
  filter(impact=="MODERATE"|impact=="HIGH")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP","nonsynonymous","inframe_indel")))%>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","inframe_indel")))%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))

ac2_maf_driver %>>%  filter(allele_num<=2,allele_num>0)%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  ggplot()+geom_histogram(aes(x=tVAF),binwidth = 0.02)+
  facet_wrap(~genotype,scales = "free")+scale_x_continuous(limits = c(0,2))

ac2_maf_driver %>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  (?.%>>%count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()




###### somatic gene conversion 確認用 ##############
all_maf%>>%filter(ascat_minor==1,ascat_major==1)%>>%
  inner_join(sample_list%>>%filter(purity>0.5,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75) %>>%
  dplyr::select(sample_id,cancer_type,chr,start,ref,alt,gene,variant_type,consequence,impact,
                purity,screening,mutect_dcv_posi,mutect_mut_num,tVAF)%>>%
  write_df("candidate_conversion_variants.tsv")
  


##############################################################################
############################### 解析 #########################################
##############################################################################
driver_conversion = driver_LOH %>>%
  filter(ascat_minor>0)
driver_conversion%>>%
  dplyr::select(patient_id,sample_id,cancer_type,chr,start,ref,alt,gene,role,variant_classification,variant_type,
                t_depth,t_alt,vaf,HE_staining,MAX,n_depth,n_alt,ascat_major,ascat_minor)%>>%
  write_df(paste0("~/Dropbox/work/somatic_gene_conversion/driver_conversion_",purity_class,".tsv"))
LOH_mutation %>>%filter(ascat_minor>0)%>>%
  count(cancer_type,sample_id)%>>%left_join(all_maf%>>%count(sample_id)%>>%dplyr::rename(alln=n))%>>%
  filter(alln!=1)%>>%
  ggplot()+geom_histogram(aes(x=n/alln),bins=100)+facet_wrap( ~ cancer_type)



