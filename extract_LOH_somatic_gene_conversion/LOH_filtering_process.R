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
#purity_class = "MAX"
purity_class = "HE_staining"
patient_list = read_tsv("patient_list.tsv")
dcv_table=read_tsv(paste0("patients_mutect/dcv_status_",purity_class,".tsv"))
all_maf = read_tsv(paste0("all_pass_with_dist_position_",purity_class,".maf.gz"))
all_maf%>>%count(patient_id)
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)

###################################### allele count == 2 ######################################
ac2_maf = all_maf %>>% filter(allele_num<=2)%>>%
  inner_join(dcv_table%>>%filter(around25>0.5)%>>%dplyr::select(sample_id,median_dcv))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNP","indel"))%>>%
  mutate(purity = get(purity_class))%>>%filter(purity>0.5)
ac2_maf %>>%count(patient_id,sample_id)%>>%(?.)%>>%
  count(patient_id)
ac2_maf%>>%count(variant_type)
ac2_maf%>>%count(ascat_minor,variant_type)

## LOH mutation
ac2_binom_test=function(.tbl){
  pbinom( q = .tbl$t_alt*.tbl$purity,size = .tbl$t_depth,
          prob=.tbl$ascat_major/.tbl$allele_num, lower.tail = F)
}
LOH_mutation_ac2 = ac2_maf%>>%mutate(vaf=t_alt/t_depth/purity)%>>%
  mutate(allele_num=ifelse(allele_num==0,1,allele_num))%>>%
  filter(vaf>(allele_num-0.5)/allele_num)%>>%
  nest(t_depth,t_alt,purity,allele_num,ascat_major)%>>%
  mutate(p_value=purrr::map(data,~ac2_binom_test(.)))%>>%
  unnest()%>>%filter(ascat_minor==0 | p_value<0.01)
LOH_mutation_ac2%>>%
  count(ascat_minor,variant_type)
## check ASCAT copy number by depth correlation value ###
conversion=LOH_mutation_ac2 %>>%
  mutate(mutect_dcv_posi_ratio=mutect_dcv_posi/mutect_mut_num)%>>%
  mutate(dcv_check = ifelse(mutect_dcv_posi_ratio<0.1,"low",
                            ifelse(mutect_dcv_posi_ratio>0.9,"dup","correct")))%>>%
  filter(t_depth>20)%>>%
  filter(ascat_minor==1,dcv_check=="correct")
  #count(ascat_minor,dcv_check,variant_type)
###################################### all ######################################
focal_maf = all_maf %>>%
  inner_join(dcv_table%>>%filter(around25>0.5)%>>%dplyr::select(sample_id,median_dcv))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNP","indel"),
         allele_count = ifelse(allele_num<=2,"<=2",">2"),
         genotype = ifelse(ascat_minor==0,"homo","hetero"))%>>%
  mutate(purity = get(purity_class))%>>%filter(purity>0.5)
focal_maf %>>%count(patient_id,sample_id)%>>%(?.)%>>%
  count(patient_id)
focal_maf%>>%count(allele_count,variant_type)
focal_maf %>>% count(allele_count,genotype,variant_type)
## LOH mutation
LOH_mutation = focal_maf%>>%mutate(vaf=t_alt/t_depth/purity)%>>%
  mutate(allele_num=ifelse(allele_num==0,1,allele_num))%>>%
  filter(vaf>(allele_num-0.5)/allele_num)%>>%
  nest(t_depth,t_alt,purity,allele_num,ascat_major)%>>%
  mutate(p_value=purrr::map(data,~ac2_binom_test(.)))%>>%
  unnest()%>>%filter(ascat_minor==0 | p_value<0.01)
LOH_mutation %>>% count(allele_count,genotype,variant_type)
## check ASCAT copy number by depth correlation value ###
LOH_mutation %>>%
  mutate(mutect_dcv_posi_ratio=mutect_dcv_posi/mutect_mut_num)%>>%
  mutate(dcv_check = ifelse(mutect_dcv_posi_ratio<0.1,"low",
                            ifelse(mutect_dcv_posi_ratio>0.8,"dup","correct")))%>>%
  count(allele_count,genotype,dcv_check,variant_type)%>>%
  write_df(paste0("~/Dropbox/work/somatic_gene_conversion/focal_mutations_",purity_class,".tsv"))

#################################### focus only driver gene ###################################
driver_maf = focal_maf %>>% left_join(driver_gene) %>>%
  mutate(cdg=ifelse(is.na(role),"no","cancer_driver_gene"),
         variant_type=ifelse(impact=="HIGH","truncating",
                             ifelse(impact=="MODERATE","missense","other"))) %>>%
  filter(cdg != "no",variant_type!="other")
driver_maf %>>%count(patient_id,sample_id)%>>%(?.)%>>%count(patient_id)
driver_maf%>>%count(variant_type)
driver_maf%>>%count(allele_count,variant_type)
driver_maf %>>% count(allele_count,genotype,variant_type)
## LOH mutation
driver_LOH = driver_maf%>>%mutate(vaf=t_alt/t_depth/purity)%>>%
  mutate(allele_num=ifelse(allele_num==0,1,allele_num))%>>%
  filter(vaf>(allele_num-0.5)/allele_num)%>>%
  nest(t_depth,t_alt,purity,allele_num,ascat_major)%>>%
  mutate(p_value=purrr::map(data,~ac2_binom_test(.)))%>>%
  unnest()%>>%filter(ascat_minor==0 | p_value<0.01)
driver_LOH %>>% count(allele_count,genotype,variant_type)
## check ASCAT copy number by depth correlation value ###
driver_LOH %>>%
  mutate(mutect_dcv_posi_ratio=mutect_dcv_posi/mutect_mut_num)%>>%
  mutate(dcv_check = ifelse(mutect_dcv_posi_ratio<0.1,"low",
                            ifelse(mutect_dcv_posi_ratio>0.8,"dup","correct")))%>>%
  count(allele_count,genotype,dcv_check,variant_type)%>>%
  write_df(paste0("~/Dropbox/work/somatic_gene_conversion/driver_mutations_",purity_class,".tsv"))


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



