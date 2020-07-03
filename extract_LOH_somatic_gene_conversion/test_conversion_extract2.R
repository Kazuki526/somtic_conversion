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
purity_class = "CPE"
#patient_list = read_tsv("patient_list.tsv")
sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_",purity_class,".maf.gz"))
all_maf%>>%count(patient_id)
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)

###################################### allele count == 2 ######################################
ac2_maf = all_maf %>>% filter(allele_num<=2)%>>%
  inner_join(dcv_table%>>%filter(around25>0.5)%>>%dplyr::select(sample_id,median_dcv))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNP","indel"))%>>%
  mutate(purity = get(purity_class))%>>%filter(purity>0.5)


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
## check ASCAT copy number by depth correlation value ###
conversion=LOH_mutation_ac2 %>>%
  mutate(mutect_dcv_posi_ratio=mutect_dcv_posi/mutect_mut_num)%>>%
  mutate(dcv_check = ifelse(mutect_dcv_posi_ratio<0.1,"low",
                            ifelse(mutect_dcv_posi_ratio>0.9,"dup","correct")))%>>%
  filter(t_depth>20)%>>%
  filter(ascat_minor==1,dcv_check=="correct")

ac2_correct= ac2_maf%>>%filter(ascat_minor==1,t_depth>20)%>>%
  mutate(vaf=t_alt/t_depth/purity)%>>%
  mutate(mutect_dcv_posi_ratio=mutect_dcv_posi/mutect_mut_num)%>>%
  filter(mutect_dcv_posi_ratio>0.1,mutect_dcv_posi_ratio<0.9)


ac2_correct %>>%
  left_join(conversion%>>%dplyr::select(sample_id,chr,start,ref,alt)%>>%
              mutate(conversion="yes"))%>>%
  mutate(conversion=ifelse(is.na(conversion),"non_conversion","conversion"))%>>%
  mutate(window_start=trunc(start/1000000))%>>%
  count(chr,window_start,conversion)%>>%tidyr::spread(conversion,n)%>>%
  mutate(conversion_rate=ifelse(is.na(conversion),0,conversion/non_conversion))%>>%
  ggplot(aes(x=window_start,conversion))+
  geom_bar(stat = "identity")+
  theme_bw()+
  facet_wrap( ~ chr,scales = "free_x")

########################################################################################  
ascat = read_tsv("/Volumes/areca42TB/tcga/CNA/all_patient/all_patient_ascat.tsv.gz")
chrarm = read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")%>>%dplyr::rename(armstart=start,armend=end)

ascat %>>%mutate(chr=paste0("chr",chr))%>>%
  left_join(chrarm,by="chr")%>>%
  filter((startpos>armstart & startpos<armend)|(endpos>armstart & endpos<armend))%>>%
  mutate(startpos=ifelse(startpos<armstart,armstart,startpos),
         endpos  =ifelse(endpos  >armend  ,armend  ,endpos))%>>%
  group_by(sample,chr,arm)%>>%
  filter(n()<10)%>>%
  mutate(region_length=endpos-startpos+1)%>>%
  filter(region_length>100000)

ac2_maf%>>%filter(ascat_minor==0,ascat_major==1,t_depth>20)%>>%
  filter(mutect_dcv_posi_ratio>0.1,mutect_dcv_posi_ratio<0.9)%>>%
  mutate(vaf=t_alt/t_depth/())

tibble(a=c(0,0,0, 1,1,1, 2,2,2, 3,3,3, 4,4,4),b=c(-1,0,1, 0,1,2, 1,2,3, 2,3,4, 3,4,5))%>>%
  ggplot(aes(x=as.factor(a),y=b))+
  geom_violin(scale="count")+
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  geom_abline(aes(intercept=0,slope=1))


