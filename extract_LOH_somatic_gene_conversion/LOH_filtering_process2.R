library(tidyverse)
library(pipeR)
library(gridExtra)
library(ggsignif)
loadNamespace('cowplot')
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
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)#%>>%filter(str_detect(role,"TSG"))
#dcv_tbl = read_tsv("patients_mutect/dcv_table_CPE.tsv")

###################################### allele count == 2 ######################################
ac2_maf = all_maf %>>% filter(allele_num<=2,allele_num>0)%>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))
all_maf %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  (?.%>>%count(patient_id)%>>%ungroup()%>>%summarise(mean(n),sd(n))%>>%as.data.frame())%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>%
  count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10)%>>%print_tbl.df()



tVAFscore=0.75
### extruct truncal LOH mutation (tVAF > )
ac2_maf %>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  ggplot()+geom_histogram(aes(x=tVAF),binwidth = 0.01)+
  facet_wrap(~genotype,scales = "free")+scale_x_continuous(limits = c(0,2))

ac2_maf %>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()


#####################################################################################################
#################################### focus the gene mutation ########################################
#####################################################################################################
all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP",
                      ifelse(variant_classification=="Silent","synonymous","nonsynonymous"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","synonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%(?.%>>%count(patient_id)%>>%count())%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%View
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)






ac2_maf_TSG = all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  inner_join(driver_gene%>>%filter(str_detect(role,"TSG"))) %>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP",
                      ifelse(variant_classification=="Silent","synonymous","nonsynonymous"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","synonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count(patient_id)%>>%count())%>>%(?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))

ac2_maf_TSG %>>%  filter(allele_num<=2,allele_num>0)%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  ggplot()+geom_histogram(aes(x=tVAF),binwidth = 0.02)+
  facet_wrap(~genotype,scales = "free")+scale_x_continuous(limits = c(0,2))

ac2_maf_TSG %>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
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
# gene conversin rate of all gene VS TSG
TSG_all_GCrate=all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%filter(genotype=="AB")%>>%
  left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  mutate(role=ifelse(is.na(role),"non TSG",role))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP",
                            ifelse(variant_classification=="Silent","synonymous","nonsynonymous"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","synonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  group_by(genotype,role,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, role, variant_type,bef,)%>>%mutate(ratio=n/bef)%>>%(?.)%>>%
  ggplot()+geom_bar(aes(x=role,y=ratio),stat="identity")+facet_wrap(.~ variant_type)+
  geom_signif(stat="identity",data=data.frame(x=1,xend=2,y=0.045,yend=0.045,p="***",variant_type="truncating"),
              aes(x=x,y=y,xend=xend,yend=yend,annotation=p))+
  theme_bw()+ylab("Gene conversion rate")+theme(axis.title.x = element_blank())
TSG_all_GCrate  
  
TSG_all_LOHfreq=all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%
  left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  mutate(role=ifelse(is.na(role),"non TSG",role))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP",
                                                                ifelse(variant_classification=="Silent","synonymous","nonsynonymous"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","synonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  group_by(genotype,role,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, role, variant_type,bef,)%>>%mutate(ratio=n/bef)%>>%(?.)%>>%
  ggplot()+geom_bar(aes(x=role,y=n,fill=genotype),color="black",stat="identity",position = "fill")+facet_wrap(.~ variant_type)+
  theme_bw()+ylab("Frequency in LOH mutation")+theme(axis.title.x = element_blank())
TSG_all_LOHfreq  

GC_rate = all_maf %>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%filter(genotype=="AB")%>>%
  #left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  #mutate(role=ifelse(!is.na(role),role,ifelse(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent","Gene region","Whole genome")))%>>%
  mutate(role=ifelse(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent","Gene region","Whole genome"))%>>%
  group_by(genotype,role)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, role,bef)%>>%mutate(ratio=n/bef)%>>%(?.)%>>%
  ggplot()+geom_bar(aes(x=role,y=ratio),stat="identity")+
  geom_signif(stat="identity",data=data.frame(x=1,xend=2,y=0.048,yend=0.048,p="***"),
              aes(x=x,y=y,xend=xend,yend=yend,annotation=p))+
  theme_bw()+ylab("Gene conversion rate")+theme(axis.title.x = element_blank())
GC_rate

cowplot::plot_grid(TSG_all_GCrate,
                   cowplot::plot_grid(TSG_all_LOHfreq,GC_rate,nrow=1,rel_widths=c(3,1),labels=c("B","C")),
                   nrow = 2,labels=c("A",""))
cowplot::plot_grid(TSG_all_GCrate,TSG_all_LOHfreq,GC_rate,nrow = 1,rel_widths = c(10,1,10))
ggsave("~/Dropbox/work/somatic_gene_conversion/TSGvsAllgene.pdf",height = 8,width = 10)


























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

ac2_maf_onco = all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>% 
  inner_join(driver_gene%>>%filter(str_detect(role,"oncogene"))) %>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","truncating",ifelse(variant_type=="SNP",
                                                                ifelse(variant_classification=="Silent","synonymous","nonsynonymous"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("truncating","nonsynonymous","synonymous","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count(patient_id)%>>%count())%>>%(?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))


ac2_maf_onco %>>%
  filter(variant_type!="inframe_indel")%>>%
  (?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  (?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75)%>>%
  count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%print_tbl.df()

