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
chr_length = read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")%>>%
  group_by(chr)%>>%summarise(start=min(start),end=max(end))
purity_cutoff=0.7
tvaf_cutoff=0.8
###################################### allele count == 2 ######################################
ac2_all_maf = all_maf %>>% filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  #filter(tVAF>0.8)%>>%
  mutate(gene_conversion=ifelse(p>=0.95 & tVAF>tvaf_cutoff,1,0))


## 1Mb window density
freq=ac2_all_maf %>>%
  mutate(window=start %/% 10^6)%>>%
  group_by(chr,window)%>>%summarise(gene_conversion=sum(gene_conversion),all=n())%>>%
  ungroup()%>>% filter(all>50)%>>%
  mutate(gene_conversion_rate=gene_conversion/all,
         chr=factor(chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11",
                          "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")))%>>%
  ggplot()+geom_bar(aes(x=window,y=gene_conversion_rate),stat="identity")+
  ylab("Gene conversion rate")+xlab("Position (Mb)")+
  facet_wrap(~chr,scale="free_x")+theme_bw()
ggsave("~/Dropbox/work/somatic_gene_conversion/gene_conversion_density.pdf",freq,height = 8,width=16)


### recombination rate
recombination_rate=read_tsv("~/Dropbox/work/grch38datas/recomb-hg38/genetic_map_GRCh38.tsv")
recomb_conv=recombination_rate%>>%
  mutate(window=pos %/% (10^6))%>>%
  group_by(chr,window)%>>%filter(pos==max(pos))%>>%ungroup()%>>%
  group_by(chr)%>>%
  mutate(recomb_rate=(pos_cm-ifelse(is.na(lag(pos_cm)),0,lag(pos_cm)))/(pos-ifelse(is.na(lag(pos)),0,lag(pos))))%>>%
  left_join(ac2_all_maf %>>%
              mutate(window=start %/% (10^6))%>>%
              group_by(chr,window)%>>%summarise(gene_conversion=sum(gene_conversion),all=n())%>>%
              ungroup()%>>% filter(all>50))%>>%
  filter(!is.na(all))%>>%mutate(gene_conversion_rate=gene_conversion/all)%>>%
  filter(gene_conversion_rate>0)

cor.test(x=recomb_conv$gene_conversion_rate,y=recomb_conv$recomb_rate) # p-value = 5.452e-06
r=cor(x=recomb_conv$gene_conversion_rate,y=recomb_conv$recomb_rate)
#cor_plot=
  recomb_conv%>>%
  ggplot()+geom_point(aes(x=gene_conversion_rate,y=recomb_rate))+
  geom_text(data=tibble(x=0,y=7e-6,label=paste0("r = ",signif(r,3))),
            aes(x=x,y=y,label=label),hjust=0,size=5)+
  theme_bw()+ylab("Recombination rate")+xlab("Gene conversion rate")

#cowplot::plot_grid(freq,cor_plot,nrow = 1,rel_widths = c(2,1),labels = c("A","B"))
ggsave("~/Dropbox/work/somatic_gene_conversion/recomb_conv_corr.pdf",cor_plot,height = 4,width=6)
