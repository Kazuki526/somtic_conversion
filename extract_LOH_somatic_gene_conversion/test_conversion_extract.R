library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/")
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


candidate = read_tsv("somatic_conversion/ascat_candidate.tsv")
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role=`Role in Cancer`)
purity = read_tsv("/Volumes/areca42TB2/gdc/purity/Dvir_purity_data.tsv")
mutect_maf =read_tsv("candidate_patient_all_CT.maf.gz")
patient_info = mutect_maf%>>%
  count(patient_id,purity,comment)

#p<0.01
tbl001 = candidate %>>%
  filter(p_value<0.01) %>>%
  tidyr::separate(tdepth,into=c("t_depth","t_ref","t_alt"),sep=":",convert=T)%>>%
  tidyr::separate(ndepth,into=c("n_depth","n_ref","n_alt"),sep=":",convert=T)%>>%
  tidyr::separate(ascat_allele,into=c("n_major","n_minor"),sep=":",convert=T)%>>%
  mutate(allele_num = n_major+n_minor)%>>%
  left_join(driver_gene%>>%dplyr::select(gene,role))

# narrow down mutedct_maf
mutect_maf_narrow = mutect_maf %>>%
  filter(allele_num!=0)%>>%
  mutate(vaf_correct = t_depth/n_depth/(2*(1-purity)+purity*allele_num))%>>%
  dplyr::select(patient_id,vaf_correct)

.tbl = candidate_comp$data[[1]]
vaf_dist_posi = function(.tbl){
  all_count = filter(mutect_maf_narrow,patient_id==.tbl$pid[1])%>>%count()%>>%{.$n}
  larger_count = filter(mutect_maf_narrow,patient_id==.tbl$pid[1],vaf_correct>.tbl$candidate_vaf_correct[1])%>>%count()%>>%{.$n}
  tibble(mutect_mutation_num=all_count,mutect_vafcor_larger_count=larger_count)
}

candidate_comp = tbl001 %>>%
  mutate(candidate_vaf_correct = t_depth/n_depth/(2*(1-purity)+purity*allele_num),pid=patient) %>>%
  nest(-patient,-chr,-start,-ref,-alt)%>>%#head(1)%>>%
  mutate(counts=purrr::map(data,~vaf_dist_posi(.)))%>>%
  unnest()%>>%dplyr::select(-pid)

candidate_comp %>>% 
  mutate(vaf_correct_dist_posi = mutect_vafcor_larger_count/mutect_mutation_num)%>>%
  write_df("somatic_conversion/ascat_candidate_vaf_correct.tsv")

test=candidate_comp %>>% 
  mutate(vaf_correct_dist_posi = mutect_vafcor_larger_count/mutect_mutation_num)%>>%
  filter(vaf_correct_dist_posi>0.1,vaf_correct_dist_posi<0.9)%>>%
  filter(!is.na(role))%>>%count(IMPACT,role)
  View


