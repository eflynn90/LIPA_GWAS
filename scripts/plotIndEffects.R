#!/usr/bin/Rscript
##
##
##  plotIndEffects.R
##
##  EDF 2/19/21
##

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

setwd("~/projects/LIPA_GWAS/")

GTs_raw=read.table("GTs/all_LIPA_region_GTs.vcf",header = TRUE, comment.char = '', stringsAsFactors = FALSE)
eQTLs=read.table("eQTL/LIPA.eQTLs.sig_egenes",header=TRUE)
tiss_col = read.table("../EFF_EXPR/input_data/tissue_translation_colors_v8.txt",header=TRUE,sep='\t')

GTs_snps = GTs_raw %>%
  filter(ID %in% c('chr10_89243088_C_T_b38', 'chr10_89243047_C_T_b38'))
DS_snps = GTs_snps %>%
  pivot_longer(cols=10:847) %>%
  separate(value,c('A1','A2'),sep='[/]') %>%
  mutate(DS=as.numeric(A1)+as.numeric(A2)) %>%
  select(ID,name,DS) %>%
  pivot_wider(names_from=ID, values_from=DS)



expr_files=list.files(path="expr/by_tiss/", 
                      pattern="*ENSG00000107798*", 
                      full.names = TRUE)
# all_lms_nocov = do.call('rbind', lapply(expr_files, function(e_file) {
#   
#   print(e_file)
#   tissi=strsplit(e_file,'[/.]')[[1]][4]
#   
#   expr=read.table(e_file,header=TRUE,comment.char='')
#   expr_long=expr %>%
#     pivot_longer(cols = 5:ncol(expr))
#   
#   tiss_data = merge(DS_snps, expr_long, by.x='name', by.y='name') %>%
#     merge(cov_tiss, by.x='name', by.y='row.names')
#   
#   lm_PU1 = summary(lm( value ~ chr10_89243088_C_T_b38, tiss_data))
#   lm_GWAS = summary(lm( value ~ chr10_89243047_C_T_b38, tiss_data))
#   lm_add = summary(lm( value ~ chr10_89243088_C_T_b38 + chr10_89243047_C_T_b38, tiss_data))
#   lm_int = summary(lm( value ~ chr10_89243088_C_T_b38 * chr10_89243047_C_T_b38, tiss_data))
#   lm_PU1_noGWAS = summary(lm( value - chr10_89243047_C_T_b38 * lm_GWAS$coefficients[2,1] ~ chr10_89243088_C_T_b38, tiss_data))
#   lm_GWAS_noPU1 = summary(lm( value - chr10_89243088_C_T_b38 * lm_PU1$coefficients[2,1] ~ chr10_89243047_C_T_b38, tiss_data))
# 
#   do.call('rbind', lapply( list(lm_PU1$coefficients,lm_GWAS$coefficients,
#                               lm_add$coefficients ,lm_int$coefficients,
#                               lm_PU1_noGWAS$coefficients ,lm_GWAS_noPU1$coefficients), 
#                            function(lmi) {
#                              print(lmi)
#             data.frame(b = lmi[-1,1],
#                        se = lmi[-1,2],
#                        p = lmi[-1,4])
#           }) ) %>%
#     mutate(
#       model = factor(c("lm_PU1",'lm_GWAS','lm_add' ,'lm_add' ,'lm_int','lm_int','lm_int','lm_PU1_noGWAS' ,'lm_GWAS_noPU1'),
#                      levels=c('lm_GWAS','lm_PU1','lm_add','lm_int','lm_PU1_noGWAS','lm_GWAS_noPU1')),
#       SNP = c('PU1','GWAS','PU1','GWAS','PU1','GWAS','int','PU1','GWAS'),
#       tiss=tissi)
#     
# }))
# 
# write.table(all_lms_nocov, "ind_eff/all_lms_nocov.txt",
#             col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

all_lms_cov = do.call('rbind', lapply(expr_files, function(e_file) {
  
  print(e_file)
  tissi=strsplit(e_file,'[/.]')[[1]][4]
  
  cov_tiss = read.table(paste0("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/", tissi, ".v8.covariates.txt"),
                        header=TRUE,sep='\t') %>%
    column_to_rownames(var='ID') %>%
    t()
  
  expr=read.table(e_file,header=TRUE,comment.char='')
  expr_long=expr %>%
    pivot_longer(cols = 5:ncol(expr))
  
  tiss_data = merge(DS_snps, expr_long, by.x='name', by.y='name') %>%
    merge(cov_tiss, by.x='name', by.y='row.names')
  cov_col = names(tiss_data)[ seq(which(names(tiss_data) == 'PC1'), length(names(tiss_data)) )]
  
  lm_PU1 = summary(lm( as.formula(paste("value ~ chr10_89243088_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  lm_GWAS = summary(lm( as.formula(paste("value ~ chr10_89243047_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  lm_add = summary(lm( as.formula(paste("value ~ chr10_89243088_C_T_b38 + chr10_89243047_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  lm_int = summary(lm( as.formula(paste("value ~ chr10_89243088_C_T_b38 * chr10_89243047_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  lm_PU1_noGWAS = summary(lm( as.formula(paste("value - chr10_89243047_C_T_b38 * lm_GWAS$coefficients[2,1] ~ chr10_89243088_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  lm_GWAS_noPU1 = summary(lm( as.formula(paste("value - chr10_89243088_C_T_b38 * lm_PU1$coefficients[2,1] ~ chr10_89243047_C_T_b38 +", paste(cov_col, collapse='+'))), tiss_data))
  
  do.call('rbind', lapply( list(lm_PU1$coefficients,lm_GWAS$coefficients,
                                lm_add$coefficients ,lm_int$coefficients,
                                lm_PU1_noGWAS$coefficients ,lm_GWAS_noPU1$coefficients), 
                           function(lmi) {
                             print(lmi)
                             data.frame(term = row.names(lmi[-1,]),
                                        b = lmi[-1,1],
                                        se = lmi[-1,2],
                                        p = lmi[-1,4])
                           }) ) %>%
    filter(! term %in% cov_col) %>%
    mutate(
      model = factor(c("lm_PU1",'lm_GWAS','lm_add' ,'lm_add' ,'lm_int','lm_int','lm_int','lm_PU1_noGWAS' ,'lm_GWAS_noPU1'),
                     levels=c('lm_GWAS','lm_PU1','lm_add','lm_int','lm_PU1_noGWAS','lm_GWAS_noPU1')),
      SNP = c('PU1','GWAS','PU1','GWAS','PU1','GWAS','int','PU1','GWAS'),
      tiss=tissi)
  
}))

write.table(all_lms_cov, "ind_eff/all_lms_cov.rs1412445.txt",
            col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)


ggplot(all_lms_cov,aes(tiss, b)) +
  geom_hline(yintercept=0) +
  geom_point(aes(fill=SNP, shape=SNP, alpha=ifelse(-log10(p)>3,3,-log10(p)))) +
  geom_point(aes(fill=SNP, shape=SNP, alpha=0)) +
  scale_fill_manual(values=c('darkorchid3','gray','chartreuse3')) +
  scale_shape_manual(values=c(23,21,21)) +
  facet_wrap(~model) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))

all_lms_cov %>%
  filter(tiss=='Whole_Blood') %>%
  ggplot(aes(model, b)) +
  geom_hline(yintercept=0) +
  geom_col(aes(group=SNP, width=c(.05,.05,.1,.1,.15,.15,.15,.05,.05)),
           width=.03,
           position=position_dodge(.2)) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=ifelse(-log10(p)>10,10,-log10(p)),
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=0,
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('darkorchid3','gray','chartreuse3')) +
  scale_shape_manual(values=c(23,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))



all_lms_cov %>%
  filter(tiss=='Whole_Blood') %>%
  filter(model %in% c('lm_GWAS','lm_PU1','lm_int')) %>%
  ggplot(aes(model, b)) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(ymin=0, ymax=b, group=SNP),
           position=position_dodge(.2)) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=ifelse(-log10(p)>10,10,-log10(p)),
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=0,
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('darkorchid3','gray','chartreuse3')) +
  scale_shape_manual(values=c(23,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))



all_lms_cov %>%
  filter(tiss=='Whole_Blood') %>%
  filter(model %in% c('lm_GWAS','lm_PU1','lm_int','lm_PU1_noGWAS','lm_GWAS_noPU1')) %>%
  ggplot(aes(model, b)) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(ymin = b-se, ymax = b+se, group=SNP),
           position=position_dodge(.2)) +
  geom_point(aes(fill=SNP, shape=SNP,
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('darkorchid3','gray','chartreuse3')) +
  scale_shape_manual(values=c(23,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))






############ GWAS results
orig_gwas = read.table("GWAS/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.LIPA_region.txt",
                       header=TRUE,sep='\t')
cojo_gtex = read.table("ind_eff/cojo/test_all_GTEx.cma.cojo",
                       header=TRUE, sep='\t')

lms_to_plot_eqtl = all_lms_cov %>%
  filter(tiss == 'Whole_Blood',
         model %in% c('lm_GWAS','lm_PU1','lm_int','lm_PU1_noGWAS')) %>%
  mutate(dataset='Blood eQTL')

lms_to_plot_orig_gwas = orig_gwas %>%
  filter(variant_id %in% c('rs1320496','rs1412445')) %>%
  mutate(term=panel_variant_id,
         b=effect_size, se=standard_error, p=pvalue,
         model=c('lm_PU1','lm_GWAS'),
         SNP=c('PU1','GWAS'),
         tiss='GWAS',
         dataset='CAD GWAS') %>%
  select(names(lms_to_plot_eqtl))

lms_to_plot_cond_gwas = cojo_gtex %>%
  filter(SNP == 'chr10_89243088_C_T_b38') %>%
  mutate(term=SNP,
         b=bC, se=bC_se, p=pC,
         model=c('lm_PU1_noGWAS'),
         SNP=c('PU1'),
         tiss='GWAS',
         dataset='CAD GWAS') %>%
  select(names(lms_to_plot_eqtl))

lms_to_plot = rbind(lms_to_plot_eqtl, lms_to_plot_orig_gwas) %>%
  rbind(lms_to_plot_cond_gwas)
  




lms_to_plot %>%
  ggplot(aes(model, b)) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(ymin=0, ymax=b, group=SNP),
                 position=position_dodge(.2)) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=ifelse(-log10(p)>10,10,-log10(p)),
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  geom_point(aes(fill=SNP, shape=SNP, 
                 alpha=0,
                 group=SNP),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('darkorchid3','gray','chartreuse3')) +
  scale_shape_manual(values=c(23,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  facet_wrap(~dataset, scales = 'free')



############ CURRENT PLOT TO SHARE##################################
lms_to_plot %>%
  mutate(model_lab = factor(ifelse(model=='lm_PU1', "rs1320496 only",
                              ifelse(model=='lm_GWAS', "rs1412445 only", 
                                     ifelse(model=='lm_int', "interaction",
                                            ifelse(model=='lm_PU1_noGWAS',"rs1320496\nconditional",'')))),
                              levels=c("rs1412445 only",'rs1320496 only','interaction','rs1320496\nconditional'))) %>%
  mutate(`SNP/term`= ifelse(SNP=='GWAS','rs1412445',
                         ifelse(SNP=='PU1','rs1320496',
                                ifelse(SNP=='int','interaction','')))) %>%
  ggplot(aes(model_lab, b)) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(ymin = b-se, ymax = b+se, group=`SNP/term`),
                 position=position_dodge(.2)) +
  geom_point(aes(fill=`SNP/term`, shape=`SNP/term`,
                 group=`SNP/term`),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('gray','chartreuse3','red')) +
  scale_shape_manual(values=c(21,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  facet_wrap(~dataset, scales = 'free') +
  ylab("Effect Size") +
  xlab("Model")

lms_to_plot %>%
  mutate(model_lab = factor(ifelse(model=='lm_PU1', "rs1320496 only",
                                   ifelse(model=='lm_GWAS', "rs1412444 only", 
                                          ifelse(model=='lm_int', "interaction",
                                                 ifelse(model=='lm_PU1_noGWAS',"rs1320496\nconditional",'')))),
                            levels=c("rs1412444 only",'rs1320496 only','interaction','rs1320496\nconditional'))) %>%
  mutate(`SNP/term`= factor(ifelse(SNP=='GWAS','rs1412444',
                            ifelse(SNP=='PU1','rs1320496',
                                   ifelse(SNP=='int','interaction',''))),
                            levels=c('rs1412444','rs1320496','interaction'))) %>%
  ggplot(aes(model_lab, b)) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(ymin = b-se, ymax = b+se, group=`SNP/term`),
                 position=position_dodge(.2)) +
  geom_point(aes(fill=`SNP/term`, shape=`SNP/term`,
                 group=`SNP/term`),
             position = position_dodge(.2),
             size=2) +
  scale_fill_manual(values=c('purple','chartreuse3','gray')) +
  scale_shape_manual(values=c(23,21,21)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  facet_wrap(~dataset, scales = 'free') +
  ylab("Effect Size") +
  xlab("Model")



#### Whole Blood LIPA expression
e_file=expr_files[49]
print(e_file)
tissi=strsplit(e_file,'[/.]')[[1]][4]

expr=read.table(e_file,header=TRUE,comment.char='')
expr_long=expr %>%
  pivot_longer(cols = 5:ncol(expr))

tiss_data = merge(DS_snps, expr_long, by.x='name', by.y='name')
min_expr = min(tiss_data$value)
max_expr = max(tiss_data$value)

tiss_data %>%
  filter(!is.na(chr10_89243088_C_T_b38), !is.na(chr10_89243047_C_T_b38)) %>%
  mutate(snp88 = factor( ifelse(chr10_89243088_C_T_b38==0, 'C/C',
                               ifelse(chr10_89243088_C_T_b38==1, 'C/T',
                                      ifelse(chr10_89243088_C_T_b38==2, 'T/T', ''))),
                         levels=c('C/C','C/T','T/T'))) %>%
  mutate(snp70 = factor( ifelse(chr10_89243047_C_T_b38==0, 'C/C',
                                ifelse(chr10_89243047_C_T_b38==1, 'C/T',
                                       ifelse(chr10_89243047_C_T_b38==2, 'T/T', ''))),
                         levels=c('T/T','C/T','C/C'))) %>%
  group_by(snp88, snp70) %>%
  mutate(med = median(value),
         min = min(value),
         max = max(value),
         q1 = quantile(value, 0.25),
         q3 = quantile(value, 0.75)) %>%
  ggplot() +
  geom_hline(aes(yintercept=med),col='darkorchid3') +
  #geom_hline(aes(yintercept=q1)) +
  #geom_hline(aes(yintercept=q3)) +
  geom_vline(aes(xintercept=med),col='darkorchid3') +
  #geom_vline(aes(xintercept=q1)) +
  #geom_vline(aes(xintercept=q3)) +
  #geom_abline(slope=1,intercept=.5) +
  #geom_abline(slope=1,intercept=-.5) +
  geom_point(aes(x=value, y=value),
             position=position_jitter(width=.2,height=.2)) +
  geom_point(aes(x=med,y=med), col='darkorchid3') +
  #geom_segment(aes(x=q1+.5, y=q1-.5, xend=q3+.5, yend=q3-.5), col='darkorchid3') +
  #geom_segment(aes(x=q1-.5, y=q1+.5, xend=q3-.5, yend=q3+.5), col='darkorchid3') +
  #geom_segment(aes(x=q1+.5, y=q1-.5, xend=q1-.5, yend=q1+.5), col='darkorchid3') +
  #geom_segment(aes(x=q3-.5, y=q3+.5, xend=q3+.5, yend=q3-.5), col='darkorchid3') +
  #geom_segment(aes(x=med-.5, y=med+.5, xend=med+.5, yend=med-.5), col='darkorchid3') +
  #geom_segment(aes(x=min, y = min, xend=q1, yend=q1), col='darkorchid3') +
  #geom_segment(aes(x=max, y = max, xend=q3, yend=q3), col='darkorchid3') +
  geom_freqpoly(aes(x=value,
                    y=(..count..)*(max_expr-min_expr)*4/sum(..count..) + min_expr), 
                col='black', fill=NA, bins=20) +
  geom_freqpoly(aes(y=value,
                    x=(..count..)*(max_expr-min_expr)*4/sum(..count..) * -1 + max_expr), 
                col='black', fill=NA, bins=20) +
  facet_grid(cols=vars(snp88),
             rows=vars(snp70),
             switch='y') +
  scale_y_continuous(position='right') +
  xlab('Whole Blood LIPA Expression\n(DESeq log2)') +
  ylab('Whole Blood LIPA Expression\n(DESeq log2)')





#### Spleen LIPA expression

e_file=expr_files[43]
print(e_file)
tissi=strsplit(e_file,'[/.]')[[1]][4]

expr=read.table(e_file,header=TRUE,comment.char='')
expr_long=expr %>%
  pivot_longer(cols = 5:ncol(expr))

tiss_data = merge(DS_snps, expr_long, by.x='name', by.y='name')
min_expr = min(tiss_data$value)
max_expr = max(tiss_data$value)

tiss_data %>%
  filter(!is.na(chr10_89243088_C_T_b38), !is.na(chr10_89243047_C_T_b38)) %>%
  mutate(snp88 = factor( ifelse(chr10_89243088_C_T_b38==0, 'C/C',
                                ifelse(chr10_89243088_C_T_b38==1, 'C/T',
                                       ifelse(chr10_89243088_C_T_b38==2, 'T/T', ''))),
                         levels=c('C/C','C/T','T/T'))) %>%
  mutate(snp70 = factor( ifelse(chr10_89243047_C_T_b38==0, 'C/C',
                                ifelse(chr10_89243047_C_T_b38==1, 'C/T',
                                       ifelse(chr10_89243047_C_T_b38==2, 'T/T', ''))),
                         levels=c('T/T','C/T','C/C'))) %>%
  group_by(snp88, snp70) %>%
  mutate(med = median(value),
         min = min(value),
         max = max(value),
         q1 = quantile(value, 0.25),
         q3 = quantile(value, 0.75)) %>%
  ggplot() +
  geom_hline(aes(yintercept=med),col='darkorchid3') +
  #geom_hline(aes(yintercept=q1)) +
  #geom_hline(aes(yintercept=q3)) +
  geom_vline(aes(xintercept=med),col='darkorchid3') +
  #geom_vline(aes(xintercept=q1)) +
  #geom_vline(aes(xintercept=q3)) +
  #geom_abline(slope=1,intercept=.5) +
  #geom_abline(slope=1,intercept=-.5) +
  geom_point(aes(x=value, y=value),
             position=position_jitter(width=.2,height=.2)) +
  geom_point(aes(x=med,y=med), col='darkorchid3') +
  #geom_segment(aes(x=q1+.5, y=q1-.5, xend=q3+.5, yend=q3-.5), col='darkorchid3') +
  #geom_segment(aes(x=q1-.5, y=q1+.5, xend=q3-.5, yend=q3+.5), col='darkorchid3') +
  #geom_segment(aes(x=q1+.5, y=q1-.5, xend=q1-.5, yend=q1+.5), col='darkorchid3') +
  #geom_segment(aes(x=q3-.5, y=q3+.5, xend=q3+.5, yend=q3-.5), col='darkorchid3') +
  #geom_segment(aes(x=med-.5, y=med+.5, xend=med+.5, yend=med-.5), col='darkorchid3') +
  #geom_segment(aes(x=min, y = min, xend=q1, yend=q1), col='darkorchid3') +
  #geom_segment(aes(x=max, y = max, xend=q3, yend=q3), col='darkorchid3') +
  geom_freqpoly(aes(x=value,
                    y=(..count..)*(max_expr-min_expr)*4/sum(..count..) + min_expr), 
                col='black', fill=NA, bins=20) +
  geom_freqpoly(aes(y=value,
                    x=(..count..)*(max_expr-min_expr)*4/sum(..count..) * -1 + max_expr), 
                col='black', fill=NA, bins=20) +
  facet_grid(cols=vars(snp88),
             rows=vars(snp70),
             switch='y') +
  scale_y_continuous(position='right') +
  xlab('Spleen LIPA Expression\n(DESeq log2)') +
  ylab('Spleen LIPA Expression\n(DESeq log2)')

