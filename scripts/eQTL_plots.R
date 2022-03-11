#!/usr/bin/Rscript
##
##  eQTL_plots.R
##
##  EDF 2/4/2021
##

setwd("~/projects/LIPA_GWAS/")

library(dplyr)
library(seqminer)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(scales)
library(qqman)



## Read in GWAS table
cad_gwas = read.table("GWAS/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.LIPA_region.txt",
                      header=TRUE, sep = '\t')



## Read in LD table (calculated by plink)
ld_gwas_table = read.table("GTs/plink.ld",
                           header=TRUE) %>%
  mutate(R2_cat = cut(as.numeric(R2), breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),
                              labels=c(0.1,0.3,0.5,0.7,0.9,1.0),
                              include.lowest = FALSE,
                              right=FALSE)) %>%
  select(-c(1,2,3,4))


## Annotate GWAS table w/ LD to lead SNP
cad_gwas_annot = merge(cad_gwas, ld_gwas_table,
                       by.x=c('panel_variant_id','position'),
                       by.y=c('SNP_B','BP_B'),
                       all.x=TRUE) %>%
  mutate(r_2_cat = factor( ifelse(is.na(R2_cat), 0.1, as.character(R2_cat)),
                           levels=c(0.9,0.7,0.5,0.3,0.1,1.0)))


## list of eQTL tissues
eQTL_tiss = read.table("eQTL/LIPA.eQTLs.sig_egenes", header=TRUE,sep='\t') %>% pull(tissue) %>% as.character()
#eQTL_files = paste0("eQTL/results/",eQTL_tiss,"allpairs.txt.gz")

## Read in eQTL files by tiss
eQTL_results = lapply(eQTL_tiss, function(tissi) {
  print(tissi)
  eQTL_file = paste0("eQTL/results/",tissi,".allpairs.txt.gz")
  tabix.read.table(eQTL_file,
                   "chr10:88000000-91000000") %>%
    filter(gene_id == "ENSG00000107798.17") %>%
    mutate(tiss=tissi)
})


## Combine eQTL tissues into single table and label snps/ld levels
eQTL_full_results = do.call('rbind', eQTL_results) %>%
  merge(ld_gwas_table,
        by.x=c('variant_id','pos'), 
        by.y=c('SNP_B','BP_B'),
        all.x=TRUE) %>%
  mutate(snp_name = ifelse(pos==89243047, 'rs1412445',
                         ifelse(pos==89243088, 'rs1320496',
                                ifelse(pos==89243170, 'rs1412444', '')))) %>%
  mutate(r_2_cat = factor( ifelse(is.na(R2_cat), 0.1, as.character(R2_cat)),
                           levels=c(0.9,0.7,0.5,0.3,0.1,1.0)))


# eQTL_results[[1]] %>%
#   ggplot(aes(pos,-log10(pval_nominal))) +
#   geom_point(aes(fill=r_2_cat,
#                  shape=ifelse(r_2_cat==1,'lead','linked'))) +
#   theme_classic() +
#   scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
#   scale_shape_manual(values=c(23,21)) +
#   geom_text_repel(aes(label=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'), 
#                                    as.character(variant_id), '')),
#                   nudge_x=.1)



## For a given tissue, plot the eQTL pvalue vs. the SNP position, with some annotations
tissi='Lung'
eQTL_full_results %>%
  filter(tiss==tissi) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(
    ## SNP color by LD
    fill=r_2_cat,
    ## SNP shape by lead/not
    shape=ifelse(r_2_cat==1,'lead','linked'))) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  ## SNP name label
  geom_text_repel(aes(label=ifelse(snp_name%in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(snp_name), '')),
                  nudge_x=.1) +
  ggtitle(paste(tissi,'eQTL'))
eQTL_full_results %>%
  filter(tiss==tissi) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(r_2_cat==1,'lead','linked'))) +
  theme_classic() +
  xlim(89200000,89300000) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(snp_name), '')),
                  nudge_x=.1) +
  ggtitle(paste(tissi,'eQTL'))


#####################################
########## Don't worry about the rest

eQTL_full_results %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(r_2_cat==1,'lead','linked'))) +
  theme_classic() +
  xlim(89200000,89300000) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(snp_name), '')),
                  nudge_x=.1) +
  facet_wrap(~tiss,scales='free_y')

eQTL_full_results %>%
  arrange(r_2_cat) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  facet_wrap(~tiss,ncol=3,scales='free_y') +
  theme(legend.position='none') +
  scale_x_continuous(labels=c('89.2','89.225','89.25','89.275','89.3'),
                     limits=c(89200000,89300000)) +
  xlab("Position on chromosome 10 (Mb)")


eQTL_enloc = do.call('rbind', lapply(eQTL_tiss, function(tissi) {
  enloc_i = read.table(paste0("ENLOC/enloc_results/",tissi,"_w_Coronary_Artery_Disease_enloc_output.txt.gz"),
             header=TRUE,sep='\t') %>%
    filter(gene_id == "ENSG00000107798.17") %>%
    mutate(tiss=tissi) } ) )

max_logp = eQTL_full_results %>%
  group_by(tiss) %>%
  summarize(max_logp=max(-log10(pval_nominal)))
enloc_labels <- eQTL_enloc %>%
  group_by(tiss) %>%
  summarize(rcp_max=max(rcp), label=paste('rcp =',rcp_max)) %>%
  merge(max_logp)

eQTL_full_results %>%
  arrange(r_2_cat) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  facet_wrap(~tiss,ncol=3,scales='free_y') +
  theme(legend.position='none') +
  scale_x_continuous(labels=c('89.2','89.225','89.25','89.275','89.3'),
                     limits=c(89200000,89300000)) +
  xlab("Position on chromosome 10 (Mb)") +
  ylab("-log10(eQTL p-value)") +
  geom_text(x = 89280000, aes(y=max_logp*.8, label = label), data = enloc_labels)









#other_gene_files = list.files("eQTL/","*region.txt",full.names = TRUE)

gene_names = read.table("expr/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                        header=TRUE,sep='\t') %>%
  select(c(1,2))
  
eQTL_other_gene_results = do.call('rbind', lapply(c('Spleen','Whole_Blood'), function(tissi) {
  print(tissi)
  eQTL_file = paste0("eQTL/",tissi,".allpairs.region.txt")
  read.table(eQTL_file,
             header=FALSE,sep='\t') %>%
    mutate(tiss=tissi)
}) )
names(eQTL_other_gene_results) <- c('gene','var','dtss','ma_samp','ma_count','maf','pval_nom','slope','slope_se','chr','pos','tiss')

eQTL_other_gene = eQTL_other_gene_results %>%
  merge(ld_gwas_table,
        by.x=c('var','chr','pos'), 
        by.y=c('variant','chr','pos')) %>%
  mutate(snp_name = ifelse(pos==89243047, 'rs1412445',
                           ifelse(pos==89243088, 'rs1320496',
                                  ifelse(pos==89243170, 'rs1412444', '')))) %>%
  merge(gene_names)

eQTL_enloc_all = do.call('rbind', lapply(c('Spleen','Whole_Blood'), function(tissi) {
  enloc_i = read.table(paste0("ENLOC/enloc_results/",tissi,"_w_Coronary_Artery_Disease_enloc_output.txt.gz"),
                       header=TRUE,sep='\t') %>%
    mutate(tiss=tissi) } ) )

max_logp = eQTL_other_gene %>%
  filter(pval_nom > 10^-100) %>%
  group_by(tiss) %>%
  summarize(max_logp=max(-log10(pval_nom)))
enloc_labels_all <- eQTL_enloc_all %>%
  rename(gene=gene_id) %>%
  group_by(tiss,gene) %>%
  summarize(rcp_max=max(rcp), label=paste('rcp =',rcp_max)) %>%
  merge(max_logp) %>%
  filter(gene %in% unique(as.character(eQTL_other_gene$gene))) %>%
  merge(gene_names)


eQTL_other_Spleen_toplot = eQTL_other_gene %>%
  filter(tiss=='Spleen') %>%
  group_by(gene) %>%
  filter(min(pval_nom) < 0.00001) %>%
  ungroup() %>%
  mutate(gene=factor(gene, levels=unique(as.character(gene)))) %>%
  arrange(r_2)
enloc_other_Spleen_toplot = enloc_labels_all %>% 
  filter(tiss=='Spleen') %>%
  filter(gene %in% unique(as.character(eQTL_other_Spleen_toplot$gene)))
eQTL_other_Spleen_toplot %>%
  ggplot(aes(pos,-log10(pval_nom))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text(x = 88100000, aes(y=max_logp*.8, label = label), 
            data = enloc_other_Spleen_toplot,
            hjust=0) +
  facet_wrap(~paste(description,gene,sep='\n')) +
  theme(legend.position='none') +
  scale_x_continuous(labels=c('88','88.5','89','89.5','90'),
                     limits=c(88000000,90000000)) +
  xlab("Position on chromosome 10 (Mb)") +
  ylab("-log10(eQTL p-value)") +
  ggtitle("GTEx eQTLs in Spleen")
  

eQTL_other_Blood_toplot = eQTL_other_gene %>%
  filter(tiss=='Whole_Blood',
         pval_nom > 10^-100) %>%
  group_by(gene) %>%
  filter(min(pval_nom) < 0.00001) %>%
  ungroup() %>%
  mutate(gene=factor(gene, levels=unique(as.character(gene)))) %>%
  arrange(r_2)
enloc_other_Blood_toplot = enloc_labels_all %>% 
  filter(tiss=='Whole_Blood') %>%
  filter(gene %in% unique(as.character(eQTL_other_Blood_toplot$gene)))
eQTL_other_Blood_toplot %>%
  ggplot(aes(pos,-log10(pval_nom))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(snp_name %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text(x = 88100000, aes(y=max_logp*.8, label = label), 
            data = enloc_other_Blood_toplot,
            hjust=0) +
  facet_wrap(~paste(description,gene,sep='\n')) +
  theme(legend.position='none') +
  scale_x_continuous(labels=c('88','88.5','89','89.5','90'),
                     limits=c(88000000,90000000)) +
  xlab("Position on chromosome 10 (Mb)") +
  ylab("-log10(eQTL p-value)") +
  ggtitle("GTEx eQTLs in Whole Blood")




##############################################
blueprint_eqtls = read.table("eQTL/blueprint_wp10_qtl_query_result_2021-03-13T00_53_02+01_00.tsv",
                             header=TRUE,sep='\t', comment.char = '')

blueprint_eqtls %>%
  add_row(snp_id='empty',p.bonferroni=1) %>%
  mutate(`SNP` = factor(snp_id, levels=c(as.character(blueprint_eqtls%>%arrange(pos)%>%pull(snp_id)%>%unique()),'empty')),
        `Cell Type`=cell_type) %>%
  ggplot(aes(SNP, -log10(p.bonferroni))) +
  geom_hline(yintercept=-log10(0.05)) +
  geom_point(aes(col=`Cell Type`)) +
  geom_label_repel(aes(label=X..an_group),
                   nudge_x=.3,nudge_y=1,
                   size=3) +
  theme_classic() +
  ylab('-log10(Bonferroni-corrected eQTL p-value)') +
  ggtitle("BluePrint LIPA eQTLs")


# cad_gwas_full = read.table("GWAS/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz",
#                            header=TRUE,sep='\t')
# cad_gwas_full %>%
#   mutate(chr = as.numeric(gsub('chr','',as.character(chromosome)))) %>%
#   manhattan(chr="chr", bp='position', p='pvalue',snp='variant_id')




