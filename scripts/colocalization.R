#!/usr/bin/Rscript
##
##  colocalization.R
##
##  EDF 2/4/2021
##

library(dplyr)
library(seqminer)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggrepel)


## GWAS file with p-values
cad_gwas = read.table("GWAS/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.LIPA_region.txt",
                      header=TRUE, sep = '\t')


### Code to get LD with lead SNP 
### (takes a bit to run, best to run once and then write to file)
# gts_44 = tabix.read.table("GTs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz",
#                           "chr10:89243170-89243170") %>%
#   select(10:847) %>%
#   t() %>%
#   as.data.frame() %>%
#   separate(V1,"/",into=c('D1','D2')) %>%
#   mutate(DS = as.numeric(D1) + as.numeric(D2))
# ld_gwas = lapply(as.character(cad_gwas$panel_variant_id), function(var_i) {
#   print(var_i)
#   chr='chr10'
#   pos=strsplit(var_i,"_")[[1]][2]
#   tabix_i = tabix.read.table("GTs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz",
#                    paste0(chr,":",pos,"-",pos))
#   if(nrow(tabix_i) > 0) {
#   gts_i = tabix_i %>%
#     filter(ID == as.character(var_i)) %>%
#     select(10:847) %>%
#     t() %>%
#     as.data.frame() %>%
#     separate(V1,"/",into=c('D1','D2')) %>%
#     mutate(DS = as.numeric(D1) + as.numeric(D2))
#   r = cor(gts_44$DS,gts_i$DS,use='na.or.complete')
#   } else {
#     r=NA
#   }
#   data.frame(variant=var_i,chr,pos,r,r_2=r*r)
# })

## Read in saved LD table
ld_gwas_table = read.table("GTs/LD.topsnp.txt",
                           header=TRUE,sep='\t') %>%
  mutate(r_2_cat = factor(cut(as.numeric(r_2), breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),
      labels=c(0.1,0.3,0.5,0.7,0.9,1.0),
      include.lowest = FALSE,
      right=FALSE),
      levels=c(0.9,0.7,0.5,0.3,0.1,1.0)))


## COmbine GWAS table with lead SNP LD
cad_gwas_annot = merge(cad_gwas, ld_gwas_table,
                       by.x=c('panel_variant_id','chromosome','position'),
                       by.y=c('variant','chr','pos'))


## Read in tissue gene eQTLs
lipa_blood_eqtl = tabix.read.table("eQTL/Whole_Blood.allpairs.txt.gz",
                             "chr10:88000000-91000000") %>%
  filter(gene_id == "ENSG00000107798.17")

lipa_spleen_eqtl = tabix.read.table("eQTL/Spleen.allpairs.txt.gz",
                             "chr10:88000000-91000000") %>%
  filter(gene_id == "ENSG00000107798.17")

lipa_lung_eqtl = tabix.read.table("eQTL/Lung.allpairs.txt.gz",
                                    "chr10:88000000-91000000") %>%
  filter(gene_id == "ENSG00000107798.17")

lipa_artery_eqtl = tabix.read.table("eQTL/Artery_Aorta.allpairs.txt.gz",
                                  "chr10:88000000-91000000") %>%
  filter(gene_id == "ENSG00000107798.17")


## Read in ENLOC files
eQTL_enloc = do.call('rbind', lapply(c('Whole_Blood','Spleen','Lung','Artery_Aorta'), function(tissi) {
  enloc_i = read.table(paste0("ENLOC/enloc_results/",tissi,"_w_Coronary_Artery_Disease_enloc_output.txt.gz"),
                       header=TRUE,sep='\t') %>%
    filter(gene_id == "ENSG00000107798.17") %>%
    mutate(tiss=tissi) } ) )





## Merge CAD GWAS and LIPA Blood eQTL data
gwas_blood = merge(cad_gwas_annot, lipa_blood_eqtl,
                   by.x=c('panel_variant_id','chromosome','position'),
                   by.y=c('variant_id','chr','pos')) %>%
  rename(gwas_p = pvalue, eqtl_p = pval_nominal)

## Plot GWAS vs eQTL p-values
gwas_blood %>%
  ggplot() +
  geom_point(aes(-log10(gwas_p),-log10(eqtl_p),
                 col=r_2)) +
  theme_classic()

gwas_blood %>%
  arrange(r_2) %>% 
  ggplot(aes(-log10(gwas_p),-log10(eqtl_p))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  ylim(0,70) + xlim(0,13) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(variant_id), '')),
                  nudge_x=.1) +
  geom_text(x=3,y=60,
            label=paste('ENLOC rcp',max(eQTL_enloc %>% filter(tiss=='Whole_Blood') %>% pull(rcp)))) +
  theme(legend.position = 'none')



gwas_spleen = merge(cad_gwas_annot, lipa_spleen_eqtl,
                   by.x=c('panel_variant_id','chromosome','position'),
                   by.y=c('variant_id','chr','pos')) %>%
  rename(gwas_p = pvalue, eqtl_p = pval_nominal)
gwas_spleen %>%
  ggplot() +
  geom_point(aes(-log10(gwas_p),-log10(eqtl_p),
                 col=r_2)) +
  theme_classic()

gwas_spleen %>%
  arrange(r_2) %>% 
  ggplot(aes(-log10(gwas_p),-log10(eqtl_p))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  ylim(0,25) + xlim(0,13) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(variant_id), '')),
                  nudge_x=.1) +
  geom_text(x=3,y=20,
            label=paste('ENLOC rcp',max(eQTL_enloc %>% filter(tiss=='Spleen') %>% pull(rcp)))) +
  theme(legend.position = 'none')



gwas_lung = merge(cad_gwas_annot, lipa_lung_eqtl,
                   by.x=c('panel_variant_id','chromosome','position'),
                   by.y=c('variant_id','chr','pos')) %>%
  rename(gwas_p = pvalue, eqtl_p = pval_nominal)
gwas_lung %>%
  ggplot() +
  geom_point(aes(-log10(gwas_p),-log10(eqtl_p),
                 col=r_2)) +
  theme_classic()

gwas_lung %>%
  arrange(r_2) %>% 
  ggplot(aes(-log10(gwas_p),-log10(eqtl_p))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  ylim(0,NA) + xlim(0,13) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(variant_id), '')),
                  nudge_x=.1) +
  geom_text(x=3,y=20,
            label=paste('ENLOC rcp',max(eQTL_enloc %>% filter(tiss=='Lung') %>% pull(rcp)))) +
  theme(legend.position = 'none')



gwas_artery = merge(cad_gwas_annot, lipa_artery_eqtl,
                  by.x=c('panel_variant_id','chromosome','position'),
                  by.y=c('variant_id','chr','pos')) %>%
  rename(gwas_p = pvalue, eqtl_p = pval_nominal)
gwas_artery %>%
  ggplot() +
  geom_point(aes(-log10(gwas_p),-log10(eqtl_p),
                 col=r_2)) +
  theme_classic()

gwas_artery %>%
  arrange(r_2) %>% 
  ggplot(aes(-log10(gwas_p),-log10(eqtl_p))) +
  geom_point(aes(fill=r_2_cat,
                 shape=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'),
                              'lead','linked'))) +
  theme_classic() +
  ylim(0,NA) + xlim(0,13) +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(23,21)) +
  geom_text_repel(aes(label=ifelse(variant_id %in% c('rs1320496','rs1412444','rs1412445'), 
                                   as.character(variant_id), '')),
                  nudge_x=.1) +
  geom_text(x=3,y=7,
            label=paste('ENLOC rcp',max(eQTL_enloc %>% filter(tiss=='Artery_Aorta') %>% pull(rcp)))) +
  theme(legend.position = 'none')

