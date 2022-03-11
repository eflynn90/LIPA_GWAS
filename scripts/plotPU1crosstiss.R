#!/usr/bin/Rscript
##
##  plotPU1crosstiss.R
##
##  EDF 3/12/2021
##

library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(tibble)

setwd("~/projects/LIPA_GWAS/")


tiss_colors = read.table("tissue_translation_colors_v8.txt",
                         header=TRUE,sep='\t')


spi1_expr = read.table("expr/ENSG00000066336.rnaseqc.median.txt",
                       header=TRUE,sep='\t') %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('tiss') %>%
  slice(3:n())
names(spi1_expr) <- c('tiss','med_tpm')

stat1_expr = read.table("expr/ENSG00000115415.rnaseqc.median.txt",
                       header=TRUE,sep='\t') %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('tiss') %>%
  slice(3:n())
names(stat1_expr) <- c('tiss','med_tpm')

lipa_afc = read.table("afcs/LIPA.afcs.cov.txt",
                      header=TRUE,sep='\t')

lipa_rs096_afc = lipa_afc %>%
  filter(pos==89243088) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('tiss') %>%
  slice(5:n())
names(lipa_rs096_afc) <- c('tiss','afc')

both_rs096 = merge(lipa_rs096_afc, spi1_expr) %>%
  merge(tiss_colors, by.x='tiss', by.y='TISSUE_ABBRV') %>%
  mutate(afc = as.numeric(as.character(afc)),
         med_tpm = as.numeric(as.character(med_tpm)))

sp_rho = cor.test(both_rs096$med_tpm, both_rs096$afc,
                  method='spearman')$estimate
sp_p = cor.test(both_rs096$med_tpm, both_rs096$afc,
                  method='spearman')$p.value

both_rs096 %>%
  ggplot(aes(log10(med_tpm),afc)) +
  geom_hline(yintercept=0) +
  geom_point(col=both_rs096$TISSUE_RCOL) +
  geom_text_repel(aes(label=ifelse(tiss=='WHLBLD','Blood',
                                   ifelse(tiss=='SPLEEN','Spleen',
                                          ifelse(tiss=='LUNG','Lung',''))))) +
  theme_classic() +
  xlab('SPI1 level\nlog10(median TPM)') +
  ylab('rs1320496 effect on LIPA expression\nlog2 allelic fold change') +
  geom_label(aes(x=-.5,y=0.65),
             label=paste0('rho=',round(sp_rho,2),'\np=',round(sp_p,4)),
             hjust=0)



lipa_rs444_afc = lipa_afc %>%
  filter(pos==89243047) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('tiss') %>%
  slice(5:n())
names(lipa_rs444_afc) <- c('tiss','afc')

both_rs444 = merge(lipa_rs444_afc, spi1_expr) %>%
  merge(tiss_colors, by.x='tiss', by.y='TISSUE_ABBRV') %>%
  mutate(afc = as.numeric(as.character(afc)),
         med_tpm = as.numeric(as.character(med_tpm)))

both_rs444 %>%
  ggplot(aes(log10(med_tpm),afc)) +
  geom_hline(yintercept=0) +
  geom_point(col=both_rs444$TISSUE_RCOL) +
  geom_text_repel(aes(label=ifelse(tiss=='WHLBLD','Blood',
                                   ifelse(tiss=='SPLEEN','Spleen',
                                          ifelse(tiss=='LUNG','Lung',''))))) +
  theme_classic() +
  xlab('SPI1 level\nlog10(median TPM)') +
  ylab('rs1412444 effect on LIPA expression\nlog2 allelic fold change')






both_rs444_stat = merge(lipa_rs444_afc, stat1_expr) %>%
  merge(tiss_colors, by.x='tiss', by.y='TISSUE_ABBRV') %>%
  mutate(afc = as.numeric(as.character(afc)),
         med_tpm = as.numeric(as.character(med_tpm)))

sp_rho = cor.test(both_rs444_stat$med_tpm, both_rs444_stat$afc,
                  method='spearman')$estimate
sp_p = cor.test(both_rs444_stat$med_tpm, both_rs444_stat$afc,
                method='spearman')$p.value

both_rs444_stat %>%
  ggplot(aes(log10(med_tpm),afc)) +
  geom_hline(yintercept=0) +
  geom_point(col=both_rs444$TISSUE_RCOL) +
  geom_text_repel(aes(label=ifelse(tiss=='WHLBLD','Blood',
                                   ifelse(tiss=='SPLEEN','Spleen',
                                          ifelse(tiss=='LUNG','Lung',''))))) +
  theme_classic() +
  xlab('STAT1 level\nlog10(median TPM)') +
  ylab('rs1412444 effect on LIPA expression\nlog2 allelic fold change') +
  geom_label(aes(x=.8,y=0.65),
             label=paste0('rho=',round(sp_rho,2),'\np=',round(sp_p,2)),
             hjust=0)






