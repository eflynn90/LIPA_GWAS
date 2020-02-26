#!/usr/bin/Rscript
##
##
##  calcIndependentEffectsLIPA.R
##
##  EDF 12/17/18
##

library(gplots)
library(ggplot2)
library(reshape2)

setwd("~/projects/LIPA_GWAS/")


GTs_raw=read.table("GTs/all_LIPA_region_GTs.vcf",header = TRUE, comment.char = '', stringsAsFactors = FALSE)
eQTLs=read.table("eQTL/LIPA.eQTLs.sig_egenes",header=TRUE)
tiss_col = read.table("../EFF_EXPR/input_data/tissue_translation_colors_v8.txt",header=TRUE,sep='\t')


apply(GTs_raw, 1, function(row) {table(row[10:ncol(GTs_raw)])})

apply(GTs_raw, 2, function(col) {if ('./.' %in% col) {FALSE} else {TRUE}})

DS = as.data.frame( apply( GTs_raw[10:ncol(GTs_raw)], 1, function(row) {
                   sapply(row, function(el) { if (el=='./.') {NA}
                     else if (el=='0/0') {0} else if (el=='0/1') {1} else if (el=='1/1') {2} })
                 }) )
names(DS) <- t(GTs_raw['ID'])
DS$sample <- names(GTs_raw[10:ncol(GTs_raw)])

cor(DS[,-c(ncol(DS))], use='complete.obs')
heatmap.2(cor(DS[,-c(ncol(DS))], use='complete.obs')^2, symm=TRUE, 
        dendrogram='none', Rowv = FALSE, trace='none',margins = c(10,10))
heatmap.2(cor(DS[,-c(ncol(DS))], use='complete.obs')^2, symm=TRUE, 
          trace='none',margins = c(10,10))


## just select one SNP from each group to avoid singularity in LM (ugh, manually??)
SNPs_to_keep = c("chr10_89243047_C_T_b38",
                 "chr10_89249116_C_T_b38",
                 "chr10_89243088_C_T_b38",
                 "chr10_89245956_G_C_b38",
                 "chr10_89245251_C_T_b38",
                 "chr10_89252003_C_G_b38",
                 "chr10_89250704_A_T_b38",
                 "chr10_89246728_C_T_b38",
                 "chr10_89248248_G_A_b38",
                 "chr10_89250191_C_T_b38",
                 "chr10_89249904_A_T_b38",
                 "chr10_89244178_C_CA_b38"
)

DS_SNPs = DS[,names(DS) %in% SNPs_to_keep]
DS_SNPs$sample <- names(GTs_raw[10:ncol(GTs_raw)])


heatmap.2(cor(DS_SNPs[,-c(ncol(DS_SNPs))], use='complete.obs')^2, symm=TRUE, 
          trace='none',margins = c(10,10))

expr_files=list.files(path="expr/by_tiss/",full.names = TRUE)
lapply(expr_files, function(e_file) {
  
  print(e_file)
  
  expr=read.table(e_file,header=TRUE,comment.char='')
  expr_long=melt(expr,id.vars=c('X.chr','start','end','gene_id'))
  
  tiss_data = merge(DS_SNPs, expr_long, by.x='sample', by.y='variable')
  
  lin_mod = lm( value ~ ., tiss_data[c(2:(length(SNPs_to_keep)+1),ncol(tiss_data))] )
  
  print(summary(lin_mod))
  
  
}) -> trash





more_GTs_raw=read.table("GTs/larger_LIPA_region_GTs.vcf",header = TRUE, comment.char = '', stringsAsFactors = FALSE)

apply(more_GTs_raw, 1, function(row) {table(row[10:ncol(GTs_raw)])})

apply(more_GTs_raw, 2, function(col) {if ('./.' %in% col) {FALSE} else {TRUE}})

DS_more = as.data.frame( apply( more_GTs_raw[10:ncol(more_GTs_raw)], 1, function(row) {
  sapply(row, function(el) { if (el=='./.') {NA}
    else if (el=='0/0') {0} else if (el=='0/1') {1} else if (el=='1/1') {2} })
}) )
names(DS_more) <- t(more_GTs_raw['ID'])
DS_more$sample <- names(more_GTs_raw[10:ncol(more_GTs_raw)])

DS_more_corrs = cor(DS_more[-c(ncol(DS_more))], use='complete.obs', method='pearson')
DS_more_corrs_noNA = DS_more_corrs[-which(rowSums(is.na(DS_more_corrs)) == 2389),-which(rowSums(is.na(DS_more_corrs)) == 2389)]
DS_more_corrs_r2 = apply(DS_more_corrs_noNA,c(1,2), function(en){en^2})

which(row.names(DS_more_corrs_r2) == 'chr10_89243088_C_T_b38')
which(row.names(DS_more_corrs_r2) == 'chr10_89243047_C_T_b38')
row_labels = c( rep('',2171), 'GWAS SNP', 'PU.1 SNP', rep('',2214-2173))

heatmap.2(DS_more_corrs_r2, symm=TRUE, labRow = row_labels,
          dendrogram='none', Rowv = FALSE, trace='none')
heatmap.2(DS_more_corrs_r2, symm=TRUE, 
          trace='none', labRow = row_labels)

PU.1_linked = which(DS_more_corrs_r2[row.names(DS_more_corrs_r2) == 'chr10_89243088_C_T_b38',] > 0.1)
GWAS_linked = which(DS_more_corrs_r2[row.names(DS_more_corrs_r2) == 'chr10_89243047_C_T_b38',] > 0.1)
new_snp_set = unique(c(PU.1_linked, GWAS_linked))
DS_more_corrs_r2_linked_only = DS_more_corrs_r2[new_snp_set,new_snp_set]
which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243088_C_T_b38')
which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243047_C_T_b38')
row_labels_linked = c( rep('',150), 'GWAS SNP', 'PU.1 SNP', rep('',262-152))

heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none')
heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none',labRow = row_labels_linked)

eQTL_snps = unique(as.character(eQTLs$variant_id))
row_labels_eQTLs = sapply(row.names(DS_more_corrs_r2_linked_only), function(snp) { 
  if (snp %in% eQTL_snps) {snp} else {''}})
heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none',labRow = row_labels_eQTLs, margins = c(8,5))


PU.1_linked = which(DS_more_corrs_r2[row.names(DS_more_corrs_r2) == 'chr10_89243088_C_T_b38',] > 0.2)
GWAS_linked = which(DS_more_corrs_r2[row.names(DS_more_corrs_r2) == 'chr10_89243047_C_T_b38',] > 0.2)
new_snp_set = unique(c(PU.1_linked, GWAS_linked))
DS_more_corrs_r2_linked_only = DS_more_corrs_r2[new_snp_set,new_snp_set]
which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243088_C_T_b38')
which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243047_C_T_b38')
row_labels_linked = c( rep('',which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243047_C_T_b38')-1), 
                       'GWAS SNP', 'PU.1 SNP', 
                       rep('',nrow(DS_more_corrs_r2_linked_only)-which(row.names(DS_more_corrs_r2_linked_only) == 'chr10_89243088_C_T_b38')))

heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none')
heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none',labRow = row_labels_linked, margins = c(8,5))

eQTL_snps = unique(as.character(eQTLs$variant_id))
row_labels_eQTLs = sapply(row.names(DS_more_corrs_r2_linked_only), function(snp) { 
  if (snp %in% eQTL_snps) {snp} else {''}})
heatmap.2(DS_more_corrs_r2_linked_only, symm=TRUE, 
          trace='none',labRow = row_labels_eQTLs, margins = c(8,8))



## first just do lin regs with PU.1 SNP and eQTL SNP:
DS_SNPs = DS_more[,c('chr10_89243088_C_T_b38','chr10_89243047_C_T_b38')]
names(DS_SNPs) <- c('PU.1_snp','GWAS_snp')
DS_SNPs$sample <- row.names(DS_SNPs)
expr_files=list.files(path="expr/by_tiss/",full.names = TRUE)
as.data.frame(t(sapply(expr_files, function(e_file) {
  
  print(e_file)
  
  expr=read.table(e_file,header=TRUE,comment.char='')
  expr_long=melt(expr,id.vars=c('X.chr','start','end','gene_id'))
  
  tiss_data = merge(DS_SNPs, expr_long, by.x='sample', by.y='variable')
  
  lin_mod = lm( value ~ PU.1_snp * GWAS_snp, tiss_data)
  
  print(c(lin_mod$coefficients, summary(lin_mod)$coefficients[,4]))
}))) -> temp
names(temp) <- c('b_0','b_PU.1','b_GWAS','b_PU.1_GWAS','p_0','p_PU.1','p_GWAS','p_PU.1_GWAS')
row.names(temp) <- sapply(expr_files, function(file_name) {
  strsplit(strsplit(file_name,'/')[[1]][4],'[.]')[[1]][1]
})
temp$tiss_long <- sapply(expr_files, function(file_name) {
  strsplit(strsplit(file_name,'/')[[1]][4],'[.]')[[1]][1]
})

temp[which(temp$p_GWAS < 0.05),'tiss_long']
temp[which(temp$p_PU.1 < 0.05),'tiss_long']

temp$sig_info <- paste(temp$p_PU.1<0.05,temp$p_GWAS<0.05,sep='_')
temp$int_sig <- as.character(temp$p_PU.1_GWAS < 0.05)
ggplot(temp,aes(b_PU.1,b_GWAS)) + 
  geom_point(aes(col=sig_info,shape=int_sig)) + 
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)

temp_cols = merge(temp,tiss_col,by.x='tiss_long',by.y='TISSUE_NAME')
ggplot(temp_cols,aes(b_PU.1,b_GWAS)) + 
  geom_point(aes(col=as.character(TISSUE_ABBRV),shape=as.character(TISSUE_ABBRV),fill=as.character(TISSUE_ABBRV))) + 
  scale_color_manual(values=as.character(temp_cols[order(temp_cols$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(temp_cols[order(temp_cols$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7))) +
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)




