#!/usr/bin/Rscript
##
##
##  exploreLIPAData.R
##
##  EDF 12/17/18
##

library(reshape2)
library(ggplot2)

setwd("~/projects/LIPA_GWAS/")
GTs_raw=read.table("GTs/GTs.vcf",header = TRUE, comment.char = '', stringsAsFactors = FALSE)
table(unlist(GTs_raw[1,10:ncol(GTs_raw)]))
table(unlist(GTs_raw[2,10:ncol(GTs_raw)]))
table(unlist(GTs_raw[1,10:ncol(GTs_raw)]),unlist(GTs_raw[2,10:ncol(GTs_raw)]))
GTs = data.frame(names(GTs_raw[10:ncol(GTs_raw)]),
                 t(GTs_raw[2,10:ncol(GTs_raw)]),
                 t(GTs_raw[1,10:ncol(GTs_raw)]),
                 apply(GTs_raw[,10:ncol(GTs_raw)], 2, function(col) {paste(col[2],col[1],collapse='',sep='_')}))
names(GTs) = c('sample','PU.1_SNP_GT','GWAS_SNP_GT','GTs')
GTs$sample <- as.character(GTs$sample)



WHLBLD_expr=read.table("expr/by_tiss/Whole_Blood.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
WHLBLD_expr_long=melt(WHLBLD_expr,id.vars = c('X.chr','start','end','gene_id'))

WHLBLD_data = merge(GTs,WHLBLD_expr_long,by.x='sample',by.y='variable')
WHLBLD_data = WHLBLD_data[ ! (WHLBLD_data$GTs %in% c("./._0/0","1/1_./.") ),]

table(WHLBLD_data$PU.1_SNP_GT,WHLBLD_data$GWAS_SNP_GT)

WHLBLD_data$GTs_PU.1snp <- factor(WHLBLD_data$GTs, 
                                  levels=c('0/0_0/0',
                                           '0/1_0/0','0/1_0/1',
                                           '1/1_0/0','1/1_0/1','1/1_1/1'))
boxplot(WHLBLD_data$value ~ WHLBLD_data$GTs_PU.1snp)
## note: log2 expr

WHLBLD_data$GTs_GWASsnp <- factor(WHLBLD_data$GTs,
                                  levels=c('0/0_0/0','0/1_0/0','1/1_0/0',
                                           '0/1_0/1','1/1_0/1',
                                           '1/1_1/1'))
boxplot(WHLBLD_data$value ~ WHLBLD_data$GTs_GWASsnp)
## note: log2 expr

wilcox.test(WHLBLD_data[WHLBLD_data$GTs == '0/1_0/1',]$value, 
            WHLBLD_data[WHLBLD_data$GTs == '1/1_0/1',]$value)
t.test(WHLBLD_data[WHLBLD_data$GTs == '0/1_0/1',]$value, 
       WHLBLD_data[WHLBLD_data$GTs == '1/1_0/1',]$value)

t.test(WHLBLD_data[WHLBLD_data$GTs == '0/0_0/0',]$value, 
       WHLBLD_data[WHLBLD_data$GTs == '0/1_0/0',]$value)
t.test(WHLBLD_data[WHLBLD_data$GTs == '0/1_0/0',]$value, 
       WHLBLD_data[WHLBLD_data$GTs == '1/1_0/0',]$value)
t.test(WHLBLD_data[WHLBLD_data$GTs == '0/0_0/0',]$value, 
       WHLBLD_data[WHLBLD_data$GTs == '1/1_0/0',]$value)


WHLBLD_data$PU.1_SNP_DS <- '0'
WHLBLD_data[WHLBLD_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- '1'
WHLBLD_data[WHLBLD_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- '2'

WHLBLD_data$GWAS_SNP_DS <- '0'
WHLBLD_data[WHLBLD_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- '1'
WHLBLD_data[WHLBLD_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- '2'


summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) + as.numeric(WHLBLD_data$GWAS_SNP_DS)))
summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) * as.numeric(WHLBLD_data$GWAS_SNP_DS)))

summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$GWAS_SNP_DS)))
summary(lm( WHLBLD_data$value - (as.numeric(WHLBLD_data$GWAS_SNP_DS) * 0.32665) ~
              as.numeric(WHLBLD_data$PU.1_SNP_DS)))

summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS)))
summary(lm( WHLBLD_data$value - (as.numeric(WHLBLD_data$PU.1_SNP_DS) * 0.29413) ~
              as.numeric(WHLBLD_data$GWAS_SNP_DS)))

summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$GWAS_SNP_DS) | as.numeric(WHLBLD_data$PU.1_SNP_DS)))
summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) | as.numeric(WHLBLD_data$GWAS_SNP_DS)))





SPLEEN_expr=read.table("expr/by_tiss/Spleen.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
SPLEEN_expr_long=melt(SPLEEN_expr,id.vars = c('X.chr','start','end','gene_id'))

SPLEEN_data = merge(GTs,SPLEEN_expr_long,by.x='sample',by.y='variable')
SPLEEN_data = SPLEEN_data[ ! (SPLEEN_data$GTs %in% c("./._0/0","1/1_./.") ),]

table(SPLEEN_data$PU.1_SNP_GT,SPLEEN_data$GWAS_SNP_GT)

SPLEEN_data$GTs_PU.1snp <- factor(SPLEEN_data$GTs, 
                                  levels=c('0/0_0/0',
                                           '0/1_0/0','0/1_0/1',
                                           '1/1_0/0','1/1_0/1','1/1_1/1'))
boxplot(SPLEEN_data$value ~ SPLEEN_data$GTs_PU.1snp)
## note: log2 expr

SPLEEN_data$GTs_GWASsnp <- factor(SPLEEN_data$GTs,
                                  levels=c('0/0_0/0','0/1_0/0','1/1_0/0',
                                           '0/1_0/1','1/1_0/1',
                                           '1/1_1/1'))
boxplot(SPLEEN_data$value ~ SPLEEN_data$GTs_GWASsnp)
## note: log2 expr


SPLEEN_data$PU.1_SNP_DS <- '0'
SPLEEN_data[SPLEEN_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- '1'
SPLEEN_data[SPLEEN_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- '2'

SPLEEN_data$GWAS_SNP_DS <- '0'
SPLEEN_data[SPLEEN_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- '1'
SPLEEN_data[SPLEEN_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- '2'


summary(lm(SPLEEN_data$value ~ as.numeric(SPLEEN_data$PU.1_SNP_DS) * as.numeric(SPLEEN_data$GWAS_SNP_DS)))
summary(lm(SPLEEN_data$value ~ as.numeric(SPLEEN_data$PU.1_SNP_DS) + as.numeric(SPLEEN_data$GWAS_SNP_DS)))




LCL_expr=read.table("expr/by_tiss/Cells_EBV-transformed_lymphocytes.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
LCL_expr_long=melt(LCL_expr,id.vars = c('X.chr','start','end','gene_id'))

LCL_data = merge(GTs,LCL_expr_long,by.x='sample',by.y='variable')
LCL_data = LCL_data[ ! (LCL_data$GTs %in% c("./._0/0","1/1_./.") ),]

table(LCL_data$PU.1_SNP_GT,LCL_data$GWAS_SNP_GT)

LCL_data$GTs_PU.1snp <- factor(LCL_data$GTs, 
                                  levels=c('0/0_0/0',
                                           '0/1_0/0','0/1_0/1',
                                           '1/1_0/0','1/1_0/1','1/1_1/1'))
boxplot(LCL_data$value ~ LCL_data$GTs_PU.1snp)
## note: log2 expr

LCL_data$GTs_GWASsnp <- factor(LCL_data$GTs,
                                  levels=c('0/0_0/0','0/1_0/0','1/1_0/0',
                                           '0/1_0/1','1/1_0/1',
                                           '1/1_1/1'))
boxplot(LCL_data$value ~ LCL_data$GTs_GWASsnp)
## note: log2 expr


LCL_data$PU.1_SNP_DS <- '0'
LCL_data[LCL_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- '1'
LCL_data[LCL_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- '2'

LCL_data$GWAS_SNP_DS <- '0'
LCL_data[LCL_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- '1'
LCL_data[LCL_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- '2'


summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS) * as.numeric(LCL_data$GWAS_SNP_DS)))
summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS) + as.numeric(LCL_data$GWAS_SNP_DS)))




FIBRBLS_expr=read.table("FIBRBLS.expr.txt",header=TRUE, comment.char = '')
FIBRBLS_expr_long=melt(FIBRBLS_expr,id.vars = c('X.chr','start','end','gene_id'))

FIBRBLS_data = merge(GTs,FIBRBLS_expr_long,by.x='sample',by.y='variable')
FIBRBLS_data = FIBRBLS_data[ ! (FIBRBLS_data$GTs %in% c("./._0/0","1/1_./.") ),]

table(FIBRBLS_data$PU.1_SNP_GT,FIBRBLS_data$GWAS_SNP_GT)

FIBRBLS_data$GTs_PU.1snp <- factor(FIBRBLS_data$GTs, 
                               levels=c('0/0_0/0',
                                        '0/1_0/0','0/1_0/1',
                                        '1/1_0/0','1/1_0/1','1/1_1/1'))
boxplot(FIBRBLS_data$value ~ FIBRBLS_data$GTs_PU.1snp)
## note: log2 expr

FIBRBLS_data$GTs_GWASsnp <- factor(FIBRBLS_data$GTs,
                               levels=c('0/0_0/0','0/1_0/0','1/1_0/0',
                                        '0/1_0/1','1/1_0/1',
                                        '1/1_1/1'))
boxplot(FIBRBLS_data$value ~ FIBRBLS_data$GTs_GWASsnp)
## note: log2 expr


FIBRBLS_data$PU.1_SNP_DS <- '0'
FIBRBLS_data[FIBRBLS_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- '1'
FIBRBLS_data[FIBRBLS_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- '2'

FIBRBLS_data$GWAS_SNP_DS <- '0'
FIBRBLS_data[FIBRBLS_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- '1'
FIBRBLS_data[FIBRBLS_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- '2'


summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$PU.1_SNP_DS) * as.numeric(FIBRBLS_data$GWAS_SNP_DS)))





LUNG_expr=read.table("LUNG.expr.txt",header=TRUE, comment.char = '')
LUNG_expr_long=melt(LUNG_expr,id.vars = c('X.chr','start','end','gene_id'))

LUNG_data = merge(GTs,LUNG_expr_long,by.x='sample',by.y='variable')
LUNG_data = LUNG_data[ ! (LUNG_data$GTs %in% c("./._0/0","1/1_./.") ),]

table(LUNG_data$PU.1_SNP_GT,LUNG_data$GWAS_SNP_GT)

LUNG_data$GTs_PU.1snp <- factor(LUNG_data$GTs, 
                               levels=c('0/0_0/0',
                                        '0/1_0/0','0/1_0/1',
                                        '1/1_0/0','1/1_0/1','1/1_1/1'))
boxplot(LUNG_data$value ~ LUNG_data$GTs_PU.1snp)
## note: log2 expr

LUNG_data$GTs_GWASsnp <- factor(LUNG_data$GTs,
                               levels=c('0/0_0/0','0/1_0/0','1/1_0/0',
                                        '0/1_0/1','1/1_0/1',
                                        '1/1_1/1'))
boxplot(LUNG_data$value ~ LUNG_data$GTs_GWASsnp)
## note: log2 expr


LUNG_data$PU.1_SNP_DS <- '0'
LUNG_data[LUNG_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- '1'
LUNG_data[LUNG_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- '2'

LUNG_data$GWAS_SNP_DS <- '0'
LUNG_data[LUNG_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- '1'
LUNG_data[LUNG_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- '2'


summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$PU.1_SNP_DS) * as.numeric(LUNG_data$GWAS_SNP_DS)))



tiss_trans=read.table('tissue_translation_colors_v8.txt',header=TRUE,sep='\t')
aFCs = read.table("aFCs_topeQTLs.txt",header=TRUE)
this_aFC=t(aFCs[aFCs$gene_id == 'ENSG00000107798.17',-c(1,2,3,4,54,55)])
row.names(this_aFC) <- sub('[.]','-',row.names(this_aFC))
these_aFCs = merge(tiss_trans, this_aFC,
                   by.x='TISSUE_NAME',by.y='row.names')
names(these_aFCs) <- c(names(these_aFCs[-c(6)]),'aFC')
cell_types=t(read.table("xCell_tpm.v8_res_all.tissues_res.m.txt",header=TRUE))
macro=merge(these_aFCs, cell_types,
            by.x='TISSUE_ABBRV',by.y='row.names')

ggplot(macro, aes(Macrophages,aFC)) +
  geom_hline(yintercept = 0) + 
  geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
  theme_classic() +
  scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))
ggplot(macro, aes(Macrophages_M1,aFC)) +
  geom_hline(yintercept = 0) + 
  geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
  theme_classic() +
  scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))
ggplot(macro, aes(Macrophages_M2,aFC)) +
  geom_hline(yintercept = 0) + 
  geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
  theme_classic() +
  scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))
ggplot(macro, aes(Monocytes,aFC)) +
  geom_hline(yintercept = 0) + 
  geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
  theme_classic() +
  scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))

apply(macro[c('Macrophages','Macrophages_M1','Macrophages_M2','Monocytes','B_cells','Csm_B_cells','Memory_B_cells','naive_B_cells','pro_B_cells')], 2, function(col) {
  print('NEW')
  print(cor.test(col,macro$aFC,method='spearman'))
  print(cor.test(col[order(macro$TISSUE_ABBRV)],
                 as.numeric(TF_expr[TF_expr$TF == 'PU1',which(names(TF_expr) %in% macro$TISSUE_ABBRV)]),
                 method='spearman'))
})

apply(macro[c('Macrophages','Macrophages_M1','Macrophages_M2','Monocytes')], 1, sum) -> macro$sum_macros
ggplot(macro, aes(sum_macros,aFC)) +
  geom_hline(yintercept = 0) + 
  geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
  theme_classic() +
  scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
  scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))
cor.test(macro$sum_macros,macro$aFC,method='spearman')
cor.test(macro[order(macro$TISSUE_ABBRV),]$sum_macros,
         as.numeric(TF_expr[TF_expr$TF == 'PU1',which(names(TF_expr) %in% macro$TISSUE_ABBRV)]),
         method='spearman')

            