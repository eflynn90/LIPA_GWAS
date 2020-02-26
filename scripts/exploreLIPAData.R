#!/usr/bin/Rscript
##
##
##  exploreLIPAData.R
##
##  EDF 12/17/18
##

library(reshape2)
library(ggplot2)
library(dplyr)

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
GTs$PU.1_SNP_DS <- 0
GTs[GTs$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
GTs[GTs$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

GTs$GWAS_SNP_DS <- 0
GTs[GTs$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
GTs[GTs$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2

tiss_trans = read.table("tissue_translation_colors_v8.txt",header=TRUE,sep='\t')


WHLBLD_expr=read.table("expr/by_tiss/Whole_Blood.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
WHLBLD_expr_long=melt(WHLBLD_expr,id.vars = c('X.chr','start','end','gene_id'))

WHLBLD_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt",
                                        header=TRUE))[-1,])
WHLBLD_covar = as.data.frame(apply(WHLBLD_covar, c(1,2), function(x) {as.numeric(x)}))
names(WHLBLD_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt",
                                  header=TRUE)[1])

WHLBLD_data = merge(GTs,WHLBLD_expr_long,by.x='sample',by.y='variable')
WHLBLD_data = WHLBLD_data[ ! (WHLBLD_data$GTs %in% c("./._0/0","1/1_./.") ),]

WHLBLD_data = merge(WHLBLD_data, WHLBLD_covar,
                    by.x='sample', by.y='row.names')

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


WHLBLD_data$PU.1_SNP_DS <- 0
WHLBLD_data[WHLBLD_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
WHLBLD_data[WHLBLD_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

WHLBLD_data$GWAS_SNP_DS <- 0
WHLBLD_data[WHLBLD_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
WHLBLD_data[WHLBLD_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2


summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) + as.numeric(WHLBLD_data$GWAS_SNP_DS)))
summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) * as.numeric(WHLBLD_data$GWAS_SNP_DS)))

summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$GWAS_SNP_DS)))
summary(lm( WHLBLD_data$value - (as.numeric(WHLBLD_data$GWAS_SNP_DS) * 0.32665) ~
              as.numeric(WHLBLD_data$PU.1_SNP_DS)))

summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS)))
summary(lm( WHLBLD_data$value - (as.numeric(WHLBLD_data$PU.1_SNP_DS) * 0.29413) ~
              as.numeric(WHLBLD_data$GWAS_SNP_DS)))

# summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$GWAS_SNP_DS) | as.numeric(WHLBLD_data$PU.1_SNP_DS)))
# summary(lm(WHLBLD_data$value ~ as.numeric(WHLBLD_data$PU.1_SNP_DS) | as.numeric(WHLBLD_data$GWAS_SNP_DS)))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))
summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS * GWAS_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))
summary(lm( as.formula( paste0( "value - 0.361178 * PU.1_SNP_DS ~ GWAS_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))

summary(lm( as.formula( paste0( "value ~ GWAS_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))
summary(lm( as.formula( paste0( "value - 0.383446 * GWAS_SNP_DS ~ PU.1_SNP_DS +",
                                paste(names(WHLBLD_data)[10:77],collapse='+')) ),
            data=WHLBLD_data))






SPLEEN_expr=read.table("expr/by_tiss/Spleen.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
SPLEEN_expr_long=melt(SPLEEN_expr,id.vars = c('X.chr','start','end','gene_id'))

SPLEEN_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Spleen.v8.covariates.txt",
                                        header=TRUE))[-1,])
SPLEEN_covar = as.data.frame(apply(SPLEEN_covar, c(1,2), function(x) {as.numeric(x)}))
names(SPLEEN_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Spleen.v8.covariates.txt",
                                         header=TRUE)[1])

SPLEEN_data = merge(GTs,SPLEEN_expr_long,by.x='sample',by.y='variable')
SPLEEN_data = SPLEEN_data[ ! (SPLEEN_data$GTs %in% c("./._0/0","1/1_./.") ),]

SPLEEN_data = merge(SPLEEN_data, SPLEEN_covar,
                    by.x='sample', by.y='row.names')

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


SPLEEN_data$PU.1_SNP_DS <- 0
SPLEEN_data[SPLEEN_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
SPLEEN_data[SPLEEN_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

SPLEEN_data$GWAS_SNP_DS <- 0
SPLEEN_data[SPLEEN_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
SPLEEN_data[SPLEEN_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2


summary(lm(SPLEEN_data$value ~ as.numeric(SPLEEN_data$PU.1_SNP_DS) * as.numeric(SPLEEN_data$GWAS_SNP_DS)))
summary(lm(SPLEEN_data$value ~ as.numeric(SPLEEN_data$PU.1_SNP_DS) + as.numeric(SPLEEN_data$GWAS_SNP_DS)))



summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))
summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS * GWAS_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))
summary(lm( as.formula( paste0( "value - 0.24604 * PU.1_SNP_DS ~ GWAS_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))

summary(lm( as.formula( paste0( "value ~ GWAS_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))
summary(lm( as.formula( paste0( "value - 0.25287 * GWAS_SNP_DS ~ PU.1_SNP_DS +",
                                paste(names(SPLEEN_data)[10:44],collapse='+')) ),
            data=SPLEEN_data))



LCL_expr=read.table("expr/by_tiss/Cells_EBV-transformed_lymphocytes.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
LCL_expr_long=melt(LCL_expr,id.vars = c('X.chr','start','end','gene_id'))

LCL_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_EBV-transformed_lymphocytes.v8.covariates.txt",
                                        header=TRUE))[-1,])
LCL_covar = as.data.frame(apply(LCL_covar, c(1,2), function(x) {as.numeric(x)}))
names(LCL_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_EBV-transformed_lymphocytes.v8.covariates.txt",
                                         header=TRUE)[1])

LCL_data = merge(GTs,LCL_expr_long,by.x='sample',by.y='variable')
LCL_data = LCL_data[ ! (LCL_data$GTs %in% c("./._0/0","1/1_./.") ),]

LCL_data = merge(LCL_data, LCL_covar,
                    by.x='sample', by.y='row.names')

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

wilcox.test(LCL_data[LCL_data$GTs == '0/1_0/1',]$value, 
            LCL_data[LCL_data$GTs == '1/1_0/1',]$value)
t.test(LCL_data[LCL_data$GTs == '0/1_0/1',]$value, 
       LCL_data[LCL_data$GTs == '1/1_0/1',]$value)

t.test(LCL_data[LCL_data$GTs == '0/0_0/0',]$value, 
       LCL_data[LCL_data$GTs == '0/1_0/0',]$value)
t.test(LCL_data[LCL_data$GTs == '0/1_0/0',]$value, 
       LCL_data[LCL_data$GTs == '1/1_0/0',]$value)
t.test(LCL_data[LCL_data$GTs == '0/0_0/0',]$value, 
       LCL_data[LCL_data$GTs == '1/1_0/0',]$value)


LCL_data$PU.1_SNP_DS <- 0
LCL_data[LCL_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
LCL_data[LCL_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

LCL_data$GWAS_SNP_DS <- 0
LCL_data[LCL_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
LCL_data[LCL_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2


summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS) + as.numeric(LCL_data$GWAS_SNP_DS)))
summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS) * as.numeric(LCL_data$GWAS_SNP_DS)))

summary(lm(LCL_data$value ~ as.numeric(LCL_data$GWAS_SNP_DS)))
summary(lm( LCL_data$value - (as.numeric(LCL_data$GWAS_SNP_DS) * -0.10678) ~
              as.numeric(LCL_data$PU.1_SNP_DS)))

summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS)))
summary(lm( LCL_data$value - (as.numeric(LCL_data$PU.1_SNP_DS) * 0.05192) ~
              as.numeric(LCL_data$GWAS_SNP_DS)))

# summary(lm(LCL_data$value ~ as.numeric(LCL_data$GWAS_SNP_DS) | as.numeric(LCL_data$PU.1_SNP_DS)))
# summary(lm(LCL_data$value ~ as.numeric(LCL_data$PU.1_SNP_DS) | as.numeric(LCL_data$GWAS_SNP_DS)))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))
summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS * GWAS_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))
summary(lm( as.formula( paste0( "value - 0.0695292 * PU.1_SNP_DS ~ GWAS_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))

summary(lm( as.formula( paste0( "value ~ GWAS_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))
summary(lm( as.formula( paste0( "value - -0.061606 * GWAS_SNP_DS ~ PU.1_SNP_DS +",
                                paste(names(LCL_data)[10:32],collapse='+')) ),
            data=LCL_data))






#FIBRBLS_expr=read.table("FIBRBLS.expr.txt",header=TRUE, comment.char = '')
FIBRBLS_expr=read.table("expr/by_tiss/Cells_Cultured_fibroblasts.v8.deseq_log2_expression.ENSG00000107798.17.bed",
                        header=TRUE, comment.char = '')
FIBRBLS_expr_long=melt(FIBRBLS_expr,id.vars = c('X.chr','start','end','gene_id'))

FIBRBLS_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_Cultured_fibroblasts.v8.covariates.txt",
                                        header=TRUE))[-1,])
FIBRBLS_covar = as.data.frame(apply(FIBRBLS_covar, c(1,2), function(x) {as.numeric(x)}))
names(FIBRBLS_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_Cultured_fibroblasts.v8.covariates.txt",
                                         header=TRUE)[1])

FIBRBLS_data = merge(GTs,FIBRBLS_expr_long,by.x='sample',by.y='variable')
FIBRBLS_data = FIBRBLS_data[ ! (FIBRBLS_data$GTs %in% c("./._0/0","1/1_./.") ),]

FIBRBLS_data = merge(FIBRBLS_data, FIBRBLS_covar,
                    by.x='sample', by.y='row.names')

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

wilcox.test(FIBRBLS_data[FIBRBLS_data$GTs == '0/1_0/1',]$value, 
            FIBRBLS_data[FIBRBLS_data$GTs == '1/1_0/1',]$value)
t.test(FIBRBLS_data[FIBRBLS_data$GTs == '0/1_0/1',]$value, 
       FIBRBLS_data[FIBRBLS_data$GTs == '1/1_0/1',]$value)

t.test(FIBRBLS_data[FIBRBLS_data$GTs == '0/0_0/0',]$value, 
       FIBRBLS_data[FIBRBLS_data$GTs == '0/1_0/0',]$value)
t.test(FIBRBLS_data[FIBRBLS_data$GTs == '0/1_0/0',]$value, 
       FIBRBLS_data[FIBRBLS_data$GTs == '1/1_0/0',]$value)
t.test(FIBRBLS_data[FIBRBLS_data$GTs == '0/0_0/0',]$value, 
       FIBRBLS_data[FIBRBLS_data$GTs == '1/1_0/0',]$value)


FIBRBLS_data$PU.1_SNP_DS <- 0
FIBRBLS_data[FIBRBLS_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
FIBRBLS_data[FIBRBLS_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

FIBRBLS_data$GWAS_SNP_DS <- 0
FIBRBLS_data[FIBRBLS_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
FIBRBLS_data[FIBRBLS_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2


summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$PU.1_SNP_DS) + as.numeric(FIBRBLS_data$GWAS_SNP_DS)))
summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$PU.1_SNP_DS) * as.numeric(FIBRBLS_data$GWAS_SNP_DS)))

summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$GWAS_SNP_DS)))
summary(lm( FIBRBLS_data$value - (as.numeric(FIBRBLS_data$GWAS_SNP_DS) * -0.01810) ~
              as.numeric(FIBRBLS_data$PU.1_SNP_DS)))

summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$PU.1_SNP_DS)))
summary(lm( FIBRBLS_data$value - (as.numeric(FIBRBLS_data$PU.1_SNP_DS) * -0.01654) ~
              as.numeric(FIBRBLS_data$GWAS_SNP_DS)))

# summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$GWAS_SNP_DS) | as.numeric(FIBRBLS_data$PU.1_SNP_DS)))
# summary(lm(FIBRBLS_data$value ~ as.numeric(FIBRBLS_data$PU.1_SNP_DS) | as.numeric(FIBRBLS_data$GWAS_SNP_DS)))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))
summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS * GWAS_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))
summary(lm( as.formula( paste0( "value - -9.011e-03 * PU.1_SNP_DS ~ GWAS_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))

summary(lm( as.formula( paste0( "value ~ GWAS_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))
summary(lm( as.formula( paste0( "value - -2.892e-04 * GWAS_SNP_DS ~ PU.1_SNP_DS +",
                                paste(names(FIBRBLS_data)[10:77],collapse='+')) ),
            data=FIBRBLS_data))






#LUNG_expr=read.table("LUNG.expr.txt",header=TRUE, comment.char = '')
LUNG_expr=read.table("expr/by_tiss/Lung.v8.deseq_log2_expression.ENSG00000107798.17.bed",
                     header=TRUE, comment.char = '')
LUNG_expr_long=melt(LUNG_expr,id.vars = c('X.chr','start','end','gene_id'))

LUNG_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt",
                                        header=TRUE))[-1,])
LUNG_covar = as.data.frame(apply(LUNG_covar, c(1,2), function(x) {as.numeric(x)}))
names(LUNG_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt",
                                         header=TRUE)[1])

LUNG_data = merge(GTs,LUNG_expr_long,by.x='sample',by.y='variable')
LUNG_data = LUNG_data[ ! (LUNG_data$GTs %in% c("./._0/0","1/1_./.") ),]

LUNG_data = merge(LUNG_data, LUNG_covar,
                    by.x='sample', by.y='row.names')

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

wilcox.test(LUNG_data[LUNG_data$GTs == '0/1_0/1',]$value, 
            LUNG_data[LUNG_data$GTs == '1/1_0/1',]$value)
t.test(LUNG_data[LUNG_data$GTs == '0/1_0/1',]$value, 
       LUNG_data[LUNG_data$GTs == '1/1_0/1',]$value)

t.test(LUNG_data[LUNG_data$GTs == '0/0_0/0',]$value, 
       LUNG_data[LUNG_data$GTs == '0/1_0/0',]$value)
t.test(LUNG_data[LUNG_data$GTs == '0/1_0/0',]$value, 
       LUNG_data[LUNG_data$GTs == '1/1_0/0',]$value)
t.test(LUNG_data[LUNG_data$GTs == '0/0_0/0',]$value, 
       LUNG_data[LUNG_data$GTs == '1/1_0/0',]$value)


LUNG_data$PU.1_SNP_DS <- 0
LUNG_data[LUNG_data$PU.1_SNP_GT == '0/1',]$PU.1_SNP_DS <- 1
LUNG_data[LUNG_data$PU.1_SNP_GT == '1/1',]$PU.1_SNP_DS <- 2

LUNG_data$GWAS_SNP_DS <- 0
LUNG_data[LUNG_data$GWAS_SNP_GT == '0/1',]$GWAS_SNP_DS <- 1
LUNG_data[LUNG_data$GWAS_SNP_GT == '1/1',]$GWAS_SNP_DS <- 2


summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$PU.1_SNP_DS) + as.numeric(LUNG_data$GWAS_SNP_DS)))
summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$PU.1_SNP_DS) * as.numeric(LUNG_data$GWAS_SNP_DS)))

summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$GWAS_SNP_DS)))
summary(lm( LUNG_data$value - (as.numeric(LUNG_data$GWAS_SNP_DS) * 0.15606) ~
              as.numeric(LUNG_data$PU.1_SNP_DS)))

summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$PU.1_SNP_DS)))
summary(lm( LUNG_data$value - (as.numeric(LUNG_data$PU.1_SNP_DS) * 0.13054) ~
              as.numeric(LUNG_data$GWAS_SNP_DS)))

# summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$GWAS_SNP_DS) | as.numeric(LUNG_data$PU.1_SNP_DS)))
# summary(lm(LUNG_data$value ~ as.numeric(LUNG_data$PU.1_SNP_DS) | as.numeric(LUNG_data$GWAS_SNP_DS)))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))
summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS * GWAS_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))

summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))
summary(lm( as.formula( paste0( "value - 0.14955 * PU.1_SNP_DS ~ GWAS_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))

summary(lm( as.formula( paste0( "value ~ GWAS_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))
summary(lm( as.formula( paste0( "value - 0.222388 * GWAS_SNP_DS ~ PU.1_SNP_DS +",
                                paste(names(LUNG_data)[10:77],collapse='+')) ),
            data=LUNG_data))





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

# apply(macro[c('Macrophages','Macrophages_M1','Macrophages_M2','Monocytes','B_cells','Csm_B_cells','Memory_B_cells','naive_B_cells','pro_B_cells')], 2, function(col) {
#   print('NEW')
#   print(cor.test(col,macro$aFC,method='spearman'))
#   print(cor.test(col[order(macro$TISSUE_ABBRV)],
#                  as.numeric(TF_expr[TF_expr$TF == 'PU1',which(names(TF_expr) %in% macro$TISSUE_ABBRV)]),
#                  method='spearman'))
# })
# 
# apply(macro[c('Macrophages','Macrophages_M1','Macrophages_M2','Monocytes')], 1, sum) -> macro$sum_macros
# ggplot(macro, aes(sum_macros,aFC)) +
#   geom_hline(yintercept = 0) + 
#   geom_point(aes(color=TISSUE_ABBRV,fill=TISSUE_ABBRV,shape=TISSUE_ABBRV),size=2) +
#   theme_classic() +
#   scale_color_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
#   scale_fill_manual(values=as.character(macro[order(macro$TISSUE_ABBRV),]$TISSUE_RCOL)) +
#   scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7)))
# cor.test(macro$sum_macros,macro$aFC,method='spearman')
# cor.test(macro[order(macro$TISSUE_ABBRV),]$sum_macros,
#          as.numeric(TF_expr[TF_expr$TF == 'PU1',which(names(TF_expr) %in% macro$TISSUE_ABBRV)]),
#          method='spearman')




tiss_lms = do.call('rbind', apply(tiss_trans[-c(7,24,25,31,35),], 1, function(row) {
  tiss_short = as.character(row[3])
  tiss_long = as.character(row[2])
  print(tiss_long)
  
  expr_file=paste0("expr/by_tiss/",tiss_long,".v8.deseq_log2_expression.ENSG00000107798.17.bed")
  cov_file=paste0("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/",tiss_long,".v8.covariates.txt")
  
  tiss_expr=read.table(expr_file,
                       header=TRUE, comment.char = '')
  tiss_expr_long=melt(tiss_expr,id.vars = c('X.chr','start','end','gene_id'))
  
  tiss_covar=as.data.frame(t(read.table(cov_file, header=TRUE))[-1,])
  tiss_covar = as.data.frame(apply(tiss_covar, c(1,2), function(x) {as.numeric(x)}))
  names(tiss_covar) <- unlist(read.table(cov_file,header=TRUE)[1])
  
  tiss_data = merge(GTs,tiss_expr_long,by.x='sample',by.y='variable')
  tiss_data = tiss_data[ ! (tiss_data$GTs %in% c("./._0/0","1/1_./.") ),]
  
  tiss_data = merge(tiss_data, tiss_covar,
                    by.x='sample', by.y='row.names')
  
  lin_mod=summary(lm( as.formula( paste0( "value ~ PU.1_SNP_DS + GWAS_SNP_DS +",
                                          paste(names(tiss_data)[which(names(tiss_data) == 'PC1'):which(names(tiss_data)=='sex')],collapse='+')) ),
                      data=tiss_data))
  
  data.frame(tiss=tiss_short, tiss_long=tiss_long, tiss_col=row[4],
             PU.1_B = lin_mod$coefficients[2,1], GWAS_B = lin_mod$coefficients[3,1],
             PU.1_p = lin_mod$coefficients[2,4], GWAS_p = lin_mod$coefficients[3,4])
  
}) )

tiss_lms %>% filter(PU.1_p < 0.05) %>% pull(tiss_long)
tiss_lms %>% filter(GWAS_p < 0.05) %>% pull(tiss_long)

tiss_lms = tiss_lms %>%
  mutate(sig_info=paste(PU.1_p<0.05,GWAS_p<0.05,sep='_')) %>%
  mutate(either_sig = PU.1_p <0.05 | GWAS_p < 0.05) %>%
  mutate(both_sig = PU.1_p < 0.05 & GWAS_p < 0.05) %>%
  mutate(which_sig = factor(ifelse(PU.1_p < 0.05 & GWAS_p < 0.05,'both',
                            ifelse(PU.1_p < 0.05,'PU.1',
                                   ifelse(GWAS_p < 0.05,'GWAS','neither'))),
         levels=c('both','GWAS','PU.1','neither')))

tiss_lms %>%
  ggplot(aes(PU.1_B, GWAS_B)) +
  geom_point(aes(col=which_sig,fill=which_sig,shape=which_sig)) +
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  scale_color_manual(values=c('blue','green3','darkorchid1','gray')) +
  scale_fill_manual(values=c('blue','green3','darkorchid1','gray')) +
  scale_shape_manual(values=c(22,24,25,1))

tiss_lms %>%
  ggplot(aes(PU.1_B, GWAS_B)) +
  geom_point(aes(col=as.character(tiss), shape=either_sig)) +
  scale_color_manual(values=as.character(tiss_lms$tiss_col)) +
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)

tiss_lms %>%
  ggplot(aes(PU.1_B, GWAS_B)) +
  geom_point(aes(col=as.character(tiss), shape=both_sig)) +
  scale_color_manual(values=as.character(tiss_lms$tiss_col)) +
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)

tiss_lms %>%
  ggplot(aes(PU.1_B, GWAS_B)) +
  geom_point(aes(col=tiss, fill=tiss, shape=which_sig)) +
  scale_color_manual(values=as.character(tiss_lms$tiss_col)) +
  scale_fill_manual(values=as.character(tiss_lms$tiss_col)) +
  scale_shape_manual(values=c(22,24,25,1)) +
  theme_classic() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)\


"expr/"

# temp$sig_info <- paste(temp$p_PU.1<0.05,temp$p_GWAS<0.05,sep='_')
# temp$int_sig <- as.character(temp$p_PU.1_GWAS < 0.05)
# ggplot(temp,aes(b_PU.1,b_GWAS)) + 
#   geom_point(aes(col=sig_info,shape=int_sig)) + 
#   theme_classic() +
#   geom_hline(yintercept=0) +
#   geom_vline(xintercept=0)
# 
# temp_cols = merge(temp,tiss_col,by.x='tiss_long',by.y='TISSUE_NAME')
# ggplot(temp_cols,aes(b_PU.1,b_GWAS)) + 
#   geom_point(aes(col=as.character(TISSUE_ABBRV),shape=as.character(TISSUE_ABBRV),fill=as.character(TISSUE_ABBRV))) + 
#   scale_color_manual(values=as.character(temp_cols[order(temp_cols$TISSUE_ABBRV),]$TISSUE_RCOL)) +
#   scale_fill_manual(values=as.character(temp_cols[order(temp_cols$TISSUE_ABBRV),]$TISSUE_RCOL)) +
#   scale_shape_manual(values=c(rep(c(15,16,17,18,19,23,25),7))) +
#   theme_classic() +
#   geom_hline(yintercept=0) +
#   geom_vline(xintercept=0)
            