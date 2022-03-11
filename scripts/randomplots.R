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




# class Residualizer(object):
# def __init__(self, C_t):
#   # center and orthogonalize
#   self.Q_t, _ = torch.qr(C_t - C_t.mean(0))
# self.dof = C_t.shape[0] - 2 - C_t.shape[1]
# 
# def transform(self, M_t, center=True):
#   """Residualize rows of M wrt columns of C"""
# if center:
#   M0_t = M_t - M_t.mean(1, keepdim=True)
# else:
#   M0_t = M_t
# return M_t - torch.mm(torch.mm(M0_t, self.Q_t), self.Q_t.t())  # keep original mean
# # 
# 
# 
# def calculate_interaction_nominal(genotypes_t, phenotypes_t, interaction_t, residualizer,
#                                   return_sparse=False, tstat_threshold=None):
#   """
# genotypes_t:   [num_genotypes x num_samples]
# phenotypes_t:   [num_phenotypes x num_samples]
# interaction_t: [1 x num_samples]
# """
# ng, ns = genotypes_t.shape
# nps = phenotypes_t.shape[0]
# 
# # centered inputs
# g0_t = genotypes_t - genotypes_t.mean(1, keepdim=True)
# gi_t = genotypes_t * interaction_t
# gi0_t = gi_t - gi_t.mean(1, keepdim=True)
# i0_t = interaction_t - interaction_t.mean()
# p0_t = phenotypes_t - phenotypes_t.mean(1, keepdim=True)
# 
# # residualize rows
# g0_t = residualizer.transform(g0_t, center=False)
# gi0_t = residualizer.transform(gi0_t, center=False)
# p0_t = residualizer.transform(p0_t, center=False)
# i0_t = residualizer.transform(i0_t, center=False)
# i0_t = i0_t.repeat(ng, 1)
# 
# # regression (in float; loss of precision may occur in edge cases)
# X_t = torch.stack([g0_t, i0_t, gi0_t], 2)  # ng x ns x 3
# Xinv = torch.matmul(torch.transpose(X_t, 1, 2), X_t).inverse() # ng x 3 x 3
# # Xinv = tf.linalg.inv(tf.matmul(X_t, X_t, transpose_a=True))  # ng x 3 x 3
# # p0_tile_t = tf.tile(tf.expand_dims(p0_t, 0), [ng,1,1])  # ng x np x ns
# p0_tile_t = p0_t.unsqueeze(0).expand([ng, *p0_t.shape])  # ng x np x ns
# 
# # calculate b, b_se
# # [(ng x 3 x 3) x (ng x 3 x ns)] x (ng x ns x np) = (ng x 3 x np)
# b_t = torch.matmul(torch.matmul(Xinv, torch.transpose(X_t, 1, 2)), torch.transpose(p0_tile_t, 1, 2))
# dof = residualizer.dof - 2
# if nps==1:
#   r_t = torch.matmul(X_t, b_t).squeeze() - p0_t
# rss_t = (r_t*r_t).sum(1)
# b_se_t = torch.sqrt(Xinv[:, torch.eye(3, dtype=torch.uint8).bool()] * rss_t.unsqueeze(1) / dof)
# b_t = b_t.squeeze(2)
# # r_t = tf.squeeze(tf.matmul(X_t, b_t)) - p0_t  # (ng x ns x 3) x (ng x 3 x 1)
# # rss_t = tf.reduce_sum(tf.multiply(r_t, r_t), axis=1)
# # b_se_t = tf.sqrt( tf.matrix_diag_part(Xinv) * tf.expand_dims(rss_t, 1) / dof )
# else:
#   # b_t = tf.matmul(p0_tile_t, tf.matmul(Xinv, X_t, transpose_b=True), transpose_b=True)
#   # convert to ng x np x 3??
#   r_t = torch.matmul(X_t, b_t) - torch.transpose(p0_tile_t, 1, 2)  # (ng x ns x np)
# rss_t = (r_t*r_t).sum(1)  # ng x np
# b_se_t = torch.sqrt(Xinv[:, torch.eye(3, dtype=torch.uint8).bool()].unsqueeze(-1).repeat([1,1,nps]) * rss_t.unsqueeze(1).repeat([1,3,1]) / dof)
# # b_se_t = tf.sqrt(tf.tile(tf.expand_dims(tf.matrix_diag_part(Xinv), 2), [1,1,nps]) * tf.tile(tf.expand_dims(rss_t, 1), [1,3,1]) / dof) # (ng x 3) -> (ng x 3 x np)
# 
# tstat_t = (b_t.double() / b_se_t.double()).float()  # (ng x 3 x np)



covars_only = WHLBLD_data[c(10:77)]
apply(covars_only, 2, mean)
covars_centered = as.data.frame(apply(covars_only, 2, function(col) {col - mean(col)}))
Q_t = qr(covars_centered)$qr

G_t = WHLBLD_data$GWAS_SNP_DS
G_t_c = G_t - mean(G_t)
G_t_q = G_t_c - (G_t_c %*% Q_t) %*% t(Q_t)

P_t = WHLBLD_data$value
P_t_c = P_t - mean(P_t)
P_t_q = P_t_c - (P_t_c %*% Q_t) %*% t(Q_t)

I_t = WHLBLD_data$PU.1_SNP_DS
I_t_c = I_t - mean(I_t)
I_t_q = I_t_c - (I_t_c %*% Q_t) %*% t(Q_t)

GI_t = G_t * I_t
GI_t_c = GI_t - mean(GI_t)
GI_t_q = GI_t_c - (GI_t_c %*% Q_t) %*% t(Q_t)

X_t = cbind(t(G_t_q), t(I_t_q), t(GI_t_q))
Xinv = ginv(t(X_t) %*% X_t)

b_t = (Xinv %*% t(X_t)) %*% t(P_t_q)
r_t = (X_t %*% b_t) - t(P_t_q)
rss_t = sum(r_t*r_t)
# b_se_t = sqrt(  )

lm(t(P_t_q) ~ t(G_t_q) + t(I_t_q) + t(GI_t_q))
lm(P_t ~ G_t + I_t + GI_t)
lm(P_t_c ~ G_t_c + I_t_c + GI_t_c)
lm(t(P_t_q) ~ G_t_c + I_t_c + GI_t_c)







gt = c(rep(0,100),rep(1,67),rep(2,33))
values = sample(c(rep(-1,50),rep(-.5,25),rep(0,50),rep(.5,25),rep(1,50)),size=200)
tf = values + rnorm(length(gt))/5
cov = values + rnorm(length(gt))/5
pheno = gt*2 + values*3 + rnorm(length(gt))

summary(lm(pheno ~ gt))
summary(lm(pheno ~ tf))
summary(lm(pheno ~ cov))
cor(cov, tf)
summary(lm(pheno ~ values))

summary(lm(pheno ~ gt + tf + cov))


pheno_rescov = pheno - summary(lm(pheno ~ cov))$coeff[2,1] * cov
gt_rescov = gt - summary(lm(gt ~ cov))$coeff[2,1] * cov
tf_rescov = tf - summary(lm(tf ~ cov))$coeff[2,1] * cov
summary(lm(pheno_rescov ~ gt + tf ))
summary(lm(pheno_rescov ~ gt_rescov + tf_rescov ))



read.table("~/lap_lab_folder/data/gtex/v8/expression_bytiss/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
           header=TRUE,sep='\t')














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

WHLBLD_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt",
                                        header=TRUE))[-1,])
WHLBLD_covar = as.data.frame(apply(WHLBLD_covar, c(1,2), function(x) {as.numeric(x)}))
names(WHLBLD_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt",
                                         header=TRUE)[1])

WHLBLD_data = merge(GTs,WHLBLD_expr_long,by.x='sample',by.y='variable')
WHLBLD_data = WHLBLD_data[ ! (WHLBLD_data$GTs %in% c("./._0/0","1/1_./.") ),]

#WHLBLD_data = merge(WHLBLD_data, WHLBLD_covar,
#                    by.x='sample', by.y='row.names')

table(WHLBLD_data$PU.1_SNP_GT,WHLBLD_data$GWAS_SNP_GT)

boxplot(log2(WHLBLD_data$value) ~ WHLBLD_data$PU.1_SNP_GT)






ADPSBQ_expr=read.table("expr/by_tiss/Adipose_Subcutaneous.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
ADPSBQ_expr_long=melt(ADPSBQ_expr,id.vars = c('X.chr','start','end','gene_id'))

ADPSBQ_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt",
                                        header=TRUE))[-1,])
ADPSBQ_covar = as.data.frame(apply(ADPSBQ_covar, c(1,2), function(x) {as.numeric(x)}))
names(ADPSBQ_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt",
                                         header=TRUE)[1])

ADPSBQ_data = merge(GTs,ADPSBQ_expr_long,by.x='sample',by.y='variable')
ADPSBQ_data = ADPSBQ_data[ ! (ADPSBQ_data$GTs %in% c("./._0/0","1/1_./.") ),]

#ADPSBQ_data = merge(ADPSBQ_data, ADPSBQ_covar,
#                    by.x='sample', by.y='row.names')

boxplot(log2(ADPSBQ_data$value) ~ ADPSBQ_data$PU.1_SNP_GT)




SKINNS_expr=read.table("expr/by_tiss/Skin_Not_Sun_Exposed_Suprapubic.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
SKINNS_expr_long=melt(SKINNS_expr,id.vars = c('X.chr','start','end','gene_id'))

SKINNS_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Skin_Not_Sun_Exposed_Suprapubic.v8.covariates.txt",
                                        header=TRUE))[-1,])
SKINNS_covar = as.data.frame(apply(SKINNS_covar, c(1,2), function(x) {as.numeric(x)}))
names(SKINNS_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Skin_Not_Sun_Exposed_Suprapubic.v8.covariates.txt",
                                         header=TRUE)[1])

SKINNS_data = merge(GTs,SKINNS_expr_long,by.x='sample',by.y='variable')
SKINNS_data = SKINNS_data[ ! (SKINNS_data$GTs %in% c("./._0/0","1/1_./.") ),]

#SKINNS_data = merge(SKINNS_data, SKINNS_covar,
#                    by.x='sample', by.y='row.names')

table(SKINNS_data$PU.1_SNP_GT,SKINNS_data$GWAS_SNP_GT)

boxplot(log2(SKINNS_data$value) ~ SKINNS_data$PU.1_SNP_GT)







FIBRBLS_expr=read.table("expr/by_tiss/Cells_Cultured_fibroblasts.v8.deseq_log2_expression.ENSG00000107798.17.bed",header=TRUE, comment.char = '')
FIBRBLS_expr_long=melt(FIBRBLS_expr,id.vars = c('X.chr','start','end','gene_id'))

FIBRBLS_covar=as.data.frame(t(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_Cultured_fibroblasts.v8.covariates.txt",
                                        header=TRUE))[-1,])
FIBRBLS_covar = as.data.frame(apply(FIBRBLS_covar, c(1,2), function(x) {as.numeric(x)}))
names(FIBRBLS_covar) <- unlist(read.table("~/data/GTEx_Analysis_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Cells_Cultured_fibroblasts.v8.covariates.txt",
                                         header=TRUE)[1])

FIBRBLS_data = merge(GTs,FIBRBLS_expr_long,by.x='sample',by.y='variable')
FIBRBLS_data = FIBRBLS_data[ ! (FIBRBLS_data$GTs %in% c("./._0/0","1/1_./.") ),]

#FIBRBLS_data = merge(FIBRBLS_data, FIBRBLS_covar,
#                    by.x='sample', by.y='row.names')

table(FIBRBLS_data$PU.1_SNP_GT,FIBRBLS_data$GWAS_SNP_GT)

boxplot(log2(FIBRBLS_data$value) ~ FIBRBLS_data$PU.1_SNP_GT)


library(dplyr)

fourtiss_data = rbind(WHLBLD_data %>% mutate(tiss='WHLBLD'),
      ADPSBQ_data %>% mutate(tiss='ADPSBQ')) %>%
  rbind(SKINNS_data %>% mutate(tiss='SKINNS')) %>%
  rbind(FIBRBLS_data %>% mutate(tiss='FIBRBLS'))

fourtiss_data %>%
  ggplot(aes(GWAS_SNP_GT, value)) +
  geom_boxplot(aes(col=tiss)) +
  facet_wrap(~tiss,
             scales = 'free_y') +
  theme_classic() +
  scale_color_manual(values=c('tan1','lightblue3','royalblue3','magenta')) +
  xlab("SNP genotype") +
  ylab("log2(LIPA expression)")

fourtiss_lms = fourtiss_data %>%
  group_by(tiss) %>%
  mutate(gt_val = as.numeric(GWAS_SNP_GT) - 1) %>%
  summarize(lm_int = summary(lm(value ~ gt_val))$coefficients[1,1],
            lm_slope = summary(lm(value ~ gt_val))$coefficients[2,1])

fourtiss_data %>%
  ggplot(aes(GWAS_SNP_GT, value)) +
  geom_boxplot(aes(col=tiss)) +
  facet_wrap(~tiss,
             scales = 'free_y') +
  theme_classic() +
  scale_color_manual(values=c('tan1','lightblue3','royalblue3','magenta')) +
  xlab("SNP genotype") +
  ylab("log2(LIPA expression)") +
  geom_abline(data=fourtiss_lms, 
              aes(slope=lm_slope, intercept=lm_int))



fourtiss_data %>%
  group_by(tiss) %>%
  mutate(gt_val = as.numeric(GWAS_SNP_GT) - 1) %>%
  summarize(lm_int = summary(lm(2^value ~ gt_val))$coefficients[1,1],
            lm_slope = summary(lm(2^value ~ gt_val))$coefficients[2,1],
            med_expr = median(value)) %>%
  mutate(afc = (2*lm_slope/lm_int + 1),
         log2afc = log2(afc)) %>%
  ggplot(aes(med_expr, log2afc)) +
  geom_point(aes(col=tiss),
             size=5) +
  scale_color_manual(values=c('tan1','lightblue3','royalblue3','magenta')) +
  theme_classic() +
  xlab("log2( median gene expression )") +
  ylab("aFC")
