#!/usr/bin/Rscript
##
##
##  GWASRegression.R
##


library(ggplot2)
library(dplyr)
library(tidyr)


set.seed(5)
snp_a = sample(c(0,1), 500, replace=TRUE)
snp_b = sample(c(0,1), 500, replace=TRUE)
part1=snp_b*sqrt(.5)
part2=sample(c(0,1), 500, replace=TRUE)*sqrt(.5)
snp_c = part1 + part2

cor(snp_a, snp_b)
cor(snp_a, snp_c)
cor(snp_b, snp_c)

set.seed(6)
trait = snp_a*.2 + sample(seq(0,1,.001), 500, replace=TRUE)
summary(lm(trait ~ snp_a))
summary(lm(trait ~ snp_b))


set.seed(6)
trait = snp_a*.2 + snp_b*.2 + sample(seq(0,1,.001), 500, replace=TRUE)
summary(lm(trait ~ snp_a))
summary(lm(trait ~ snp_b))


trait = snp_b*.2 + sample(seq(0,1,.001), 500, replace=TRUE)
summary(lm(trait ~ snp_b))
summary(lm(trait ~ snp_c))
summary(lm(trait ~ part1 + part2 ))

trait = snp_b*.2 + snp_c*.2 + sample(seq(0,1,.001), 500, replace=TRUE)
summary(lm(trait ~ snp_b))
summary(lm(trait ~ snp_c))
summary(lm(trait ~ part1 + part2 ))

lms_both = do.call('rbind', lapply(1:500, function(i) {
  r2 = sample(seq(.1,.9,.01),1)
  snp_b = sample(c(0,1), 5000, replace=TRUE)
  part1=snp_b
  part2=sample(c(0,1), 5000, replace=TRUE)
  parts = cbind(part1,part2)
  snp_c = unlist(sapply(1:5000, function(i) {
    parts[i, sample(c(1,2) , 1, prob=c(sqrt(r2),1-sqrt(r2)))]
  }))
  r_snps = cor(snp_b, snp_c)
  r_snps_sq = r_snps^2
  
  trait_both = snp_b*.1 + snp_c*.1 + sample(seq(0,1,.001), 5000, replace=TRUE)
  
  lm_both_b = summary(lm(trait_both ~ snp_b))
  lm_both_c = summary(lm(trait_both ~ snp_c))
  lm_both_csplit = summary(lm(trait_both ~ part1 + part2))
  
  trait_b = snp_b*.2 + sample(seq(0,1,.001), 5000, replace=TRUE)
  
  lm_b_b = summary(lm(trait_b ~ snp_b))
  lm_b_c = summary(lm(trait_b ~ snp_c))
  lm_b_csplit = summary(lm(trait_b ~ part1 + part2))
  
  data.frame( rep = i, run=c('both','bonly'), r2, r_snps, r_snps_sq,
              lm_b_B = c(lm_both_b$coefficients[2,1],lm_b_b$coefficients[2,1]), 
              lm_b_Bse = c(lm_both_b$coefficients[2,2],lm_b_b$coefficients[2,2]), 
              lm_b_p = c(lm_both_b$coefficients[2,4],lm_b_b$coefficients[2,4]),
              lm_c_B = c(lm_both_c$coefficients[2,1],lm_b_c$coefficients[2,1]), 
              lm_c_Bse = c(lm_both_c$coefficients[2,2],lm_b_c$coefficients[2,2]), 
              lm_c_p = c(lm_both_c$coefficients[2,4],lm_b_c$coefficients[2,4]),
              lm_c_split1_B = c(lm_both_csplit$coefficients[2,1],lm_b_csplit$coefficients[2,1]), 
              lm_c_split1_Bse = c(lm_both_csplit$coefficients[2,2],lm_b_csplit$coefficients[2,2]), 
              lm_c_split1_p = c(lm_both_csplit$coefficients[2,4],lm_b_csplit$coefficients[2,4]), 
              lm_c_split2_B = c(lm_both_csplit$coefficients[3,1],lm_b_csplit$coefficients[3,1]), 
              lm_c_split2_Bse = c(lm_both_csplit$coefficients[3,2],lm_b_csplit$coefficients[3,2]), 
              lm_c_split2_p = c(lm_both_csplit$coefficients[3,4],lm_b_csplit$coefficients[3,4]) )
  
}))

# lms_b = do.call('rbind', lapply(1:100, function(i) {
#   r2 = sample(seq(.1,.9,.01),1)
#   snp_b = sample(c(0,1), 5000, replace=TRUE)
#   part1=snp_b*r2
#   part2=sample(c(0,1), 5000, replace=TRUE)*(1-r2)
#   snp_c = unlist(sapply(1:5000, function(i) {
#     parts[i, sample(c(1,2) , 1, prob=c(sqrt(r2),1-sqrt(r2)))]
#   }))
#   r_snps = cor(snp_b, snp_c)
#   r_snps_sq = r_snps^2
#   
#   trait_b = snp_b*.1 + sample(seq(0,1,.001), 5000, replace=TRUE)
#   
#   lm_b = summary(lm(trait_b ~ snp_b))
#   lm_c = summary(lm(trait_b ~ snp_c))
#   lm_c_split = summary(lm(trait_b ~ part1 + part2))
#   
#   data.frame( run = i, r2, r_snps, r_snps_sq,
#               lm_b_B = lm_b$coefficients[2,1], lm_b_p = lm_b$coefficients[2,4], lm_b_r2 = lm_b$r.squared,
#               lm_c_B = lm_c$coefficients[2,1], lm_c_p = lm_c$coefficients[2,4], lm_c_r2 = lm_c$r.squared,
#               lm_c_split1_B = lm_c_split$coefficients[2,1], lm_c_split1_p = lm_c_split$coefficients[2,4], 
#               lm_c_split2_B = lm_c_split$coefficients[3,1], lm_c_split2_p = lm_c_split$coefficients[3,4], 
#               lm_c_split_r2 = lm_c_split$r.squared )
#   
# }))

## this is the effect I want to be able to find, only using the r_sq and the summ stats from first two
ggplot(lms_both) + geom_point(aes(r_snps_sq, lm_c_split2_B)) + facet_wrap(~run)
ggplot(lms_both) + geom_point(aes(r_snps_sq, -log10(lm_c_split2_p))) + facet_wrap(~run)

temp_lms = lms_both %>%
  pivot_longer(c(lm_b_B,lm_c_B,lm_b_p, lm_c_p)) %>%
  separate(name,into=c(NA,"lm","stat"),sep="_") %>%
  pivot_wider(names_from=stat, values_from=value)
temp_lms %>%
  ggplot() +
  geom_point(aes(r_snps, B)) +
  facet_grid(~ run + lm)
temp_lms %>%
  ggplot() +
  geom_point(aes(r_snps, -log10(p))) +
  facet_grid(~ run + lm)

lms_both %>%
  select(rep,run,r_snps,r_snps_sq,lm_c_B) %>%
  pivot_wider(names_from=run,values_from=lm_c_B) %>%
  ggplot(aes(bonly,both)) + 
  geom_point(aes(color=r_snps_sq)) +
  geom_abline(slope=1, intercept=0)

lms_both %>%
  select(rep,run,r_snps,r_snps_sq,lm_c_p) %>%
  pivot_wider(names_from=run,values_from=lm_c_p) %>%
  ggplot(aes(-log10(bonly),-log10(both))) + 
  geom_point(aes(color=r_snps_sq)) +
  geom_abline(slope=1, intercept=0)

lms_both %>%
  pivot_wider(id_cols=c(rep,r_snps,r_snps_sq),
              names_from=run, values_from=c(lm_b_B,lm_b_p,lm_c_B,lm_c_p)) %>%
  mutate(exp_lm_c_B_bonly = r_snps * lm_b_B_bonly) %>%
  ggplot(aes(exp_lm_c_B_bonly, lm_c_B_bonly)) +
  geom_point() + geom_abline(slope=1, intercept=0)

lms_both %>%
  pivot_wider(id_cols=c(rep,r_snps,r_snps_sq),
              names_from=run, values_from=c(lm_b_B,lm_b_p,lm_c_B,lm_c_p)) %>%
  mutate(exp_lm_c_B_both = r_snps * lm_b_B_both) %>%
  ggplot(aes(exp_lm_c_B_both, lm_c_B_both)) +
  geom_point() + geom_abline(slope=1, intercept=0)

lms_both %>%
  pivot_wider(id_cols=c(rep,r_snps,r_snps_sq),
              names_from=run, values_from=c(lm_b_B,lm_b_p,lm_c_B,lm_c_p)) %>%
  mutate(exp_lm_c_B_both = r_snps * lm_b_B_both,
         diff_lm_c_B_both = lm_c_B_both - exp_lm_c_B_both) %>%
  ggplot(aes(r_snps_sq,diff_lm_c_B_both)) +
  geom_point() + geom_hline(yintercept=0)

lms_both %>%
  pivot_wider(id_cols=c(rep,r_snps,r_snps_sq),
              names_from=run, values_from=c(lm_b_B,lm_b_p,lm_c_B,lm_c_p)) %>%
  mutate(exp_lm_c_B_both = r_snps * lm_b_B_both,
         per_diff_lm_c_B_both = (lm_c_B_both - exp_lm_c_B_both)/lm_c_B_both) %>%
  ggplot(aes(r_snps,per_diff_lm_c_B_both)) +
  geom_point() + geom_hline(yintercept=0) +
  geom_abline(slope=-1,intercept=1)


## Find residual B, p val of c snp when b snp is removed
newvals = lms_both %>%
  mutate(exp_lm_c_B = r_snps * lm_b_B,
         diff_lm_c_B = (lm_c_B - exp_lm_c_B),
         per_diff_lm_c_B = diff_lm_c_B/lm_c_B,
         lm_c_newt = diff_lm_c_B / lm_c_Bse,
         lm_c_newp = pt(-abs(lm_c_newt), 498)*2,
         lm_c_ptest = pt(lm_c_newt, 498*2) )

newvals %>%
  ggplot(aes(lm_c_B, diff_lm_c_B)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=r_snps)) +
  facet_wrap(~run)

newvals %>%
  ggplot(aes(lm_c_B, -log10(lm_c_newp))) +
  geom_hline(yintercept=-log10(0.05)) +
  geom_point(aes(col=r_snps_sq)) +
  facet_wrap(~run)

newvals %>%
  filter(run == 'bonly') %>%
  ggplot() +
  geom_histogram(aes(lm_c_newp))
newvals %>%
  filter(run == 'bonly') %>%
  ggplot() +
  geom_histogram(aes(lm_c_ptest))
sum(newvals %>% filter(run=='bonly') %>% pull(lm_c_newp) < 0.05)


## USE this math on the CAD GWAS
cad_gwas = read.table("GWAS/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.LIPA_region.txt",
                      header=TRUE, sep = '\t')

snp_pu1 = cad_gwas %>% filter(variant_id == 'rs1320496')
snp_gwas = cad_gwas %>% filter(variant_id == 'rs1412445')

snp_gwas_nopu1 = snp_gwas %>%
  mutate(exp_B = sqrt(0.48) * snp_pu1$effect_size,
         B_residual = effect_size - exp_B,
         per_diff_B = B_residual/effect_size,
         newt = B_residual / standard_error,
         newp = pt(-newt, sample_size)*2 )

snp_pu1_nogwas = snp_pu1 %>%
  mutate(exp_B = sqrt(0.45) * snp_gwas$effect_size,
         B_residual = effect_size - exp_B,
         per_diff_B = B_residual/effect_size,
         newt = B_residual / standard_error,
         newp = pt(-newt, sample_size)*2 )

rbind(snp_gwas_nopu1, snp_pu1_nogwas) %>%
  select(variant_id, effect_size, standard_error, zscore, pvalue, B_residual, newt, newp)



