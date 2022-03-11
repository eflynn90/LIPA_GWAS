#!/usr/bin/Rscript
##
##  TF_binding.R
##
##  EDF 2/12/21
##

library(dplyr)
library(seqminer)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

final_tfs = read.table("../TFi-eQTL/input_files/TFs.info.curated.txt",
                       header=TRUE,sep='\t')

chipseq_info = read.table("tf_chipseq_20200110/info_files/final_chipseqs.list",
                          header=FALSE, sep='\t') %>%
  mutate(tabix_file = gsub(pattern = 'processed', replacement = 'tabix',V1))
chipseq_info_2 = read.table("tf_chipseq_20210212/tfs.ct.list",
                            header=FALSE, sep='\t')

chipseq_overlap = do.call('rbind', apply(chipseq_info, 1, function(rowi) {
  tabix_file = paste0("tf_chipseq_20200110/",rowi[5])
  tabix.read.table(tabix_file, "chr10:89243000-89243500") %>%
    mutate(tf=rowi[3], ct=rowi[2])
}) )
chipseq_overlap_2 = do.call('rbind', apply(chipseq_info_2, 1, function(rowi) {
  tabix_file = paste0("tf_chipseq_20210212/merge_files/",rowi[1],'.',rowi[2],".bed.gz")
  tabix.read.table(tabix_file, "chr10:89243000-89243500") %>%
    mutate(tf=rowi[1], ct=rowi[2])
}) ) %>% select(V1,V2,V3,tf,ct)

chipseq_overlap_all = rbind(chipseq_overlap, chipseq_overlap_2)
chipseq_overlap_all %>%
  ggplot() +
  geom_segment(aes(x=tf,xend=tf,y=V2,yend=V3)) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  coord_flip() +
  theme_classic()

motif_files = list.files("motif_overlap/motif_overlap_vars", "*.txt.gz$",
                         full.names=TRUE)

motif_overlap = do.call('rbind', lapply(motif_files, function(filei) {
  print(filei)
  tabix.read.table(filei, "chr10:89243000-89243500") %>%
          mutate(tf=strsplit(filei,"[./]")[[1]][3])
}))

motif_overlap %>%
  filter(V8=='Y') %>%
  ggplot() +
  geom_point(aes(x=tf,y=V3)) +
  theme_classic() +
  coord_flip()

both_overlap = merge(chipseq_overlap_all, motif_overlap %>% filter(V8=='Y'), 
      by='tf', all = TRUE)
both_overlap_gen = merge(chipseq_overlap_all %>% mutate(tf_gen = gsub('[[:digit:]]+.*','',tf)), 
                         motif_overlap %>% filter(V8=='Y') %>% mutate(tf_gen = gsub('[[:digit:]]+.*','',tf)), 
                     by='tf_gen', all = TRUE)

both_overlap %>%
  ggplot() +
  geom_segment(aes(x=tf,xend=tf,y=V2.x,yend=V3.x)) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf,y=V3.y)) +
  coord_flip() +
  theme_classic()

both_overlap %>%
  filter(tf %in% c(as.character(final_tfs$TF),'SPI1','STAT1','MEF2A','MEF2B')) %>%
  ggplot() +
  geom_segment(aes(x=tf,xend=tf,y=V2.x,yend=V3.x)) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf,y=V3.y)) +
  coord_flip() +
  theme_classic()

both_overlap_gen %>%
  ggplot() +
  geom_linerange(aes(x=tf_gen, ymin=V2.x, ymax=V3.x, group=tf.x, col=ct), 
                 position = position_dodge(width=.8), size=2) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf_gen,y=V3.y)) +
  coord_flip() +
  theme_classic()

both_overlap_gen %>%
  filter(!is.na(tf.x), !is.na(tf.y)) %>%
  ggplot() +
  geom_linerange(aes(x=tf_gen, ymin=V2.x, ymax=V3.x, group=tf.x, col=ct), 
                 position = position_dodge(width=.8), size=2) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf_gen,y=V3.y)) +
  coord_flip() +
  theme_classic()

manual_motif_info_2013 = data.frame(
  pos=c(rep(89243047,3),
        rep(89243088,12),
        rep(89243170,4)),
  tf=c('ETS','MAF','MEF2',
       'BCL','ETS','IRF','MAF','NFKB','NKX','PLAG1','SPI1','PAX5','RXRA','STAT','P300',
       'ETS','FEV','PAX5','STAT')
)

both_overlap_manual = merge(chipseq_overlap_all %>% mutate(tf_gen = gsub('[[:digit:]]+.*','',tf)), 
                            manual_motif_info_2013 %>% mutate(tf_gen = gsub('[[:digit:]]+.*','',tf)), 
                            by.x='tf_gen', by.y='tf_gen', all = TRUE) %>%
  mutate(tf_ct = paste(tf.x,ct))

########################################################################################################################
#### THIS IS PLOT TO SHARE!!! ######
both_overlap_manual %>%
  mutate(`ChIPseq Cell Type`=ct) %>%
  ggplot() +
  geom_linerange(aes(x=tf_gen, ymin=V2,ymax=V3, group=tf_ct, col=`ChIPseq Cell Type`), position = position_dodge(width=.8), size=2) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf_gen,y=pos), size=2) +
  coord_flip() +
  theme_classic() +
  xlab('Transcription Factor Family') +
  ylab('Position on chr10')
########################################################################################################################

both_overlap_manual %>% 
  filter(tf_gen %in% c('SPI','MEF','ETS','BCL','STAT')) %>%
  ggplot() +
  geom_linerange(aes(x=tf_gen, ymin=V2,ymax=V3, group=tf_ct, col=ct), position = position_dodge(width=.8), size=2) +
  geom_hline(yintercept=c(89243047,89243088,89243170)) +
  geom_point(aes(x=tf_gen,y=pos), size=2) +
  scale_color_manual(values=c(hue_pal()(4)[1],hue_pal()(4)[2],hue_pal()(4)[3]), na.value='gray') +
  coord_flip() +
  theme_classic()










asb_ADASTRA = read.table("temp.ADASTRA.txt", fill=TRUE, comment.char = '',
                         header=TRUE, sep='\t') %>%
  separate(X.chr, c(NA,NA,NA,NA,NA,NA,'tf'), sep='[/_]') %>%
  mutate(rat_mostsig_ref1 = altc_mostsig_ref/refc_mostsig_ref, 
         rat_mostsig_ref2 = refc_mostsig_ref/altc_mostsig_ref,
         rat_mostsig_ref = ifelse(rat_mostsig_ref1 > rat_mostsig_ref2,
                                  rat_mostsig_ref1, rat_mostsig_ref2)) %>%
  mutate(rat_mostsig_alt1 = altc_mostsig_alt/refc_mostsig_alt, 
         rat_mostsig_alt2 = refc_mostsig_alt/altc_mostsig_alt,
         rat_mostsig_alt = ifelse(rat_mostsig_alt1 > rat_mostsig_alt2,
                                  rat_mostsig_alt1, rat_mostsig_alt2))

asb_ADASTRA %>%
  ggplot(aes(tf)) +
  geom_point(aes(y=rat_mostsig_ref, fill=ID, shape=ID), position = position_nudge(-.1)) +
  geom_point(aes(y=rat_mostsig_alt, fill=ID, shape=ID), position = position_nudge(.1)) +
  scale_fill_manual(values=c('darkorchid3','red')) +
  scale_shape_manual(values=c(23,21)) +
  geom_text(aes(y=BAD_mostsig_ref, label="--"), position = position_nudge(-.1)) +
  geom_text(aes(y=BAD_mostsig_alt, label="--"), position = position_nudge(.1)) +
  theme_classic() +
  ylab("Allellic Ratio") +
  ylim(1,NA)

asb_ADASTRA %>%
  mutate(SNP=ID) %>%
  ggplot(aes(tf)) +
  geom_point(aes(y=log2(rat_mostsig_ref), fill=SNP, shape=SNP), position = position_nudge(-.1)) +
  geom_point(aes(y=log2(rat_mostsig_alt), fill=SNP, shape=SNP), position = position_nudge(.1)) +
  scale_fill_manual(values=c('darkorchid3','red')) +
  scale_shape_manual(values=c(23,21)) +
  geom_text(aes(y=log2(BAD_mostsig_ref), label="--"), position = position_nudge(-.1)) +
  geom_text(aes(y=log2(BAD_mostsig_alt), label="--"), position = position_nudge(.1)) +
  theme_classic() +
  ylab("Log2 Allellic Ratio") +
  ylim(0,NA)



PU1_ENCODE_manual = data.frame(samp=rep('GM12891 LCL',6),
                               rep=c(1,1,1,2,2,2),
                               snp=factor(rep(c('rs1412445','rs1320496','rs1412444'),2),
                                          levels=c('rs1412445','rs1320496','rs1412444')),
                               ref_c=c(3,  4,1,  3,  3,1),
                               alt_c=c(19,  27,1,  18,  11,11))

PU1_ENCODE_manual %>%
  pivot_longer(cols=c(ref_c,alt_c), 
               names_to='count',
               values_to='read_count') %>%
  mutate(allele=ifelse(count=='ref_c','C',
                       ifelse(count=='alt_c','T',''))) %>%
  ggplot(aes(snp,read_count)) +
  geom_col(aes(fill=allele, group=allele),position=position_dodge()) +
  theme_classic() +
  scale_fill_manual(values=c('blue','red')) +
  facet_wrap(~paste(samp,'rep',rep),
             ncol=1) +
  ylab('read count') +
  xlab('SNP') +
  ggtitle("SPI1 allele-specific ChIPseq binding")
  



STAT1_ADASTRA_table = asb_ADASTRA %>%
  filter(tf == 'STAT1') %>%
  pivot_longer(cols=c(refc_mostsig_ref,altc_mostsig_ref),
               names_to='count',
               values_to='read_count') %>%
  mutate(allele = ifelse(count=='refc_mostsig_ref','C',
                       ifelse(count=='altc_mostsig_ref','T','')),
         sample='CD14+ Monocytes',
         SNP=factor(ID,
                    levels=c('rs1412445','rs1320496','rs1412444')))

STAT1_ADASTRA_table %>%
  ggplot(aes(SNP,read_count)) +
  geom_col(aes(fill=allele, group=allele),position=position_dodge()) +
  theme_classic() +
  scale_fill_manual(values=c('blue','red')) +
  ylab('read count') +
  xlab('SNP')  +
  ggtitle("STAT1 allele-specific ChIPseq binding") +
  facet_wrap(~ sample, ncol=1)




########### MAIN FIG ################
## combine reps for SPI1
## show nothing for STAT1?

PU1_ENCODE_manual %>%
  group_by(snp) %>%
  summarize(ref_c = sum(ref_c),
            alt_c = sum(alt_c)) %>%
  pivot_longer(cols=c(ref_c,alt_c), 
               names_to='count',
               values_to='read_count') %>%
  mutate(allele=ifelse(count=='ref_c','C',
                       ifelse(count=='alt_c','T',''))) %>%
  filter(snp=='rs1320496') %>%
  ggplot(aes(snp,read_count)) +
  geom_col(aes(fill=allele, group=allele),position=position_dodge()) +
  theme_classic() +
  scale_fill_manual(values=c('blue','red')) +
  ylab('read count') +
  xlab('') +
  ggtitle("SPI1 allele-specific ChIPseq binding")














