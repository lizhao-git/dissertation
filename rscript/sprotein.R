# LOAD LIBRARIES ####

## DATAFORMAT TRANSFORM ####
library(dplyr)
library(tidyr)
library(stringr)

## DATA PROCESSING ####
library(MatchIt)

## PLOT LIBRARIES ####
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(patchwork)

# LOAD DATA ####

## LOAD SPEPTIDE MATRIX ####
speptide <- read.csv("data/sprotein/speptide.csv")

## LOAD GENEINFO MATRIX ####
geneinfo <- read.csv("data/geneinfo/LncBookv2.0_GENCODEv34_GRCh38.filtered.gene.info.tsv", sep = "\t") %>% dplyr::select(c(geneid, symbol, gene_len))


# STAT ####
length(unique(speptide$smprotid)) # TOTAL 34012 speptides
length(unique((speptide %>% dplyr::filter(evidence == 'Ribo-seq'))$smprotid)) # TOTAL 33679 speptides from Ribo-seq
length(unique((speptide %>% dplyr::filter(evidence == 'Mass spectrometry'))$smprotid)) # TOTAL 333 speptides from Mass spectrometry
length(unique(speptide$geneid)) # 5743 lncRNAs with at least 1 speptide
length(unique(speptide$transcriptid)) # 5743 lncRNAs with at least 1 speptide

# PLOT ####
speptide.count.in.lncrna.gene <- speptide %>% 
  dplyr::group_by(geneid) %>%
  dplyr::summarize(smprotid_count = n_distinct(smprotid)) %>% 
  dplyr::count(smprotid_count)

speptide.count.large <- speptide.count.in.lncrna.gene %>% 
  dplyr::filter(smprotid_count > 10) %>%
  dplyr::summarize(smprotid_count = ">10", n = sum(n))

speptide.count.middle <- speptide.count.in.lncrna.gene %>% 
  dplyr::filter(smprotid_count > 3, smprotid_count <= 10) %>%
  dplyr::summarize(smprotid_count = "3-10", n = sum(n))


rbind(speptide.count.in.lncrna.gene %>% dplyr::filter(smprotid_count <= 3), speptide.count.middle, speptide.count.large) %>% 
  dplyr::mutate(ratio = n / sum(n)) %>% 
  ggplot(aes(x = "", y = ratio, fill = as.factor(smprotid_count))) + 
  geom_bar(stat = "identity") + 
  geom_text(
    aes(label = paste(smprotid_count, "\n", "(", round(ratio * 100, 1), "%)")), 
    position = position_stack(vjust = 0.5)
  ) + 
  scale_fill_manual(values = c('1' = '#')) + 
  coord_polar("y", start = 0) + 
  labs(title = "") + 
  theme_void() + 
  theme(legend.position = 'none')


speptide %>% 
  dplyr::group_by(geneid) %>%
  dplyr::summarize(smprotid_count = n_distinct(smprotid)) %>% 
  dplyr::inner_join(geneinfo, by = 'geneid') %>% 
  dplyr::mutate(count = log10(smprotid_count), gene_len = log10(gene_len/1000)) %>% 
  ggplot(aes(x = gene_len, y = count)) + 
  geom_point() + 
  labs(x = 'Log10(gene length)', y = 'Log10(#small protein site)') + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )
