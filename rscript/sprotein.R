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
geneinfo <- read.delim("data/geneinfo/gene_info.csv", sep = ",") %>% 
  dplyr::select(c(geneid, symbol, length)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(symbol = if_else(symbol == 'N/A', geneid, symbol))

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


speptide.ration.p <- rbind(speptide.count.in.lncrna.gene %>% dplyr::filter(smprotid_count <= 3), speptide.count.middle, speptide.count.large) %>% 
  dplyr::mutate(ratio = n / sum(n)) %>% 
  ggplot(aes(x = "", y = ratio, fill = as.factor(smprotid_count))) + 
  geom_bar(stat = "identity") + 
  geom_text(
    aes(label = paste(smprotid_count, "\n", "(", round(ratio * 100, 1), "%)")), 
    position = position_stack(vjust = 0.5)
  ) + 
  scale_fill_manual(values = c('>10' = '#6D9578', '3-10' = '#1A9E74', '3' = '#B8DCB8', '2' = '#BCE2BF', '1' = '#DFF2E0')) + 
  coord_polar("y", start = 0) + 
  labs(title = "") + 
  theme_publish() + 
  theme_void() + 
  theme(legend.position = 'none')


speptide.number.vs.length.p <- speptide %>% 
  dplyr::group_by(geneid) %>%
  dplyr::summarize(smprotid_count = n_distinct(smprotid)) %>% 
  dplyr::inner_join(geneinfo, by = 'geneid') %>% 
  dplyr::mutate(count = log10(smprotid_count), gene_len = log10(length/1000)) %>% 
  ggplot(aes(x = gene_len, y = count)) + 
  geom_point(colour = '#B8DCB8') + 
  labs(x = 'Log10(gene length)', y = 'Log10(#small protein site)') + 
  theme_publish() + 
  theme(
    aspect.ratio = 1
  )

speptide.rank.p <- speptide %>% 
  dplyr::group_by(geneid) %>%
  dplyr::summarize(smprotid_count = n_distinct(smprotid)) %>% 
  dplyr::inner_join(geneinfo, by = 'geneid') %>% 
  dplyr::mutate(count = log10(smprotid_count), gene_len = length/1000) %>% 
  dplyr::mutate(average = count/gene_len) %>% 
  dplyr::arrange(desc(average)) %>% 
  dplyr::mutate(rank = row_number()) %>% 
  dplyr::mutate(symbol = if_else(row_number() > 10, NA_character_, symbol)) %>% 
  ggplot(aes(x = rank, y = average, label = symbol)) + 
  geom_point(colour = '#B8DCB8') +
  geom_text_repel(nudge_y = 2, color = '#B8DCB8') +
  labs(x = NULL, y = "Average(#small protein site)") +
  theme_publish()

patched.p <- (speptide.ration.p | speptide.number.vs.length.p | speptide.rank.p) + 
  plot_layout(nrow = 4, ncol = 3) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/speptide.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

