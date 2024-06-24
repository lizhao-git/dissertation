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

## LOAD DIFF MATRIX ####
de.hiv <- read.csv("data/featured/virus/de_hiv_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'HIV') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.hbv <- read.csv("data/featured/virus/de_hbv_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'HBV') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.hcv <- read.csv("data/featured/virus/de_hcv_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'HCV') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.covid <- read.csv("data/featured/virus/de_sars2_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'COVID') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.combined <- rbind(de.hiv, de.hbv, de.hcv, de.covid)

de.combined %>% 
  dplyr::group_by(disease, status) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(count = if_else(status =='down', -count, count)) %>% 
  ggplot(aes(x= reorder(disease, ), y = count)) + 
  geom_bar(
    aes(fill = status),
    stat = 'identity',
    width = 0.6
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.6, size = 3) +
  scale_fill_manual(values = c("up" = "#C84E6F", "down" = "#43649E")) + 
  #scale_x_discrete(labels = c('Coronary', 'ESCC', 'Pancreatic', 'Hepatocellular', 'Colorectal')) + 
  scale_y_continuous(limits = c(-1000, 500)) + 
  labs(x = NULL, y = '#LncRNA gene') + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

de.combined %>% 
  dplyr::count(geneid) %>% 
  dplyr::count(n)
