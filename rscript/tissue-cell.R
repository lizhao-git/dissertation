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

## LOAD METHYLATION MATRIX ####
tau.hpa <- read.csv("data/featured/hpa/tau_hpa.csv") %>% 
  dplyr::filter(max_value >= 10, (tau <= 0.35) | (tau >= 0.9)) %>% 
  dplyr::mutate(status = if_else(tau >= 0.9, 'specific', 'consistent'))
  
tau.encode <- read.csv("data/featured/encode/tau_encode.csv") %>% 
  dplyr::filter(max_value >= 10, (tau <= 0.35) | (tau >= 0.9)) %>% 
  dplyr::mutate(status = if_else(tau >= 0.9, 'specific', 'consistent'))
  
tau.ccle <- read.csv("data/featured/ccle/tau_ccle.csv") %>% 
  dplyr::filter(max_value >= 10, (tau <= 0.35) | (tau >= 0.9)) %>% 
  dplyr::mutate(status = if_else(tau >= 0.9, 'specific', 'consistent'))

# PLOT ####
tau.hpa.p<- tau.hpa %>% 
  dplyr::group_by(max_name) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(-count) %>% 
  dplyr::mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank, y = count)) + 
  geom_point(
    color = '#B8DCB8',
    stat = 'identity',
    width = 0.8
  ) + 
  geom_label_repel(
    aes(x = rank, y = count, label = max_name),
    direction = 'y',
    size = 3
  ) + 
  labs(x = NULL, y = '#LncRNA') + 
  theme_publish() + 
  theme(
    aspect.ratio = 1
  )

tau.encode.p <- tau.encode %>% 
  dplyr::group_by(max_name) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(-count) %>% 
  dplyr::mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank, y = count)) + 
  geom_point(
    color = '#B8DCB8',
    stat = 'identity',
    width = 0.8
  ) + 
  geom_label_repel(
    aes(x = rank, y = count, label = max_name),
    direction = 'y',
    size = 3
  ) + 
  labs(x = NULL, y = '#LncRNA') + 
  theme_publish() + 
  theme(
    aspect.ratio = 1
  )

tau.ccle.p <- tau.ccle %>% 
  dplyr::group_by(max_name) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(-count) %>% 
  dplyr::mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank, y = count)) + 
  geom_point(
    color = '#B8DCB8',
    stat = 'identity',
    width = 0.8
  ) + 
  geom_label_repel(
    aes(x = rank, y = count, label = max_name),
    direction = 'y',
    size = 3
  ) + 
  labs(x = NULL, y = '#LncRNA') + 
  theme_publish() + 
  theme(
    aspect.ratio = 1
  )

patched.p <- (tau.hpa.p | tau.encode.p | tau.ccle.p) + 
  plot_layout(nrow = 3, ncol = 3) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/tissue.cell.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()
