# LOAD LIBRARIES ####

## ENVIRONMENT SETTINGS ####
library(here)

## DATAFORMAT TRANSFORM ####
library(dplyr)
library(tidyr)
library(stringr)

## DATA PROCESSING ####

## PLOT LIBRARIES ####
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(patchwork)

## STYLES ####
library(envalysis)

## SETTING WORKING DIRECTORY ####
setwd("../../Volumes/Extreme SSD/Projects/dissertation/")

## PRESETTING ####
up.color <- "#E48BAA"
down.color <- "#86B7D3"
default.color <- "#B8DCB8"

# LOAD DATA ####

## LOAD GENEINFO MATRIX ####
geneinfo <- read.csv("data/geneinfo/LncBookv2.0_GENCODEv34_GRCh38.filtered.gene.info.tsv", sep = "\t") %>% dplyr::select(c(geneid, symbol))

## LOAD PATTERN MATRIX ####
pattern.differentiation <- read.csv("data/featured/differentiation/pattern_srp168391.csv") %>% 
  dplyr::inner_join(geneinfo, by = 'geneid') %>% 
  dplyr::filter(grepl("HSALNG", geneid))

## LOAD TPM MATRIX ####
tpm.differentiation <- read.csv("data/tpm/differentiation/ge_srp168391.csv")

# PLOT ####

## GENE NUMBER OF 4 CLUSTERS ####
differentiation.pattern.stat.p <- pattern.differentiation %>% 
  dplyr::filter(str_detect(geneid, 'HSALNG')) %>% 
  dplyr::group_by(k_4) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(status = if_else(k_4 %in% c(2, 3), "up", "down")) %>% 
  ggplot(aes(x = k_4, y = count)) + 
  geom_bar(
    aes(fill = status),
    stat = 'identity',
    width = 0.5
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.6, size = 3) +
  scale_fill_manual(values = c("up" = up.color, "down" = down.color)) + 
  scale_y_continuous(limits = c(0, 200)) + 
  labs(x = 'Cluster', y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    legend.position = "none"
  )


## SPECIFY FACTOR ####
time.order <- c('zero', 'one', 'two', 'three', 'four', 
                'five', 'six', 'seven', 'eight', 'nine', 
                'ten', 'elven', 'twelve', 'thirteen', 'fourteen', 'fifteen')

## ZSCORE TRANSFORM ####
zscore.differentiation <- tpm.differentiation %>% 
  dplyr::select(-c(id)) %>% 
  dplyr::filter(geneid %in% pattern.differentiation$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(zero:fifteen, ~ (.-mean(c_across(zero:fifteen))) / sd(c_across(zero:fifteen))))


## PLOT ####

dynamic.plot <- function(zscored.df, cluster.df, levels, cluster, color, labels) {
  p <- zscored.df %>% 
    tidyr::pivot_longer(cols = -geneid, names_to = 'Day', values_to = 'TPM') %>% 
    dplyr::mutate(Day = factor(Day, levels = levels)) %>% 
    dplyr::inner_join(cluster.df %>% dplyr::select(c(geneid, k_4)), by = 'geneid') %>% 
    dplyr::filter(k_4 == cluster) %>% 
    ggplot(aes(x = Day, y = TPM, group = geneid)) + 
    geom_point(
      color = color,
      size = 1
    ) + 
    geom_line(
      color = color,
      alpha = 0.6
    ) + 
    scale_x_discrete(labels = labels) + 
    labs(x = 'Day', y = 'Normalized TPM') + 
    theme_publish() + 
    theme(
      legend.position = "none"
    )
  
  return(p)
  
}

C1.p <- dynamic.plot(
  zscored.df = zscore.differentiation, 
  cluster.df = pattern.differentiation, 
  levels = time.order, 
  cluster = 1, 
  color = down.color, 
  labels = as.character(0:15)
)

C2.p <- dynamic.plot(
  zscored.df = zscore.differentiation, 
  cluster.df = pattern.differentiation, 
  levels = time.order, 
  cluster = 2, 
  color = up.color, 
  labels = as.character(0:15)
)

C3.p <- dynamic.plot(
  zscored.df = zscore.differentiation, 
  cluster.df = pattern.differentiation, 
  levels = time.order, 
  cluster = 3, 
  color = up.color, 
  labels = as.character(0:15)
)

C4.p <- dynamic.plot(
  zscored.df = zscore.differentiation, 
  cluster.df = pattern.differentiation, 
  levels = time.order, 
  cluster = 4, 
  color = down.color, 
  labels = as.character(0:15)
)

patched.p <- (differentiation.pattern.stat.p | C1.p | C2.p) / (C3.p | C4.p | plot_spacer()) + 
  plot_layout(nrow = 4, guides = "collect") + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/differentiation.pdf", width = 10, height = 11.69)
print(patched.p)
dev.off()

## PICK SIGNIFICANT GENES ####

