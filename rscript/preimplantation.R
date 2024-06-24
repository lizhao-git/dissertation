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
library(ggvenn)
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
pattern.srp011546 <- read.csv("data/featured/preimplantation/pattern_srp011546.csv") %>% dplyr::inner_join(geneinfo, by = 'geneid') %>% dplyr::filter(grepl("HSALNG", geneid))
pattern.srp061636 <- read.csv("data/featured/preimplantation/pattern_srp061636.csv") %>% dplyr::inner_join(geneinfo, by = 'geneid') %>% dplyr::filter(grepl("HSALNG", geneid))

dynamic.lncrnas.srp011546 <- unique(pattern.srp011546$geneid)
dynamic.lncrnas.srp061636 <- unique(pattern.srp061636$geneid)
length(dynamic.lncrnas.srp011546) # 957
length(dynamic.lncrnas.srp061636) # 1249

venn.p <- ggvenn(
  data = list(SRP011546 = dynamic.lncrnas.srp011546, SRP061636 = dynamic.lncrnas.srp061636),
  fill_color = c(default.color, "#E9CB7E"),
  set_name_size = 4
  ) + 
  theme_publish() + 
  theme_void() + 
  theme(
    
  )

## STAT LNCRNAS IN 4 CLUSTERS FROM 2 INDENPENDENT DATASETS ####
srp011546.stat.p <- pattern.srp011546 %>% 
  dplyr::filter(str_detect(geneid, 'HSALNG')) %>% 
  dplyr::group_by(k_4) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(status = if_else(k_4 %in% c(3, 4), "up", "down")) %>% 
  ggplot(aes(x = k_4, y = count)) + 
  geom_bar(
    aes(fill = status),
    stat = 'identity',
    width = 0.5
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.6, size = 3) +
  scale_fill_manual(values = c("up" = up.color, "down" = down.color)) + 
  scale_y_continuous(limits = c(0, 700)) + 
  labs(x = 'Cluster', y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    legend.position = "none"
  )

srp061636.stat.p <- pattern.srp061636 %>% 
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
  scale_y_continuous(limits = c(0, 800)) + 
  labs(x = 'Cluster', y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    legend.position = "none"
  )

## LOAD TPM MATRIX ####
zscored.srp011546 <- read.csv("data/tpm/preimplantation/ge_srp011546.csv") %>% 
  dplyr::filter(geneid %in% pattern.srp011546$geneid) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    oocyte = mean(c_across(starts_with("oocyte"))),
    zygote = mean(c_across(starts_with("zygote"))),
    two_cell = mean(c_across(starts_with("two_cell"))),
    four_cell = mean(c_across(starts_with("four_cell"))),
    eight_cell = mean(c_across(starts_with("eight_cell"))),
    morulae = mean(c_across(starts_with("morulae"))),
    lateblastocyst = mean(c_across(starts_with("lateblastocyst")))
  ) %>% 
  dplyr::select(c(geneid, oocyte, zygote, two_cell, 
                  four_cell, eight_cell, morulae, lateblastocyst)) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(oocyte:lateblastocyst, ~ (.-mean(c_across(oocyte:lateblastocyst))) / sd(c_across(oocyte:lateblastocyst))))

zscored.srp061636 <- read.csv("data/tpm/preimplantation/ge_srp061636.csv") %>% 
  dplyr::filter(geneid %in% pattern.srp061636$geneid) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    oocyte = mean(c_across(starts_with("oocyte"))),
    zygote = mean(c_across(starts_with("zygote"))),
    two_cell = mean(c_across(starts_with("two_cell"))),
    four_cell = mean(c_across(starts_with("four_cell"))),
    eight_cell = mean(c_across(starts_with("eight_cell"))),
    morulae = mean(c_across(starts_with("morulae"))),
    earlyblastocyst = mean(c_across(starts_with("earlyblastocyst"))),
    middleblastocyst = mean(c_across(starts_with("middleblastocyst"))),
    lateblastocyst = mean(c_across(starts_with("lateblastocyst")))
  ) %>% 
  dplyr::select(c(geneid, oocyte, zygote, two_cell, 
                  four_cell, eight_cell, morulae, earlyblastocyst, 
                  middleblastocyst, lateblastocyst)) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(oocyte:lateblastocyst, ~ (.-mean(c_across(oocyte:lateblastocyst))) / sd(c_across(oocyte:lateblastocyst))))

srp011546.levels <- c('oocyte', 'zygote', 'two_cell', 'four_cell', 'eight_cell', 'morulate', 'lateblastocyst')
srp061636.levels <- c('oocyte', 'zygote', 'two_cell', 'four_cell', 'eight_cell', 'morulate', 'earlyblastocyst', 'middleblastocyst', 'lateblastocyst')
srp011546.labels <- c('Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 'Morulate', 'Late blastocyst')
srp061636.labels <- c('Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 'Morulate', 'Early blastocyst', 'Middle blastocyst', 'Late blastocyst')


dynamic.plot <- function(zscored.df, cluster.df, levels, cluster, color, labels) {
  p <- zscored.df %>% 
    tidyr::pivot_longer(cols = -geneid, names_to = 'Period', values_to = 'TPM') %>% 
    dplyr::mutate(Period = factor(Period, levels = levels)) %>% 
    dplyr::inner_join(cluster.df %>% dplyr::select(c(geneid, k_4)), by = 'geneid') %>% 
    dplyr::filter(k_4 == cluster) %>% 
    ggplot(aes(x = Period, y = TPM, group = geneid)) + 
    geom_point(
      color = color,
      size = 1
    ) + 
    geom_line(
      color = color,
      alpha = 0.6
    ) + 
    scale_x_discrete(labels = labels) + 
    labs(x = 'Period', y = 'Normalized TPM') + 
    theme_publish() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    )
  
  return(p)
  
}

srp011546.c1.p <- dynamic.plot(zscored.df = zscored.srp011546, cluster.df = pattern.srp011546, levels = srp011546.levels, cluster = 1, color = down.color, labels = srp011546.labels)
srp011546.c2.p <- dynamic.plot(zscored.df = zscored.srp011546, cluster.df = pattern.srp011546, levels = srp011546.levels, cluster = 2, color = down.color, labels = srp011546.labels)
srp011546.c3.p <- dynamic.plot(zscored.df = zscored.srp011546, cluster.df = pattern.srp011546, levels = srp011546.levels, cluster = 3, color = up.color, labels = srp011546.labels)
srp011546.c4.p <- dynamic.plot(zscored.df = zscored.srp011546, cluster.df = pattern.srp011546, levels = srp011546.levels, cluster = 4, color = up.color, labels = srp011546.labels)

srp061636.c1.p <- dynamic.plot(zscored.df = zscored.srp061636, cluster.df = pattern.srp061636, levels = srp061636.levels, cluster = 1, color = down.color, labels = srp061636.labels)
srp061636.c2.p <- dynamic.plot(zscored.df = zscored.srp061636, cluster.df = pattern.srp061636, levels = srp061636.levels, cluster = 2, color = up.color, labels = srp061636.labels)
srp061636.c3.p <- dynamic.plot(zscored.df = zscored.srp061636, cluster.df = pattern.srp061636, levels = srp061636.levels, cluster = 3, color = up.color, labels = srp061636.labels)
srp061636.c4.p <- dynamic.plot(zscored.df = zscored.srp061636, cluster.df = pattern.srp061636, levels = srp061636.levels, cluster = 4, color = down.color, labels = srp061636.labels)

patched.p <- (venn.p + (srp011546.stat.p + srp061636.stat.p + plot_layout(axes = "collect", axis_titles = "collect", ncol = 2)) + plot_spacer() + plot_layout(nrow = 1, widths = c(1,2,1))) /
(srp011546.c1.p + srp011546.c2.p + srp011546.c3.p + srp011546.c4.p + plot_layout(axes = "collect", axis_titles = "collect", nrow = 1, ncol = 4)) /
  (srp061636.c1.p + srp061636.c2.p + srp061636.c3.p + srp061636.c4.p + plot_layout(axes = "collect", axis_titles = "collect", nrow = 1, ncol = 4)) + 
  plot_layout(nrow = 4, widths = c(1,1,1,1)) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/preimplantation.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

## PICK SIGNIFICANT LNCRNAS ####
pattern.srp011546.ug <- pattern.srp011546 %>% dplyr::filter(k_4 %in% c(3,4))
pattern.srp011546.dg <- pattern.srp011546 %>% dplyr::filter(!k_4 %in% c(3,4))
pattern.srp061636.ug <- pattern.srp061636 %>% dplyr::filter(k_4 %in% c(2,3))
pattern.srp061636.dg <- pattern.srp061636 %>% dplyr::filter(!k_4 %in% c(2,3))

union.lncrnas <- union(intersect(pattern.srp011546.ug$geneid, pattern.srp061636.ug$geneid), 
      intersect(pattern.srp011546.dg$geneid, pattern.srp061636.dg$geneid))

picked.lncrnas <- geneinfo %>% 
  dplyr::filter(geneid %in% union.lncrnas) %>% 
  dplyr::filter(!grepl("HSALNG", symbol))

