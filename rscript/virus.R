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

length(unique(de.combined$geneid))

virus.de.p <- de.combined %>% 
  dplyr::group_by(disease, status) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(count = if_else(status =='down', -count, count)) %>% 
  ggplot(aes(x= disease, y = count)) + 
  geom_bar(
    aes(fill = status),
    stat = 'identity',
    width = 0.5
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.6, size = 3) +
  scale_fill_manual(values = c("up" = up.color, "down" = down.color)) + 
  #scale_x_discrete(labels = c('Coronary', 'ESCC', 'Pancreatic', 'Hepatocellular', 'Colorectal')) + 
  scale_y_continuous(
    limits = c(-600, 500),
    breaks = seq(-600, 500, by = 100),
    labels = c('-600', '-500', '-400', '-300', '-200', '-100', '0', '100', '200', '300', '400', '500')
  ) + 
  labs(x = NULL, y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# Virus: 445 UG, 536 DG, 4 BOTH
# Exosome: 380 UG, 1151 DG, 7 BOTH

virus.stat.p <- data.frame(
  Group = c("Up-regulated", "Down-regulated", "Both"),
  Count = c(445, 536, 4)
) %>% 
  dplyr::mutate(fraction = Count / sum(Count)) %>%
  dplyr::mutate(percentage = paste0(round(fraction * 100, 1), "%")) %>% 
  dplyr::mutate(label = paste0(percentage, "\n(", Count, ")")) %>% 
  ggplot(aes(x = "", y = fraction, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(values = c("Up-regulated"=up.color, "Down-regulated" = down.color, "Both" = "#B8DCB8")) + 
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) + 
  theme_publish() + 
  theme_void() + 
  theme(
    legend.position = "top"
  )

virus.down.p <- de.combined %>% 
  dplyr::select(c('geneid', 'disease', 'status')) %>% 
  tidyr::pivot_wider(names_from = 'disease', values_from = 'status') %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    Down = sum(c_across(HIV:COVID) == "down", na.rm = TRUE),
    Up = sum(c_across(HIV:COVID) == "up", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(Up == 0, Down > 0) %>% 
  dplyr::group_by(Down) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      Down = c(3, 4),
      count = c(0, 0)
    )
  ) %>% 
  dplyr::mutate(Down = as.factor(Down)) %>% 
  ggplot(aes(x = Down, y = count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = down.color
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_y_continuous(
    limits = c(0, 600),
    breaks = seq(0, 600, by = 100), 
    labels = c('0', '100', '200', '300', '400', '500', '600')
  ) + 
  labs(x = 'Occurrence', y = '#Down-regulated lncRNAs') + 
  theme_publish()

patched.p <- (virus.de.p + exosome.de.p + plot_layout(axis_titles = "collect", axes = "collect")) /
  (virus.stat.p + (virus.up.p + virus.down.p + plot_layout(axis_titles = "collect", axes = "collect")) + plot_layout(nrow = 1, widths = c(1,2))) / 
  (exosome.stat.p + (exosome.up.p + exosome.down.p + plot_layout(axis_titles = "collect", axes = "collect")) + plot_layout(nrow = 1, widths = c(1,2))) + 
  plot_layout(nrow = 4) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/virus.exosome.de.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()


## PICK SIGNIFICANT LNCRNAS ####
picked.lncrnas <- de.combined %>% 
  dplyr::group_by(disease) %>% 
  dplyr::arrange(desc(log2foldchange)) %>% 
  dplyr::slice(c(1:3, (n()-2):n())) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct() %>% 
  tidyr::pivot_wider(names_from = 'disease', values_from = 'log2foldchange')

write_xlsx(picked.lncrnas, "data/featured/virus/picked.lncrnas.xlsx")
