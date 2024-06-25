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

## DATA READER/WRITER ####
library(writexl)

# LOAD DATA ####

## LOAD DIFF MATRIX ####
de.colorectal <- read.csv("data/featured/exosome/de_colorectal_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Colorecatal cancer') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.coronary <- read.csv("data/featured/exosome/de_coronary_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Coronary heart disease') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.escc <- read.csv("data/featured/exosome/de_escc_control.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Eearly-stage esophageal squamous cell carcinoma') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.esophagitis <- read.csv("data/featured/exosome/de_esophagitis_control.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Colorecatal cancer') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.hepatocellular <- read.csv("data/featured/exosome/de_hepatocellular_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Hepatocellular carcinoma') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.pancreatic <- read.csv("data/featured/exosome/de_pancreatic_normal.csv") %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::filter(padj<=0.05, abs(log2foldchange) >= 1) %>% 
  dplyr::mutate(disease = 'Pancreatic adenocarcinoma') %>% 
  dplyr::mutate(status = if_else(log2foldchange >=1, "up", "down"))

de.combined <- rbind(de.colorectal, de.coronary, de.escc, de.hepatocellular, de.pancreatic)

exosome.de.p <- de.combined %>% 
  dplyr::group_by(disease, status) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(count = if_else(status =='down', -count, count)) %>% 
  ggplot(aes(x= reorder(disease, count), y = count)) + 
  geom_bar(
    aes(fill = status),
    stat = 'identity',
    width = 0.6
    ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.6, size = 3) +
  scale_fill_manual(values = c("up" = up.color, "down" = down.color)) + 
  scale_x_discrete(labels = c('Coronary', 'ESCC', 'Pancreatic', 'Hepatocellular', 'Colorectal')) + 
  scale_y_continuous(limits = c(-1000, 500)) + 
  labs(x = NULL, y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )


exosome.stat.p <- data.frame(
  Group = c("Up-regulated", "Down-regulated", "Both"),
  Count = c(380, 1151, 7)
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

exosome.down.p <- de.combined %>% 
  dplyr::select(c('geneid', 'disease', 'status')) %>% 
  tidyr::pivot_wider(names_from = 'disease', values_from = 'status') %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    Down = sum(c_across('Colorecatal cancer':'Pancreatic adenocarcinoma') == "down", na.rm = TRUE),
    Up = sum(c_across('Colorecatal cancer':'Pancreatic adenocarcinoma') == "up", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(Up == 0, Down > 0) %>% 
  dplyr::group_by(Down) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      Down = c(3, 4, 5),
      count = c(0, 0, 0)
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
    limits = c(0, 1200),
    breaks = seq(0, 1200, by = 200), 
    labels = c('0', '200', '400', '600', '800', '1000', '1200')
  ) + 
  labs(x = 'Occurrence', y = '#Down-regulated lncRNAs') + 
  theme_publish()

## PICK SIGNIFICANT LNCRNAS ####
picked.lncrnas <- de.combined %>% 
  dplyr::group_by(disease) %>% 
  dplyr::arrange(desc(log2foldchange)) %>% 
  dplyr::slice_head(n = 5) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = 'disease', values_from = 'log2foldchange')

write_xlsx(picked.lncrnas, "data/featured/exosome/picked.lncrnas.xlsx")


