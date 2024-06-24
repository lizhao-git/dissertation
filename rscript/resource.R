# LOAD LIBRARIES ####

## ENVIRONMENT SETTING ####
library(here)
library(envalysis)
library(devtools)

## DATAFORMAT TRANSFORM ####
library(dplyr)
library(tidyr)
library(stringr)


## PLOT LIBRARIES ####
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(UpSetR)
library(patchwork)

## PRESETTING ####
source(here("rscript", "themes.R"))

# LOAD DATA ####

## LOAD CONVERSION MATRIX ####
conversion <- read.csv("data/resource/conversion.csv")
conversion <- conversion %>% dplyr::filter(grepl("HSALNG", geneid))

resource.number <- conversion %>% 
  dplyr::select(noncode:fantom) %>%
  dplyr::summarise_all(~ sum(. != "N/A")) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = "Resource", values_to = "Gene") %>% 
  ggplot(aes(x = reorder(Resource, -Gene), y = Gene)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = '#B8DCB8'
  ) + 
  geom_text(
    aes(label = Gene), 
    size = 3,
    vjust = -0.5
  ) +
  scale_x_discrete(
    labels = c('NONCODE', 'MITranscriptome', 'LNCipedia', 'RefLnc', 'CHESS', 'GENCODE', 'BIGTranscriptome', 'FANTOM')
  ) + 
  scale_y_continuous(
    limits = c(0, 65000),
    breaks = seq(0, 60000, by = 10000), 
    labels = c('0', '10k', '20k', '30k', '40k', '50k', '60k')
  ) + 
  labs(x = NULL, y = '#LncRNAs') + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
    plot.margin = unit(c(1,1,1,1), "cm")
  )

resource.common <- conversion %>%
  dplyr::rowwise() %>%
  dplyr::mutate(count = sum(c_across(noncode:fantom) != "N/A")) %>% 
  dplyr::count(count) %>% 
  ggplot(aes(x = reorder(count, -n), y = n)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = '#B8DCB8'
  ) + 
  geom_text(
    aes(label = n), 
    size = 3,
    vjust = -0.5
  ) +
  scale_y_continuous(
    limits = c(0, 60000),
    breaks = seq(0, 60000, by = 10000), 
    labels = c('0', '10k', '20k', '30k', '40k', '50k', '60k')
  ) + 
  labs(x = 'Occurrence', y = '#LncRNAs') + 
  theme_publish() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

resource.upset <- upset(
  conversion[,-1] %>% dplyr::mutate(across(everything(), ~ ifelse(. == "N/A", 0, 1)))
  )


patched.plot <- (resource.number | resource.common) / resource.upset + 
  plot_layout(nrow = 4, ncol = 2) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/reference/reference.pdf", width = 8.27, height = 11.69)
print(patched.plot)
dev.off()

