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
library(UpSetR)

# LOAD DATA ####

## LOAD METHYLATION MATRIX ####
apex <- read.csv("data/featured/apex/apex.csv") %>% dplyr::select(-c('kedl'))

location.breadth <- apex %>% 
  dplyr::count(breadth) %>% 
  dplyr::rename(breadth = 1, count  = 2)

## STAT ####
location.preference <- apex %>% 
  dplyr::summarise(across(erm:omm, ~ sum(. == "Yes"), .names = "count_{col}"))

mitochondrion <- apex %>% 
  dplyr::filter(mito == 'Yes' | omm == 'Yes')

er <- apex %>% 
  dplyr::filter(erm == 'Yes' | kdel == 'Yes')

# PLOT ####
occurrence.p <- apex %>% 
  dplyr::count(breadth) %>% 
  dplyr::rename(Occurrence = 1, Count = 2) %>% 
  ggplot(aes(x = Occurrence, y = Count)) + 
  geom_bar(
    fill = '#B8DCB8',
    stat = 'identity',
    width = 0.6
    ) + 
  geom_text(
    aes(label = Count), 
    size = 3,
    vjust = -0.5
  ) +
  labs(x = 'Occurrence', y = '#LncRNA') + 
  theme_publish() + 
  theme(
    legend.position = "none"
  )

location.p <- apex %>% 
  dplyr::summarise(across(erm:omm, ~ sum(. == "Yes"))) %>% 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Localization",
                      values_to = "Count") %>% 
  ggplot(aes(x = reorder(Localization, -Count), y = Count)) + 
  geom_bar(
    fill = '#B8DCB8',
    stat = 'identity',
    width = 0.6
  ) + 
  geom_text(
    aes(label = Count), 
    size = 3,
    vjust = -0.5
  ) +
  labs(x = NULL, y = '#LncRNA') + 
  scale_x_discrete(labels = c('Nucleus', 'Nuclear lamina', 'Nucleolus', 'ER membrane', 'Mitochondrial matrix', 
                              'Outer Mito. Membrane', 'Cytosol', 'ER lumen', 'Nuclear pore')) +
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

## UPSETR PLOT ####
upset(apex %>% 
  dplyr::mutate(across(everything(), ~ ifelse(. == "Yes", 1, 0))) %>% 
  dplyr::select(-c(geneid, genename, breadth)),
  order.by = "freq",
  nintersects = 10,
  mb.ratio = c(0.5, 0.5),
  point.size = 1,
  line.size = 1, 
  main.bar.color = "#2a83a2", 
  sets.bar.color = "#3b7960"
  )

patched.p <- occurrence.p + location.p + 
  plot_layout(nrow = 4, ncol = 2) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/apex.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

## MITRO PLOT ####
