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
coexpression <- read.csv("data/coexpression/coexpression.csv")


# PLOT ####
coexpression %>% 
  dplyr::count(breadth) %>% 
  dplyr::mutate(breadth = factor(breadth)) %>% 
  ggplot(aes(x = '', y = n, fill = breadth)) + 
  geom_bar(stat = 'identity', position = 'stack', width = 1) + 
  coord_polar(theta = 'y') + 
  theme_void() + 
  theme(legend.title = element_blank())


coexpression %>% 
  dplyr::filter()
