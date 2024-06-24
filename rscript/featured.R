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
featured <- read.csv("data/featured/featured.csv")

featured <- featured %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(num = sum(c_across(normal:virus) == 'True'))

featured %>% 
  dplyr::summarise(count = n())
