# LOAD LIBRARIES ####

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

## PLOT STYLES ####
library(envalysis)

# LOAD DATA ####

## LOAD METHYLATION MATRIX ####
capacity <- read.csv("data/capacity/capacity.csv")

capacity <- capacity %>% 
  dplyr::select(-c(none, low, medium, high)) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    HC_count = sum(c_across(normal:virus) == "HC"),
    MC_count = sum(c_across(normal:virus) == "MC"),
    LC_count = sum(c_across(normal:virus) == "LC"),
    NE_count = sum(c_across(normal:virus) == "NE")
  ) %>%
  dplyr::ungroup()

capacity