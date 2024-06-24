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
methylation.profile <- read.csv("data/methylation/methylation.csv")

diff.methylation.all <- read.csv("data/methylation/methbrowse.csv")
diff.methylation.promoter <- diff.methylation.all %>% dplyr::filter(type == 'promoter')
diff.methylation.body <- diff.methylation.all %>% dplyr::filter(type == 'body')


## STAT FUNCTION ####
print(length(unique(diff.methylation.all$geneid))) #TOTAL: 19543
print(length(unique(diff.methylation.promoter$geneid))) #PROMOTER: 12639
print(length(unique(diff.methylation.body$geneid))) #BODY: 9915
print(length(intersect(unique(diff.methylation.promoter$geneid), unique(diff.methylation.body$geneid)))) #COMMON: 3011

diff.methylation.promoter.stat <- diff.methylation.all %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    hyper_count = sum(c_across(where(is.character)) == "hyper"),
    hypo_count = sum(c_across(where(is.character)) == "hypo"),
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(geneid, hyper_count, hypo_count, type)) %>% 
  dplyr::filter(type == "promoter")

diff.methylation.body.stat <- diff.methylation.all %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    hyper_count = sum(c_across(where(is.character)) == "hyper"),
    hypo_count = sum(c_across(where(is.character)) == "hypo"),
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(geneid, hyper_count, hypo_count, type)) %>% 
  dplyr::filter(type == "body")

## PLOT FUNCTION ####
diff.methylation.
  


