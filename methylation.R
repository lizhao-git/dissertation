# LOAD LIBRARIES ####

## DATAFORMAT TRANSFORM ####
library(dplyr)
library(tidyr)
library(stringr)

## DATA PROCESSING ####
library(MatchIt)
library(forcats)

## PLOT LIBRARIES ####
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggvenn)
library(gghalves)
library(patchwork)

## PRESETTING ####
up.color <- "#E48BAA"
down.color <- "#86B7D3"
default.color <- "#B8DCB8"

## STYLING ####
library(envalysis)

# LOAD DATA ####

## LOAD METHYLATION MATRIX ####
methylation.profile <- read.csv("data/methylation/methylation.csv") %>% 
  dplyr::mutate(
    # 当 set 列为 TCGA 时，将 condition 列按照 _ 拆分，并将拆分后的第一个元素赋值给 set 列
    set = ifelse(set == "TCGA", sapply(str_split(condition, "_"), `[`, 1), set),
    # 将 condition 列根据条件进行修改
    condition = case_when(
      str_detect(condition, "_Normal") ~ "Normal",
      str_detect(condition, "_Cancer") ~ sapply(str_split(condition, "_"), `[`, 1),
      TRUE ~ condition
    )
  )

diff.methylation.all <- read.csv("data/methylation/methbrowse.csv")
diff.methylation.promoter <- diff.methylation.all %>% dplyr::filter(type == 'promoter')
diff.methylation.body <- diff.methylation.all %>% dplyr::filter(type == 'body')


## STAT FUNCTION ####
diff.methylation.all.lncrnas <- unique(diff.methylation.all$geneid)
diff.methylation.promoter.lncrnas <- unique(diff.methylation.promoter$geneid)
diff.methylation.body.lncrnas <- unique(diff.methylation.body$geneid)
diff.methylation.both.lncrnas <- intersect(diff.methylation.promoter.lncrnas, diff.methylation.body.lncrnas)
diff.methylation.promoter.uniq.lncrnas <- setdiff(diff.methylation.promoter.lncrnas, diff.methylation.both.lncrnas)
diff.methylation.body.uniq.lncrnas <- setdiff(diff.methylation.body.lncrnas, diff.methylation.both.lncrnas)

print(length(diff.methylation.all.lncrnas)) #TOTAL: 19543
print(length(diff.methylation.promoter.lncrnas)) #PROMOTER: 12639
print(length(diff.methylation.body.lncrnas)) #BODY: 9915
print(length(intersect(unique(diff.methylation.promoter$geneid), unique(diff.methylation.body$geneid)))) #COMMON: 3011

venn.p <- ggvenn(
  data = list(Promoter = diff.methylation.promoter.lncrnas, Body = diff.methylation.body.lncrnas),
  fill_color = c(default.color, "#E9CB7E"),
  set_name_size = 4
) + 
  theme_publish() + 
  theme_void() + 
  theme(
    
  )

## DIFFMETHYLATED IN PROMOTER ####
diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.promoter.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper > 0, Hypo > 0)

# Hyper: 2269, Hypo: 7136, Both: 224
promoter.stat.p <- data.frame(
  Group = c("Hyper", "Hypo", "Both"),
  Count = c(2269, 7136, 224)
) %>% 
  dplyr::mutate(fraction = Count / sum(Count)) %>%
  dplyr::mutate(percentage = paste0(round(fraction * 100, 1), "%")) %>% 
  dplyr::mutate(label = paste0(percentage, "\n(", Count, ")")) %>% 
  ggplot(aes(x = "", y = fraction, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(values = c("Hyper"=up.color, "Hypo" = down.color, "Both" = "#B8DCB8")) + 
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) + 
  theme_publish() + 
  theme_void() + 
  theme(
    legend.position = "top"
  )

diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.body.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper > 0, Hypo > 0)

# Hyper: 1458, Hypo: 5213, Both: 224
body.stat.p <- data.frame(
  Group = c("Hyper", "Hypo", "Both"),
  Count = c(2269, 7136, 233)
) %>% 
  dplyr::mutate(fraction = Count / sum(Count)) %>%
  dplyr::mutate(percentage = paste0(round(fraction * 100, 1), "%")) %>% 
  dplyr::mutate(label = paste0(percentage, "\n(", Count, ")")) %>% 
  ggplot(aes(x = "", y = fraction, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(values = c("Hyper"=up.color, "Hypo" = down.color, "Both" = "#B8DCB8")) + 
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) + 
  theme_publish() + 
  theme_void() + 
  theme(
    legend.position = "top"
  )

promoter.hyper.p <- diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.promoter.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper > 0, Hypo == 0) %>% 
  dplyr::group_by(Hyper) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::mutate(Hyper = as.factor(Hyper)) %>% 
  ggplot(aes(x = Hyper, y = Count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = up.color
  ) + 
  geom_text(aes(label = Count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6', '7', '8', '9')) + 
  scale_y_continuous(
    limits = c(0, 2000),
    breaks = seq(0, 2000, by = 500), 
    labels = c('0', '0.5k', '1k', '1.5k', '2k')
  ) + 
  labs(x = 'Occurrence', y = '#Hyper-methylated lncRNAs') + 
  theme_publish()

promoter.hypo.p <- diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.promoter.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper == 0, Hypo > 0) %>% 
  dplyr::group_by(Hypo) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      Hypo = c(6, 7, 8, 9),
      Count = c(0, 0, 0, 0)
    )
  ) %>% 
  dplyr::mutate(Hypo = as.factor(Hypo)) %>% 
  ggplot(aes(x = Hypo, y = Count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = down.color
  ) + 
  geom_text(aes(label = Count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6', '7', '8', '9')) + 
  scale_y_continuous(
    limits = c(0, 6000),
    breaks = seq(0, 6000, by = 1000), 
    labels = c('0', '1k', '2k', '3k', '4k', '5k', '6k')
  ) + 
  labs(x = 'Occurrence', y = '#Hypo-methylated lncRNAs') + 
  theme_publish()

body.hyper.p <- diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.body.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper > 0, Hypo == 0) %>% 
  dplyr::group_by(Hyper) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      Hyper = c(9),
      Count = c(0)
    )
  ) %>% 
  dplyr::mutate(Hyper = as.factor(Hyper)) %>% 
  ggplot(aes(x = Hyper, y = Count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = up.color
  ) + 
  geom_text(aes(label = Count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6', '7', '8', '9')) + 
  scale_y_continuous(
    limits = c(0, 1200),
    breaks = seq(0, 1200, by = 200), 
    labels = c('0', '200', '400', '600', '800', '1000', '1200')
  ) + 
  labs(x = 'Occurrence', y = '#Hyper-methylated lncRNAs') + 
  theme_publish()

body.hypo.p <- diff.methylation.all %>% 
  dplyr::filter(geneid %in% diff.methylation.body.uniq.lncrnas) %>% 
  dplyr::select(-content, -wikiindex, -id) %>% 
  dplyr::select(-type, everything(), type) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    Hyper = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hyper"),
    Hypo = sum(c_across(starts_with("autism_spectrum_disorder"):ends_with("liver_cancer")) == "hypo")
  ) %>% 
  dplyr::select(geneid, Hyper, Hypo) %>% 
  dplyr::filter(Hyper == 0, Hypo > 0) %>% 
  dplyr::group_by(Hypo) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      Hypo = c(8, 9),
      Count = c(0, 0)
    )
  ) %>% 
  dplyr::mutate(Hypo = as.factor(Hypo)) %>% 
  ggplot(aes(x = Hypo, y = Count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = down.color
  ) + 
  geom_text(aes(label = Count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6', '7', '8', '9')) + 
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 5000, by = 1000), 
    labels = c('0', '1k', '2k', '3k', '4k', '5k')
  ) + 
  labs(x = 'Occurrence', y = '#Hypo-methylated lncRNAs') + 
  theme_publish()


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
patched.p <- (venn.p + plot_spacer() + plot_layout(nrow = 1, widths = c(1,2))) /
  (promoter.stat.p + (promoter.hyper.p + promoter.hypo.p + plot_layout(axis_titles = "collect", axes = "collect")) + plot_layout(nrow = 1, widths = c(1,3))) /
  (body.stat.p + (body.hyper.p + body.hypo.p + plot_layout(axis_titles = "collect", axes = "collect")) + plot_layout(nrow = 1, widths = c(1,3))) + 
  plot_layout(nrow = 4) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/methylation.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

## BOXPLOT ####
beta.boxplot.p <- function(methylation.tbl, geneid, dataset, region, colors, wrap.length) {
  p <- methylation.tbl %>% 
    dplyr::filter(geneid == !!geneid) %>% 
    tidyr::separate_rows(value, sep = ',') %>% 
    dplyr::mutate(value = as.numeric(value)) %>% 
    dplyr::filter(set == dataset, type == region) %>% 
    dplyr::mutate(condition = fct_relevel(condition, "Normal", after = Inf)) %>% 
    ggplot(aes(x = condition, y = value)) + 
    geom_half_boxplot(
      aes(fill = condition),
      width = 0.4
    ) + 
    geom_half_dotplot(
      aes(fill = condition),
      dotsize = 1,
      position = position_nudge(x = 0.05)
    ) + 
    scale_fill_manual(values = colors) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = wrap.length)) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1), 
      labels = c('0', '0.1', '0.2', '0.3', '0.4', '0.5',
                 '0.6', '0.7', '0.8', '0.9', '1')
    ) + 
    labs(x = NULL, y = 'Beta value') + 
    theme_publish() + 
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = 'none'
    )
  
  return(p)
}

patched.p.list <- list()
counter <- 1
for (set in unique(methylation.profile$set)) {
  if (set != 'GSE79799') {
    patched.p.list[[counter]] <- beta.boxplot.p(
      methylation.tbl = methylation.profile, 
      geneid = 'HSALNG0091318', 
      dataset = set, 
      region = 'body', 
      colors = c(up.color, down.color), 
      wrap.length = 16
    )
    
    counter <- counter + 1
  }
}

patched.p <- wrap_plots(patched.p.list) + 
  plot_layout(nrow = 5, ncol = 4) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/methylation.hotair.body.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

