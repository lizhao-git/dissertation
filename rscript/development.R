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

## STYLES ####
library(envalysis)

## PRESETTING ####
up.color <- "#E48BAA"
down.color <- "#86B7D3"
stable.color <- "#B8DCB8"
# LOAD DATA ####

## LOAD GENE INFORMATION ####
gene.info <- read.delim("data/geneinfo/LncBookv2.0_GENCODEv34_GRCh38.filtered.gene.info.tsv", sep = "\t") %>% dplyr::select(c("geneid", "symbol"))

## LOAD PATTERN MATRIX ####
development.brain <- read.csv("data/featured/development/pattern_brain.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.cerebellum <- read.csv("data/featured/development/pattern_cerebellum.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.heart <- read.csv("data/featured/development/pattern_heart.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.kidney <- read.csv("data/featured/development/pattern_kidney.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.liver <- read.csv("data/featured/development/pattern_liver.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.overy <- read.csv("data/featured/development/pattern_ovary.csv") %>% dplyr::filter(grepl("HSALNG", geneid))
development.testis <- read.csv("data/featured/development/pattern_testis.csv") %>% dplyr::filter(grepl("HSALNG", geneid))

## LOAD TPM MATRIX ####
zscored.brain <- read.csv("data/tpm/development/ge_brain.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.brain$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))

zscored.cerebellum <- read.csv("data/tpm/development/ge_cerebellum.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.cerebellum$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))

zscored.heart <- read.csv("data/tpm/development/ge_heart.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.heart$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))
  
zscored.kidney <- read.csv("data/tpm/development/ge_kidney.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.kidney$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))

zscored.liver <- read.csv("data/tpm/development/ge_liver.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.liver$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))

zscored.overy <- read.csv("data/tpm/development/ge_ovary.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.overy$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))

zscored.testis <- read.csv("data/tpm/development/ge_testis.csv") %>% 
  dplyr::select(contains("week"), contains("geneid")) %>% 
  dplyr::filter(geneid %in% development.testis$geneid) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(across(ten_week:nine_week, ~ (.-mean(c_across(ten_week:nine_week))) / sd(c_across(ten_week:nine_week))))



levels.brain <- c("four_week", "five_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                  "eighteen_week", "nineteen_week", "twenty_week")

labels.brain <- c("4", "5", "7", "8", "9", "10", "11", "12", "13", "16", "18", "19", "20")

levels.cerebellum <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week")

labels.cerebellum <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16")

levels.heart <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                  "enghteen_week", "nineteek_week")

labels.heart <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16", "18", "19")

levels.kidney <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                  "enghteen_week", "nineteek_week")

labels.kidney <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16", "18", "19")

levels.liver <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                   "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                   "enghteen_week", "nineteek_week")

labels.liver <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16", "18", "19")

levels.overy <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                  "enghteen_week")

labels.overy <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16", "18")

levels.testis <- c("four_week", "five_week", "six_week", "seven_week", "eight_week", "nine_week", 
                  "ten_week", "eleven_week", "twelve_week", "thirteen_week", "sixteen_week",
                  "enghteen_week", "nineteen_week")

labels.testis <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "16", "18", "19")

# DATA PROCESSING ####

## FEATURED DF ####
development.featured <- rbind(
  development.brain %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'brain'),
  development.cerebellum %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'cerebellum'),
  development.heart %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'heart'),
  development.kidney %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'kidney'),
  development.liver %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'liver'),
  development.overy %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'overy'),
  development.testis %>% dplyr::select(c('geneid')) %>% dplyr::mutate(organ = 'testis')
)

## STAT FEATURED GENES ####
length(unique(development.featured$geneid)) # 6,483 featured genes


development.featured.p <- development.featured %>% 
  dplyr::group_by(organ) %>% 
  dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = reorder(organ, -count), y = count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = '#B8DCB8'
  ) + 
  geom_text(
    aes(label = count), 
    size = 3,
    vjust = -0.5
  ) +
  scale_x_discrete(
    labels = c('Testis', 'Cerebellum', 'Overy', 'Brain', 'Kidney', 'Liver', 'Heart')
  ) + 
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 5000, by = 1000), 
    labels = c('0', '1k', '2k', '3k', '4k', '5k')
  ) + 
  labs(x = NULL, y = '#Dynamic lncRNAs') + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
    plot.margin = unit(c(1,1,1,1), "cm")
  )

development.featured.compared <- development.featured %>% 
  dplyr::mutate(count = 1) %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "count")

development.featured.compared[is.na(development.featured.compared)] <- 0
rownames(development.featured.compared) <- development.featured.compared$geneid
development.featured.compared <- development.featured.compared[, -1]


count_pairs <- function(df) {
  # 获取列名
  columns <- colnames(df)
  
  # 初始化空的数据框，用于存储结果
  result <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(result) <- c("Column1", "Column2", "Count")
  
  # 使用嵌套的循环遍历每对列（只考虑上三角）
  for (i in 1:length(columns)) {
    for (j in i:length(columns)) {
      count <- if (i == j) NA else sum(df[[columns[i]]] == 1 & df[[columns[j]]] == 1)
      result <- rbind(result, data.frame(Column1 = columns[i], Column2 = columns[j], Count = count))
    }
  }
  
  return(result)
}

development.compared.p <- count_pairs(development.featured.compared) %>% 
  ggplot(aes(x = Column1, y = Column2, size = Count, label = Count)) + 
  geom_point(color = '#B8DCB8') + 
  geom_text(aes(label = Count, size = 10)) + 
  scale_size_area(max_size = 12) + 
  scale_x_discrete(labels = c('Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Overy', 'Testis')) + 
  scale_y_discrete(labels = c('Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Overy', 'Testis')) + 
  labs(x = NULL, y = NULL) + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none"
  )

patched.p <- (development.featured.p | development.compared.p) + 
  plot_layout(nrow = 4, ncol = 2) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/development.compared.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

development.p.function <- function(zscored.df, cluster.df, levels, cluster, color, labels) {
  p <- zscored.df %>% 
    tidyr::pivot_longer(cols = -geneid, names_to = 'Week', values_to = 'TPM') %>% 
    dplyr::mutate(Week = factor(Week, levels = levels)) %>% 
    dplyr::inner_join(cluster.df %>% dplyr::select(c(geneid, k_4)), by = 'geneid') %>% 
    dplyr::filter(k_4 == cluster) %>% 
    ggplot(aes(x = Week, y = TPM, group = geneid)) + 
    geom_point(
      color = color,
      size = 1
    ) + 
    geom_line(
      color = color,
      alpha = 0.6
    ) + 
    scale_x_discrete(labels = labels) + 
    labs(x = 'Week', y = 'Normalized TPM') + 
    theme_publish() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  return(p)
  
}

patched.p <- list()
counter <- 1
colors <- c(
  up.color, up.color, down.color, down.color,
  down.color, up.color, up.color, down.color,
  down.color, down.color, down.color, up.color,
  up.color, down.color, down.color, down.color,
  down.color, down.color, up.color, down.color,
  down.color, up.color, up.color, up.color,
  down.color, stable.color, stable.color, up.color
            )
for (organ in c('brain', 'cerebellum', 'heart', 'kidney', 'liver', 'overy', 'testis')) {
  for (cluster in 1:4) {
    patched.p[[counter]] <- development.p.function(
    zscored.df = get(paste0("zscored.", organ)), 
    cluster.df = get(paste0("development.", organ)), 
    levels = get(paste0("levels.", organ)), 
    cluster = cluster, 
    color = colors[counter], 
    labels = get(paste0("labels.", organ))
    )
    counter <- counter + 1
  }
}

development.brain.c1 <- 
patched.p <- wrap_plots(patched.p) + 
  plot_layout(nrow = 9, ncol = 4, axis_titles = "collect")

pdf("figs/featured/development.pdf", width = 10, height = 11.69)
print(patched.p)
dev.off()

## UG/DG ####
development.brain.ug <- development.brain %>% dplyr::filter(k_4 %in% c(1,2))
development.brain.dg <- development.brain %>% dplyr::filter(!k_4 %in% c(1,2))
development.cerebellum.ug <- development.cerebellum %>% dplyr::filter(k_4 %in% c(2,3))
development.cerebellum.dg <- development.cerebellum %>% dplyr::filter(!k_4 %in% c(2,3))
development.heart.ug <- development.heart %>% dplyr::filter(k_4 %in% c(4))
development.heart.dg <- development.heart %>% dplyr::filter(!k_4 %in% c(4))
development.kidney.ug <- development.kidney %>% dplyr::filter(k_4 %in% c(1))
development.kidney.dg <- development.kidney %>% dplyr::filter(!k_4 %in% c(1))
development.liver.ug <- development.liver %>% dplyr::filter(k_4 %in% c(3))
development.liver.dg <- development.liver %>% dplyr::filter(!k_4 %in% c(3))
development.overy.ug <- development.overy %>% dplyr::filter(k_4 %in% c(2,3,4))
development.overy.dg <- development.overy %>% dplyr::filter(!k_4 %in% c(2,3,4))
development.testis.ug <- development.testis %>% dplyr::filter(k_4 %in% c(4))
development.testis.dg <- development.testis %>% dplyr::filter(!k_4 %in% c(4))

development.group <- rbind(
  development.brain.ug %>% dplyr::mutate(organ = 'Brain', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.brain.dg %>% dplyr::mutate(organ = 'Brain', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.cerebellum.ug %>% dplyr::mutate(organ = 'Cerebellum', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.cerebellum.dg %>% dplyr::mutate(organ = 'Cerebellum', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.heart.ug %>% dplyr::mutate(organ = 'Heart', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.heart.dg %>% dplyr::mutate(organ = 'Heart', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.kidney.ug %>% dplyr::mutate(organ = 'Kidney', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.kidney.dg %>% dplyr::mutate(organ = 'Kidney', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.liver.ug %>% dplyr::mutate(organ = 'Liver', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.liver.dg %>% dplyr::mutate(organ = 'Liver', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.overy.ug %>% dplyr::mutate(organ = 'Overy', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.overy.dg %>% dplyr::mutate(organ = 'Overy', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.testis.ug %>% dplyr::mutate(organ = 'Testis', group = 'UG') %>% dplyr::select(c('geneid', 'organ', 'group')),
  development.testis.dg %>% dplyr::mutate(organ = 'Testis', group = 'DG') %>% dplyr::select(c('geneid', 'organ', 'group'))
) 

ug.dg.p <- development.group %>% 
  dplyr::group_by(organ, group) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(count = if_else(group =='DG', -count, count)) %>% 
  ggplot(aes(x = organ, y = count, fill = group)) + 
  geom_bar(
    aes(fill = group),
    stat = 'identity',
    width = 0.4
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_fill_manual(values = c("UG" = up.color, "DG" = down.color)) + 
  scale_y_continuous(
    limits = c(-6000, 2000),
    breaks = seq(-6000, 2000, by = 1000), 
    labels = c('6k', '5k', '4k', '3k', '2k', '1k', '0', '1k', '2k')
    ) + 
  labs(x = NULL, y = '#LncRNA gene') + 
  theme_publish() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

## COMMON LNCRNAS IN DIFFERENT ORGANS ####
development.group %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "group") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    UG = sum(c_across(Brain:Testis) == "UG", na.rm = TRUE),
    DG = sum(c_across(Brain:Testis) == "DG", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(UG == 0, DG > 0)

# 1087 UG, 5198 DG, 198 BOTH
group.stat.p <- data.frame(
  Group = c("UG", "DG", "Both"),
  Count = c(1087, 5198, 198)
) %>% 
  dplyr::mutate(fraction = Count / sum(Count)) %>%
  dplyr::mutate(percentage = paste0(round(fraction * 100, 1), "%")) %>% 
  dplyr::mutate(label = paste0(percentage, "\n(", Count, ")")) %>% 
  ggplot(aes(x = "", y = fraction, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(values = c("UG"=up.color, "DG" = down.color, "Both" = "#B8DCB8")) + 
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) + 
  theme_publish() + 
  theme_void() + 
  theme(
    legend.position = "top"
  )


ug.p <- development.group %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "group") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    UG = sum(c_across(Brain:Testis) == "UG", na.rm = TRUE),
    DG = sum(c_across(Brain:Testis) == "DG", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(UG > 0, DG == 0) %>% 
  dplyr::group_by(UG) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::bind_rows(
    data.frame(
      UG = c(3, 4, 5, 6),
      count = c(0, 0, 0, 0)
    )
  ) %>% 
  dplyr::mutate(UG = as.factor(UG)) %>% 
  ggplot(aes(x = UG, y = count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = up.color
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.5, size = 3) +
  scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6')) + 
  scale_y_continuous(
    limits = c(0, 1500),
    breaks = seq(0, 1500, by = 500), 
    labels = c('0', '0.5k', '1k', '1.5k')
  ) + 
  labs(x = 'Occurrence', y = '#UG lncRNAs') + 
  theme_publish()

dg.p <- development.group %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "group") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    UG = sum(c_across(Brain:Testis) == "UG", na.rm = TRUE),
    DG = sum(c_across(Brain:Testis) == "DG", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(DG > 0, UG == 0) %>% 
  dplyr::group_by(DG) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(DG = as.factor(DG)) %>% 
  ggplot(aes(x = DG, y = count)) + 
  geom_bar(
    stat = 'identity',
    width = 0.6,
    fill = down.color
  ) + 
  geom_text(aes(label = count), vjust = -0.4, hjust = 0.5, size = 3) + 
  scale_x_discrete(label = c(1,2,3,4,5,6)) + 
  scale_y_continuous(
    limits = c(0, 6000),
    breaks = seq(0, 6000, by = 1000), 
    labels = c('0', '1k', '2k', '3k', '4k', '5k', '6k')
  ) + 
  labs(x = 'Occurrence', y = '#DG lncRNAs') + 
  theme_publish()

ug.dg.compared.p <- development.group %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "group") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    UG = sum(c_across(Brain:Testis) == "UG", na.rm = TRUE),
    DG = sum(c_across(Brain:Testis) == "DG", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(DG > 0, UG > 0) %>% 
  dplyr::group_by(UG,DG) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(UG = as.factor(UG), DG = as.factor(DG)) %>% 
  ggplot(aes(x = UG, y = DG, size = count, label = count)) + 
  geom_point(color = '#B8DCB8') + 
  geom_text(aes(label = count, size = 10)) + 
  scale_size_area(max_size = 16) + 
  labs(x = "Occurrence", y = "Occurrence") + 
  theme_publish() + 
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )

patched.p <- (ug.dg.p | group.stat.p) / (ug.p | dg.p | ug.dg.compared.p) + 
  plot_layout(nrow = 4) + 
  plot_annotation(tag_levels = 'A')

pdf("figs/featured/development.group.pdf", width = 8.27, height = 11.69)
print(patched.p)
dev.off()

## PICK SIGNIFICANT GENES ####
picked.genes.tb <- development.group %>% 
  tidyr::pivot_wider(names_from = "organ", values_from = "group") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    UG = sum(c_across(Brain:Testis) == "UG", na.rm = TRUE),
    DG = sum(c_across(Brain:Testis) == "DG", na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(UG + DG >= 5) %>% 
  dplyr::inner_join(gene.info, by = "geneid")

