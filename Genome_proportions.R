rm(list = ls())

library(tidyverse)
library(ggpubr)
library(forcats)

#### Nivalis ####
## Load datasets
nivFiles_run1 <- list.files(path = "input/niv", pattern = "*_run1.txt",
                            full.names = TRUE)
for (f in 1:length(nivFiles_run1)) {
  df_niv <- read.delim(nivFiles_run1[f])
  # Fuse All and All/repeat together
  df_niv$Annotation[df_niv$Annotation == 'All'] <- 'Unidentified repeat'
  df_niv$Annotation[df_niv$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df_niv$Annotation[df_niv$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df_niv$Annotation[grep("All/repeat/rDNA/*", df_niv$Annotation)] <-
    'All/repeat/rDNA'
  df_niv$Proportion... <- as.numeric(df_niv$Proportion...)
  df_niv <- aggregate(Proportion... ~ Annotation, data = df_niv, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop_niv <- c("Non repetitive", 50 - sum(df_niv$Proportion...))
  df_niv <- rbind(df_niv, non_repeat_prop_niv)
  # Add a column with individual name
  column_name_niv <- paste0(str_match(nivFiles_run1[f],
                                      "Genome_proportion_(.*?)_run1.txt")[, 2])
  df_niv$Indiv <- column_name_niv
  # Set file name
  dataframe_name_niv <- str_match(nivFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name_niv, df_niv)
}

## Merge datasets
nivList <- lapply(ls(pattern = "input/niv/Genome_proportion*"), get)
nivFull <- nivList %>%
  reduce(full_join, by = c('Annotation', 'Indiv', 'Proportion...'))
nivFull$Proportion... <- as.numeric(nivFull$Proportion...)
## Rename REs
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Bel-Pao'] <- 'TE/Class_I/LTR/Bel-Pao'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'] <- 'TE/Class_I/LTR/Ty1_copia'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_I/Penelope'] <- 'TE/Class_I/Penelope'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Helitron'] <- 'TE/Class_II/Helitron'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Maverick'] <- 'TE/Class_II/Maverick'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
nivFull$Annotation[nivFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
## Change levels to change order of REs in plot
nivFull <- nivFull %>%
  mutate(Annotation = factor(Annotation, levels = c("TE/Class_I/LINE",
                                                    "TE/Class_I/LTR",
                                                    "TE/Class_I/LTR/Bel-Pao",
                                                    "TE/Class_I/LTR/Ty1_copia",
                                                    "TE/Class_I/LTR/Ty3_gypsy",
                                                    "TE/Class_I/Penelope",
                                                    "TE/Class_II/Helitron",
                                                    "TE/Class_II/Maverick",
                                                    "rDNA",
                                                    "satDNA",
                                                    "Unidentified repeat",
                                                    "Non repetitive")))

## Remove 1 sibling
nivFull <- nivFull[!(nivFull$Indiv == "saj_Saj10"),]

## Rename samples
nivFull$Indiv[nivFull$Indiv == 'gri_X0492'] <- 'GRI_1'
nivFull$Indiv[nivFull$Indiv == 'gri_X0330'] <- 'GRI_2'
nivFull$Indiv[nivFull$Indiv == 'gri_Z0016'] <- 'GRI_3'
nivFull$Indiv[nivFull$Indiv == 'gro_Glo05'] <- 'GRO_1'
nivFull$Indiv[nivFull$Indiv == 'gro_Glo08'] <- 'GRO_2'
nivFull$Indiv[nivFull$Indiv == 'gro_Glo11'] <- 'GRO_3'
nivFull$Indiv[nivFull$Indiv == 'saj_Saj04'] <- 'SAJ_1'
nivFull$Indiv[nivFull$Indiv == 'saj_Saj08'] <- 'SAJ_2'
nivFull$Indiv[nivFull$Indiv == 'sch_X1841'] <- 'SCH_1'
nivFull$Indiv[nivFull$Indiv == 'sch_X1819'] <- 'SCH_2'
nivFull$Indiv <- fct_rev(nivFull$Indiv)

## Plot genome proportions
nivPlot <- nivFull %>%
  ggplot(aes(x = Indiv, y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(c) ", italic("Erebia nivalis")))) +
  scale_y_continuous(limits = c(0, 50.001), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4",
                               "darkslategrey", "seagreen", "#b15928",
                               "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))
nivPlot

#### Tyndarus ####
## Load datasets
tynFiles_run1 <- list.files(path = "input/tyn", pattern = "*_run1.txt",
                            full.names = TRUE)
for (f in 1:length(tynFiles_run1)) {
  df_tyn <- read.delim(tynFiles_run1[f])
  # Fuse All and All/repeat together
  df_tyn$Annotation[df_tyn$Annotation == 'All'] <- 'Unidentified repeat'
  df_tyn$Annotation[df_tyn$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df_tyn$Annotation[df_tyn$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df_tyn$Annotation[grep("All/repeat/rDNA/*", df_tyn$Annotation)] <-
    'All/repeat/rDNA'
  df_tyn$Proportion... <- as.numeric(df_tyn$Proportion...)
  df_tyn <- aggregate(Proportion... ~ Annotation, data = df_tyn, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop_tyn <- c("Non repetitive", 50 - sum(df_tyn$Proportion...))
  df_tyn <- rbind(df_tyn, non_repeat_prop_tyn)
  # Add a column with individual name
  column_name_tyn <- paste0(str_match(tynFiles_run1[f],
                                      "Genome_proportion_(.*?)_run1.txt")[, 2])
  df_tyn$Indiv <- column_name_tyn
  # Set file name
  dataframe_name_tyn <- str_match(tynFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name_tyn, df_tyn)
}

## Merge datasets
tynList <- lapply(ls(pattern = "input/tyn/Genome_proportion*"), get)
tynFull <- tynList %>%
  reduce(full_join, by = c('Annotation', 'Indiv', 'Proportion...'))
tynFull$Proportion... <- as.numeric(tynFull$Proportion...)
## Rename REs
tynFull$Annotation[tynFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Bel-Pao'] <- 'TE/Class_I/LTR/Bel-Pao'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/mobile_element/Class_I/Penelope'] <- 'TE/Class_I/Penelope'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
tynFull$Annotation[tynFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
tynFull <- tynFull %>%
  mutate(Annotation = factor(Annotation, levels=c("TE/Class_I/LINE",
                                                  "TE/Class_I/LTR",
                                                  "TE/Class_I/LTR/Bel-Pao",
                                                  "TE/Class_I/LTR/Ty3_gypsy",
                                                  "TE/Class_I/Penelope",
                                                  "rDNA",
                                                  "satDNA",
                                                  "Unidentified repeat",
                                                  "Non repetitive")))

## Remove 1 sibling
tynFull <- tynFull[!(tynFull$Indiv == "grw_X2129"),]

# Rename samples
tynFull$Indiv[tynFull$Indiv == 'lau_X1470'] <- 'LAU_1'
tynFull$Indiv[tynFull$Indiv == 'lau_X1483'] <- 'LAU_2'
tynFull$Indiv[tynFull$Indiv == 'lau_X1488'] <- 'LAU_3'
tynFull$Indiv[tynFull$Indiv == 'got_X1707'] <- 'GOT_1'
tynFull$Indiv[tynFull$Indiv == 'got_X1721'] <- 'GOT_2'
tynFull$Indiv[tynFull$Indiv == 'got_X1741'] <- 'GOT_3'
tynFull$Indiv[tynFull$Indiv == 'aro_X1971'] <- 'ARO_1'
tynFull$Indiv[tynFull$Indiv == 'aro_X2006'] <- 'ARO_2'
tynFull$Indiv[tynFull$Indiv == 'aro_X2009'] <- 'ARO_3'
tynFull$Indiv[tynFull$Indiv == 'grw_X2124'] <- 'GRW_1'
tynFull$Indiv[tynFull$Indiv == 'grw_X2125'] <- 'GRW_2'
tynFull$Indiv[tynFull$Indiv == 'com_X2193'] <- 'COM_1'
tynFull$Indiv[tynFull$Indiv == 'com_X2198'] <- 'COM_2'
tynFull$Indiv[tynFull$Indiv == 'com_X2231'] <- 'COM_3'
tynFull$Indiv[tynFull$Indiv == 'san_X2264'] <- 'SAN_1'
tynFull$Indiv[tynFull$Indiv == 'san_X2268'] <- 'SAN_2'
tynFull$Indiv[tynFull$Indiv == 'san_X2274'] <- 'SAN_3'
tynFull$Indiv[tynFull$Indiv == 'grf_X2965'] <- 'GRF_1'
tynFull$Indiv[tynFull$Indiv == 'grf_X3029'] <- 'GRF_2'
tynFull$Indiv[tynFull$Indiv == 'grf_X3037'] <- 'GRF_3'
tynFull$Indiv <- fct_rev(tynFull$Indiv)

## Plot genome proportions
tynPlot <- tynFull %>%
  ggplot(aes(x = Indiv, y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(b) ", italic("Erebia tyndarus")))) +
  scale_y_continuous(limits = c(0, 50.001), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4",
                               "darkslategrey", "seagreen", "#b15928",
                               "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))
tynPlot

#### Cassioides ####
## Load datasets
casFiles_run1 <- list.files(path = "input/cas", pattern = "*_run1.txt",
                            full.names = TRUE)
for (f in 1:length(casFiles_run1)) {
  df <- read.delim(casFiles_run1[f])
  # Fuse All and All/repeat together
  df$Annotation[df$Annotation == 'All'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df$Annotation[grep("All/repeat/rDNA/*", df$Annotation)] <-
    'All/repeat/rDNA'
  df$Proportion... <- as.numeric(df$Proportion...)
  df <- aggregate(Proportion... ~ Annotation, data = df, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop <- c("Non repetitive", 50 - sum(df$Proportion...))
  df <- rbind(df, non_repeat_prop)
  # Add a column with individual name
  column_name <- paste0(str_match(casFiles_run1[f],
                                  "Genome_proportion_(.*?)_run1.txt")[, 2])
  df$Indiv <- column_name
  # Set file name
  dataframe_name <- str_match(casFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name, df)
}

## Merge datasets
casList <- lapply(ls(pattern = "input/cas/Genome_proportion*"), get)
casFull <- casList %>%
  reduce(full_join, by = c('Annotation', 'Indiv', 'Proportion...'))
casFull$Proportion... <- as.numeric(casFull$Proportion...)
## Rename REs
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/DIRS'] <- 'TE/Class_I/DIRS'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Bel-Pao'] <- 'TE/Class_I/LTR/Bel-Pao'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'] <- 'TE/Class_I/LTR/Ty1_copia'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
casFull$Annotation[casFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Maverick'] <- 'TE/Class_II/Maverick'
casFull$Annotation[casFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
casFull$Annotation[casFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
casFull <- casFull %>%
  mutate(Annotation = factor(Annotation, levels=c("TE/Class_I/DIRS",
                                                  "TE/Class_I/LINE",
                                                  "TE/Class_I/LTR",
                                                  "TE/Class_I/LTR/Bel-Pao",
                                                  "TE/Class_I/LTR/Ty1_copia",
                                                  "TE/Class_I/LTR/Ty3_gypsy",
                                                  "TE/Class_II/Maverick",
                                                  "rDNA",
                                                  "satDNA",
                                                  "Unidentified repeat",
                                                  "Non repetitive")))

# Rename samples
casFull$Indiv[casFull$Indiv == 'grf_X2856'] <- 'GRF_1'
casFull$Indiv[casFull$Indiv == 'grf_X2893'] <- 'GRF_2'
casFull$Indiv[casFull$Indiv == 'grf_X2894'] <- 'GRF_3'
casFull$Indiv[casFull$Indiv == 'grw_X2054'] <- 'GRW_1'
casFull$Indiv[casFull$Indiv == 'grw_X2454'] <- 'GRW_2'
casFull$Indiv[casFull$Indiv == 'grw_X2601'] <- 'GRW_3'
casFull$Indiv[casFull$Indiv == 'kan_X1906'] <- 'KAN_1'
casFull$Indiv[casFull$Indiv == 'kan_X1920'] <- 'KAN_2'
casFull$Indiv[casFull$Indiv == 'kan_X1945'] <- 'KAN_3'
casFull$Indiv[casFull$Indiv == 'pdm_X2707'] <- 'PDM_1'
casFull$Indiv[casFull$Indiv == 'pdm_X2715'] <- 'PDM_2'
casFull$Indiv[casFull$Indiv == 'pdm_X2729'] <- 'PDM_3'
casFull$Indiv[casFull$Indiv == 'rou_X2180'] <- 'ROU_1'
casFull$Indiv[casFull$Indiv == 'rou_X2305'] <- 'ROU_2'
casFull$Indiv[casFull$Indiv == 'rou_X2312'] <- 'ROU_3'
casFull$Indiv[casFull$Indiv == 'sch_X1825'] <- 'SCH_1'
casFull$Indiv[casFull$Indiv == 'sch_X1844'] <- 'SCH_2'
casFull$Indiv[casFull$Indiv == 'sch_X2033'] <- 'SCH_3'
casFull$Indiv[casFull$Indiv == 'sto_X1246'] <- 'STO_1'
casFull$Indiv[casFull$Indiv == 'sto_X1774'] <- 'STO_2'
casFull$Indiv[casFull$Indiv == 'sto_X1800'] <- 'STO_3'
casFull$Indiv <- fct_rev(casFull$Indiv)

## Plot genome proportions
casPlot <- casFull %>%
  ggplot(aes(x = Indiv, y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(a) ", italic("Erebia cassioides")))) +
  scale_y_continuous(limits = c(0, 50.001), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4",
                               "darkslategrey", "seagreen", "#b15928",
                               "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))
casPlot

#### Pronoe ####
## Load datasets
proFiles_run1 <- list.files(path = "input/pro", pattern = "*_run1.txt",
                            full.names = TRUE)
for (f in 1:length(proFiles_run1)) {
  df <- read.delim(proFiles_run1[f])
  # Fuse All and All/repeat together
  df$Annotation[df$Annotation == 'All'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element/Class_I'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df$Annotation[grep("All/repeat/rDNA/*", df$Annotation)] <- 'All/repeat/rDNA'
  df$Proportion... <- as.numeric(df$Proportion...)
  df <- aggregate(Proportion... ~ Annotation, data = df, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop <- c("Non repetitive", 50 - sum(df$Proportion...))
  df <- rbind(df, non_repeat_prop)
  # Add a column with individual name
  column_name <- paste0(str_match(proFiles_run1[f],
                                  "Genome_proportion_(.*?)_run1.txt")[, 2])
  df$Indiv <- column_name
  # Set file name
  dataframe_name <- str_match(proFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name, df)
}

## Merge datasets
proList <- lapply(ls(pattern = "input/pro/Genome_proportion*"), get)
proFull <- proList %>%
  reduce(full_join, by = c('Annotation', 'Indiv', 'Proportion...'))
proFull$Proportion... <- as.numeric(proFull$Proportion...)
## Rename REs
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Bel-Pao'] <- 'TE/Class_I/LTR/Bel-Pao'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'] <- 'TE/Class_I/LTR/Ty1_copia'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_I/Penelope'] <- 'TE/Class_I/Penelope'
proFull$Annotation[proFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Maverick'] <- 'TE/Class_II/Maverick'
proFull$Annotation[proFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
proFull$Annotation[proFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
proFull <- proFull %>%
  mutate(Annotation = factor(Annotation, levels = c("TE/Class_I/LINE",
                                                    "TE/Class_I/LTR",
                                                    "TE/Class_I/LTR/Bel-Pao",
                                                    "TE/Class_I/LTR/Ty1_copia",
                                                    "TE/Class_I/LTR/Ty3_gypsy",
                                                    "TE/Class_I/Penelope",
                                                    "TE/Class_II/Maverick",
                                                    "rDNA",
                                                    "satDNA",
                                                    "Unidentified repeat",
                                                    "Non repetitive")))


# Rename samples
proFull$Indiv[proFull$Indiv == 'psa_00058'] <- 'PSA_1'
proFull$Indiv[proFull$Indiv == 'psa_00076'] <- 'PSA_2'
proFull$Indiv[proFull$Indiv == 'psa_00101'] <- 'PSA_3'
proFull$Indiv[proFull$Indiv == 'psa_00056'] <- 'PSA_4'
proFull$Indiv[proFull$Indiv == 'psa_00064'] <- 'PSA_5'
proFull$Indiv[proFull$Indiv == 'ver_XDB01'] <- 'VER_1'
proFull$Indiv[proFull$Indiv == 'ver_A0049'] <- 'VER_2'
proFull$Indiv[proFull$Indiv == 'ver_00053'] <- 'VER_3'
proFull$Indiv[proFull$Indiv == 'ver_00057'] <- 'VER_4'
proFull$Indiv[proFull$Indiv == 'ver_00188'] <- 'VER_5'
proFull$Indiv[proFull$Indiv == 'ver_00221'] <- 'VER_6'
proFull$Indiv <- fct_rev(proFull$Indiv)

## Plot genome proportions
proPlot <- proFull %>%
  ggplot(aes(x = Indiv, y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(d) ", italic("Erebia pronoe")))) +
  scale_y_continuous(limits = c(0, 50.001), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4",
                               "darkslategrey", "seagreen", "#b15928",
                               "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))
proPlot

#### Plot full ####
pdf("output/Fig_Genome_proportions_pops.pdf", width = 6.65354)
ggarrange(casPlot, tynPlot, nivPlot, proPlot, ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom", heights = c(1.7,1))
dev.off()

#### Erebia ####
rm(list = ls())
## Load datasets
ereFiles_run1 <- list.files(path = "input/ere", pattern = "*_run1.txt",
                            full.names = TRUE)
for (f in 1:length(ereFiles_run1)) {
  df <- read.delim(ereFiles_run1[f])
  # Fuse All and All/repeat together
  df$Annotation[df$Annotation == 'All'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element/Class_I'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df$Annotation[grep("All/repeat/rDNA/*", df$Annotation)] <- 'All/repeat/rDNA'
  df$Proportion... <- as.numeric(df$Proportion...)
  df <- aggregate(Proportion... ~ Annotation, data = df, FUN = sum)
  # Fuse DNA transposons together (just for calculation of proportion)
  #df$Annotation[grep("All/repeat/mobile_element/Class_II/*", df$Annotation)] <- 'DNAtransposons'
  # Fuse retrotransposons together (just for calculation of proportion)
  #df$Annotation[grep("All/repeat/mobile_element/Class_I/*", df$Annotation)] <- 'retrotransposons'
  #df$Proportion... <- as.numeric(df$Proportion...)
  #df <- aggregate(Proportion... ~ Annotation, data = df, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop <- c("Non repetitive", 50 - sum(df$Proportion...))
  df <- rbind(df, non_repeat_prop)
  # Add a column with individual name
  column_name <- paste0(str_match(ereFiles_run1[f],
                                  "Genome_proportion_(.*?)_run1.txt")[, 2])
  df$Indiv <- column_name
  # Set file name
  dataframe_name <- str_match(ereFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name, df)
}

## Merge datasets
ereList <- lapply(ls(pattern = "input/ere/Genome_proportion*"), get)
ereFull <- ereList %>%
  reduce(full_join, by = c('Annotation', 'Indiv', 'Proportion...'))
ereFull$Proportion... <- as.numeric(ereFull$Proportion...)
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/DIRS'] <- 'TE/Class_I/DIRS'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Bel-Pao'] <- 'TE/Class_I/LTR/Bel-Pao'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'] <- 'TE/Class_I/LTR/Ty1_copia'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_I/Penelope'] <- 'TE/Class_I/Penelope'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Helitron'] <- 'TE/Class_II/Helitron'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_2/Maverick'] <- 'TE/Class_II/Maverick'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
ereFull$Annotation[ereFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
ereFull <- ereFull %>%
  mutate(Annotation = factor(Annotation, levels = c("TE/Class_I/DIRS",
                                                    "TE/Class_I/LINE",
                                                    "TE/Class_I/LTR",
                                                    "TE/Class_I/LTR/Bel-Pao",
                                                    "TE/Class_I/LTR/Ty1_copia",
                                                    "TE/Class_I/LTR/Ty3_gypsy",
                                                    "TE/Class_I/Penelope",
                                                    "TE/Class_II/Helitron",
                                                    "TE/Class_II/Maverick",
                                                    "rDNA",
                                                    "satDNA",
                                                    "Unidentified repeat",
                                                    "Non repetitive")))

ereFull <- ereFull %>%
  mutate(Indiv = fct_rev(factor(Indiv, levels = c("pol", "med", "epd", "tri",
                                                  "meo", "alb", "aep", "mne",
                                                  "gon", "gor", "plu", "neo",
                                                  "lef", "sti", "mon", "sci",
                                                  "pro", "mls", "sty", "cla",
                                                  "oem", "epp", "pha", "nip",
                                                  "aet", "mel", "his", "ron",
                                                  "cas", "nel", "tyn", "niv",
                                                  "cal", "cll", "gra", "ira",
                                                  "ott", "eri", "lig", "eur",
                                                  "man", "eps", "pan", "sth",
                                                  "mag", "dle", "dis"))))

## Plot genome proportions
erePlot <- ereFull %>%
  ggplot(aes(x = Indiv, y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(a) ", italic("Erebia")))) +
  scale_y_continuous(limits = c(0, 50.001), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4", "#b2df8a",
                               "#33a02c", "darkslategrey", "seagreen",
                               "#b15928", "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II",
                               "TE/Class_II/TIR", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10))
erePlot


#### Carex ####
## Load datasets
carexFiles_run1 <- list.files(path = "input/carex", pattern = "*_run1.txt",
                              full.names = TRUE)
for (f in 1:length(carexFiles_run1)) {
  df <- read.delim(carexFiles_run1[f])
  # Fuse All and All/repeat together
  df$Annotation[df$Annotation == 'All'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat'] <- 'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element'] <-
    'Unidentified repeat'
  df$Annotation[df$Annotation == 'All/repeat/mobile_element/Class_I'] <-
    'Unidentified repeat'
  # Fuse rDNA together
  df$Annotation[grep("All/repeat/rDNA/*", df$Annotation)] <- 'All/repeat/rDNA'
  # Fuse All/repeat/mobile_element/Class_II/Subclass_1/TIR together
  df$Annotation[grep("All/repeat/mobile_element/Class_II/Subclass_1/TIR*",
                     df$Annotation)] <-
    'All/repeat/mobile_element/Class_II/Subclass_1/TIR'
  # Fuse All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy together
  df$Annotation[grep("All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy*",
                     df$Annotation)] <-
    'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'
  # Fuse All/repeat/mobile_element/Class_I/LTR/Ty1_copia together
  df$Annotation[grep("All/repeat/mobile_element/Class_I/LTR/Ty1_copia*",
                     df$Annotation)] <-
    'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'
  df$Proportion... <- as.numeric(df$Proportion...)
  df$Number_of_reads <- as.numeric(df$Number_of_reads)
  df <- aggregate(cbind(Number_of_reads, Proportion...) ~ Annotation,
                  data = df, FUN = sum)
  # Fuse DNA transposons together (just for calculation of proportion)
  #df$Annotation[grep("All/repeat/mobile_element/Class_II/*", df$Annotation)] <- 'DNAtransposons'
  # Fuse retrotransposons together (just for calculation of proportion)
  #df$Annotation[grep("All/repeat/mobile_element/Class_I/*", df$Annotation)] <- 'retrotransposons'
  #df$Proportion... <- as.numeric(df$Proportion...)
  #df$Number_of_reads <- as.numeric(df$Number_of_reads)
  #df <- aggregate(cbind(Number_of_reads, Proportion...) ~ Annotation, data = df, FUN = sum)
  # Add a row with non-repeat proportion
  non_repeat_prop <- c("Non repetitive",
                       df[df$Annotation == 'Total', 'Number_of_reads'] -
                         sum(df[df$Annotation != 'Total', 'Number_of_reads']),
                       50 - sum(df[df$Annotation != 'Total', 'Proportion...']))
  df <- rbind(df, non_repeat_prop)
  # Add a column with individual name
  column_name <- paste0(str_match(carexFiles_run1[f],
                                  "Genome_proportion_(.*?)_run1.txt")[, 2])
  df$Indiv <- column_name
  # Set file name
  dataframe_name <- str_match(carexFiles_run1[f], "(.*?)_run1.txt")[, 2]
  assign(dataframe_name, df)
}

## Merge datasets
carexList <- lapply(ls(pattern = "input/carex/Genome_proportion*"), get)
carexFull <- carexList %>%
  reduce(full_join,
         by = c('Annotation', 'Indiv', 'Number_of_reads', 'Proportion...'))
carexFull$Number_of_reads <- as.numeric(carexFull$Number_of_reads)
carexFull$Proportion... <- as.numeric(carexFull$Proportion...)
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_I/LINE'] <- 'TE/Class_I/LINE'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR'] <- 'TE/Class_I/LTR'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'] <- 'TE/Class_I/LTR/Ty1_copia'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'] <- 'TE/Class_I/LTR/Ty3_gypsy'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_II'] <- 'TE/Class_II'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/mobile_element/Class_II/Subclass_1/TIR'] <- 'TE/Class_II/TIR'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/rDNA'] <- 'rDNA'
carexFull$Annotation[carexFull$Annotation == 'All/repeat/satellite'] <- 'satDNA'
carexFull <- carexFull %>%
  mutate(Annotation = factor(Annotation, levels=c("TE/Class_I/LINE",
                                                  "TE/Class_I/LTR",
                                                  "TE/Class_I/LTR/Ty1_copia",
                                                  "TE/Class_I/LTR/Ty3_gypsy",
                                                  "TE/Class_II",
                                                  "TE/Class_II/TIR",
                                                  "rDNA",
                                                  "satDNA",
                                                  "Unidentified repeat",
                                                  "Non repetitive")))

carexFull <- carexFull[(carexFull$Indiv == "ME1" | carexFull$Indiv == "ME3" |
                        carexFull$Indiv == "ME5" | carexFull$Indiv == "ME6" |
                        carexFull$Indiv == "ME7" | carexFull$Indiv == "ME10" |
                        carexFull$Indiv == "ME11" | carexFull$Indiv == "ME13" |
                        carexFull$Indiv == "ME15" | carexFull$Indiv == "ME18" |
                        carexFull$Indiv == "ME19" | carexFull$Indiv == "ME21" |
                        carexFull$Indiv == "ME23" | carexFull$Indiv == "ME25"),]
carexFull <- carexFull %>%
  mutate(Indiv = factor(Indiv, levels = c("ME3", "ME1", "ME19", "ME10", "ME15",
                                          "ME25", "ME5", "ME6", "ME21", "ME23",
                                          "ME7", "ME13", "ME18", "ME11")))

## Plot genome proportions
carexPlot_Prop <- carexFull %>% filter(Annotation != "Total") %>%
  ggplot(aes(x = fct_rev(Indiv), y = Proportion...)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Genome content (%)") +
  labs(title = expression(paste("(b) ", italic("Carex")))) +
  scale_y_continuous(limits = c(0, 50.002), expand = c(0.01,0)) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                               "#6a3d9a", "#a6cee3", "#1f78b4", "#b2df8a",
                               "#33a02c", "darkslategrey", "seagreen",
                               "#b15928", "#ffff99", "#5d5d5d", "#DCDCDC"),
                    limits = c("TE/Class_I/DIRS", "TE/Class_I/LINE",
                               "TE/Class_I/LTR", "TE/Class_I/LTR/Bel-Pao",
                               "TE/Class_I/LTR/Ty1_copia",
                               "TE/Class_I/LTR/Ty3_gypsy",
                               "TE/Class_I/Penelope", "TE/Class_II",
                               "TE/Class_II/TIR", "TE/Class_II/Helitron",
                               "TE/Class_II/Maverick", "rDNA", "satDNA",
                               "Unidentified repeat", "Non repetitive"),
                    name = "Annotation") +
  coord_flip() +
  theme(aspect.ratio = 1/2) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10))
carexPlot_Prop

## Plot genome composition in terms of number of reads
carexPlot_Nb <- carexFull %>% filter(Annotation != "Total") %>%
  ggplot(aes(x = fct_rev(Indiv), y = Number_of_reads)) +
  geom_col(aes(fill = fct_rev(Annotation)), width = 0.3) +
  xlab("") + ylab("Number of reads") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = c("#DCDCDC", "#5d5d5d", "#ffff99", "#b15928",
                               "seagreen", "darkslategrey", "#a6cee3",
                               "#6a3d9a", "#fdbf6f", "#e31a1c"),
                    guide = guide_legend(reverse = TRUE), name = "Annotation") +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 10), legend.position = "bottom")
pdf("output/Fig_Genome_content_carex.pdf", width = 6.65354, height = 3.3)
carexPlot_Nb
dev.off()

### 2 plots together
pdf("output/Fig_Genome_proportions_genus.pdf", width = 6.65354)
ggarrange(erePlot, carexPlot_Prop, ncol = 1, nrow = 2, align = "v",
          common.legend = TRUE, legend = "bottom", heights = c(2.2,1))
dev.off()
