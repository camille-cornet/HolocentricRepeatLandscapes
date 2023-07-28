rm(list = ls())

library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(janitor) # For the function adorn_totals
library(phylobase) # For the function phylo4d
library(phylosignal) # To test for a phylogenetic signal
library(vegan)
library(dendextend) # To compare 2 dendrograms
library(phylogram)
library(MCMCglmm) # for phylo glmm
library(ggtree)
library(ggpubr)
library(patchwork)

#### Erebia ####

## Load data
# Repeat dataset restricted to the annotated clusters
ere_repeats <- read.delim("input/erebia_clusters.txt")
ere_repeats$Annotation[ere_repeats$Annotation == 'All'] <- 'Unindentified repeat'
ere_repeats$Annotation[ere_repeats$Annotation == 'All/repeat'] <- 'Unindentified repeat'
ere_repeats$Annotation[ere_repeats$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
ere_repeats$Annotation[ere_repeats$Annotation == 'All/repeat/mobile_element/Class_I'] <- 'Unindentified repeat'

ere_repeats_full <- read.delim("input/erebia_clusters.txt")
ere_repeats_full <- subset(ere_repeats_full,
                           Annotation != "All/organelle/mitochondria")
ere_repeats_full <- t(ere_repeats_full)
ere_repeats_full <- as.data.frame(ere_repeats_full)
colnames(ere_repeats_full) <- ere_repeats_full[1,]
ere_repeats_full <- ere_repeats_full[-1,]
colnames(ere_repeats_full) <- make.names(colnames(ere_repeats_full))
for (i in 1:ncol(ere_repeats_full)) {
  ere_repeats_full[,i] <- as.numeric(ere_repeats_full[,i])
}
# Sum rows with same annotation
ere_repeats <- aggregate(. ~ Annotation, data = ere_repeats, FUN = sum)
# Add row with total amount of repeats (without mitochondria)
ere_repeats <-
  ere_repeats[ere_repeats$Annotation != "All/organelle/mitochondria",] %>%
  adorn_totals("row")
# Other variables
ere_variables <- read.delim("input/erebia_variables.txt")
# Use Chromosome 2N instead of N !
ere_variables$Chromosome_2N <- 2*ere_variables$Chromosome_N
# Merging both datasets
ere_repeats <- t(ere_repeats)
ere_repeats <- as.data.frame(ere_repeats)
colnames(ere_repeats) <- ere_repeats[1,]
ere_repeats <- ere_repeats[-1,]
colnames(ere_repeats) <- make.names(colnames(ere_repeats))
for (i in 1:ncol(ere_repeats)) {
  ere_repeats[,i] <- as.numeric(ere_repeats[,i])
}
ere <- cbind(ere_variables,ere_repeats)
rownames(ere) <- ere[,3]
rownames(ere_repeats_full) <- ere[,3]
write.table(ere, file = "output/erebia_dataset.txt", sep = "\t")
# Nexus file
ere_phylo_full <- read.nexus("input/erebia_nondated_phylo.nex")

## Phylogenetic regression
# Check if same species in phylo and in dataset
ere_check <- name.check(ere_phylo_full, ere)
ere_check
# Drop species in phylo that are not in dataset
ere_phylo <- drop.tip(ere_phylo_full, ere_check$tree_not_data)
# Define a variance-covariance structure
ere_pagel <- corPagel(1, ere_phylo)
# Phylogenetic Generalized Least Squares
model_DIRS <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.DIRS,
                  data = ere, correlation = ere_pagel)
summary(model_DIRS) # p 0.0965
model_LINE <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.LINE,
                  data = ere, correlation = ere_pagel)
summary(model_LINE) # p 0.3410
model_LTR <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.LTR,
                 data = ere, correlation = ere_pagel)
summary(model_LTR) # p 0.0372
model_BelPao <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.LTR.Bel.Pao,
                    data = ere, correlation = ere_pagel)
summary(model_BelPao) # p 0.4367
model_Ty1Copia <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.LTR.Ty1_copia,
                      data = ere, correlation = ere_pagel)
summary(model_Ty1Copia) # p 0.1540
model_Ty3gypsy <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.LTR.Ty3_gypsy,
                      data = ere, correlation = ere_pagel)
summary(model_Ty3gypsy) # p 0.1284
model_Penelope <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_I.Penelope,
                      data = ere, correlation = ere_pagel)
summary(model_Penelope) # p 0.4087
model_Helitron <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_II.Subclass_2.Helitron,
                      data = ere, correlation = ere_pagel)
summary(model_Helitron) # p 0.0108
model_Maverick <- gls(Chromosome_2N ~ All.repeat.mobile_element.Class_II.Subclass_2.Maverick,
                      data = ere, correlation = ere_pagel)
summary(model_Maverick) # p 0.0409
model_rDNA <- gls(Chromosome_2N ~ All.repeat.rDNA.45S_rDNA,
                  data = ere, correlation = ere_pagel)
summary(model_rDNA) # p 0.0978
model_satellite <- gls(Chromosome_2N ~ All.repeat.satellite,
                       data = ere, correlation = ere_pagel)
summary(model_satellite) # p 0.2351

## Scatterplot of repeats against chromosome number
ere_DIRS_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.DIRS)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of DIRS") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_LINE_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.LINE)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of LINE") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_LTR_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.LTR)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of LTR") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_BelPao_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.LTR.Bel.Pao)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Bel-Pao") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_Ty1Copia_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.LTR.Ty1_copia)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Ty1/Copia") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_Ty3Gypsy_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.LTR.Ty3_gypsy)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Ty3/Gypsy") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_Penelope_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_I.Penelope)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Penelope") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_Helitron_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_II.Subclass_2.Helitron)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Helitron") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_Maverick_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.mobile_element.Class_II.Subclass_2.Maverick)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Maverick") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_rDNA_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.rDNA.45S_rDNA)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of rDNA") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
ere_satDNA_plot <- ggplot(ere, aes(x = Chromosome_2N, y = All.repeat.satellite)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of satDNA") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))

pdf("output/erebia_PGLS.pdf", width = 7)
ere_scatterplot <- egg::ggarrange(ere_DIRS_plot, ere_LINE_plot, ere_LTR_plot, ere_BelPao_plot,
               ere_Ty1Copia_plot, ere_Ty3Gypsy_plot, ere_Penelope_plot,
               ere_Helitron_plot, ere_Maverick_plot, ere_rDNA_plot,
               ere_satDNA_plot,
               label.args = list(gp = grid::gpar(fontface = "plain")),
               ncol = 4, nrow = 3)
dev.off()

# Correction for multiple testing
ere_p <- c(0.0965, 0.3410, 0.0372, 0.4367, 0.1540, 0.1284,
           0.4087, 0.0108, 0.0409, 0.0978,0.2351)
ere_p_adjust <- p.adjust(ere_p, method = "fdr")
ere_p_adjust

## Testing for a phylogenetic signal
ere_nexus <- readNexus("input/erebia_nondated_phylo.nex")
ere_nexus <- phylobase::subset(ere_nexus,
                               tips.exclude = ere_check$tree_not_data)
ere_data <- cbind(ere[,4:16])
ere_phylo4d <- phylo4d(ere_nexus, ere_data)
ere_signal <- phyloSignal(ere_phylo4d, methods = c("all"),
                          reps = 9999, W = NULL)
write.table(ere_signal, file = "output/Erebia_phylo_signal.txt", sep = "\t")
## Correction for multiple testing
ere_signal_p <- c(1, 1, 1, 0.008, 0, 0.003, 0, 0.003, 0.617, 0.217, 1)
ere_signal_p_adj <- p.adjust(ere_signal_p, method = "fdr")
ere_signal_p_adj

## Plot phylogeny in front of clustering according to repeats
ere_phylo_dendro <- as.dendrogram.phylo(ere_phylo)
ere_phylo_reroot <- ape::root(ere_phylo, node = 92, resolve.root = T)
ere_plot <- ggtree(ere_phylo_reroot, branch.length = "none") +
  geom_tiplab(fontface = 3) +
  xlim(-100, 100)
pdf("output/erebia_phylo.pdf", width = 6.65354)
flip(ere_plot, 80, 67)
dev.off()

ere_phylo_dendro <- as.dendrogram.phylo(ere_phylo_reroot, method = "extend")
ere_dendro <- as.dendrogram(hclust(vegdist(ere_repeats_full, method = "bray")))

pdf("output/erebia_tanglegram.pdf", width = 6.65354)
ere_tangle <- dendlist(ere_phylo_dendro, ere_dendro) %>%
  dendextend::untangle(method = "step2side") %>%
  set("labels_cex", 1.4) %>%
  tanglegram(margin_inner = 7, highlight_branches_lwd = F, common_subtrees_color_lines = F,
             lwd = 2, axes = F, rank_branches = T, highlight_distinct_edges = F,
             color_lines = rev(c("grey", "grey", "#a00e00", "#a00e00", "#a00e00",
                             "#a00e00", "#a00e00", "#d04e00", "#d04e00",
                             "#d04e00", "#d04e00", "#d04e00", "#d04e00",
                             "#d04e00", "#d04e00", "#d04e00", "#d04e00",
                             "#092a6e", "#092a6e", "grey", "grey", "#f6c200",
                             "#f6c200", "#f6c200", "#f6c200", "#f6c200",
                             "#f6c200", "#0086a8", "#0086a8", "#0086a8",
                             "#0086a8", "#0086a8", "grey", "#6600ff", "#6600ff",
                             "#6600ff", "#6600ff", "#6600ff", "#6600ff",
                             "#6600ff", "#6600ff", "#6600ff", "#6600ff",
                             "#6600ff", "#dfefb2", "#dfefb2", "#a1c198")),
             main_left = expression(paste("(a) ", italic("Erebia"))))
dev.off()


#### Carex ####
## Load data
# Repeat dataset restricted to the annotated clusters
carex_repeats <- read.delim("input/carex_clusters.txt")

carex_repeats_full <- read.delim("input/carex_clusters.txt")
carex_repeats_full <- subset(carex_repeats_full,
                             Annotation != "All/organelle/mitochondria")
carex_repeats_full <- subset(carex_repeats_full,
                             Annotation != "All/organelle/plastid")
carex_repeats_full <- t(carex_repeats)
carex_repeats_full <- as.data.frame(carex_repeats_full)
colnames(carex_repeats_full) <- carex_repeats_full[1,]
carex_repeats_full <- carex_repeats_full[-1,]
colnames(carex_repeats_full) <- make.names(colnames(carex_repeats_full))
for (i in 1:ncol(carex_repeats_full)) {
  carex_repeats_full[,i] <- as.numeric(carex_repeats_full[,i])
}

# Add row with total amount of repeats (without mitochondria)
carex_repeats <-
  carex_repeats[carex_repeats$Annotation != "All/organelle/mitochondria" &
                  carex_repeats$Annotation != "All/organelle/plastid",] %>%
  adorn_totals("row")

# Fuse repeats together
carex_repeats$Annotation[carex_repeats$Annotation == 'All' | carex_repeats$Annotation == 'All/repeat' | carex_repeats$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
carex_repeats$Annotation[grep("All/repeat/mobile_element/Class_II/Subclass_1/TIR/*", carex_repeats$Annotation)] <- 'All/repeat/mobile_element/Class_II/Subclass_1/TIR'
carex_repeats$Annotation[grep("All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy/*", carex_repeats$Annotation)] <- 'All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy'
carex_repeats$Annotation[grep("All/repeat/mobile_element/Class_I/LTR/Ty1_copia/*", carex_repeats$Annotation)] <- 'All/repeat/mobile_element/Class_I/LTR/Ty1_copia'
carex_repeats$Annotation[grep("All/repeat/rDNA/*", carex_repeats$Annotation)] <- 'All/repeat/rDNA'

# Sum rows with same annotation
carex_repeats <- aggregate(. ~ Annotation, data = carex_repeats, FUN = sum)

# Other variables
carex_variables <- read.delim("input/carex_variables.txt")
# Merging both datasets
carex_repeats <- t(carex_repeats)
carex_repeats <- as.data.frame(carex_repeats)
colnames(carex_repeats) <- carex_repeats[1,]
carex_repeats <- carex_repeats[-1,]
colnames(carex_repeats) <- make.names(colnames(carex_repeats))
for (i in 1:ncol(carex_repeats)) {
  carex_repeats[,i] <- as.numeric(carex_repeats[,i])
}
rownames(carex_variables) <- carex_variables[,1]
carex_variables <- carex_variables[,-1]
carex <- merge(carex_variables, carex_repeats, by = 'row.names', all = TRUE)
rownames(carex) <- carex[,2]
rownames(carex_repeats_full) <- carex[,2]
carex$rate_of_chromevol <- as.factor(carex$rate_of_chromevol)
write.table(carex, file = "output/carex_dataset.txt", sep = "\t")
# Nexus file
carex_phylo <- read.tree("input/carex_nondated_phylo.nex")

## Phylogenetic regression
# Check is same species in phylo and in dataset
carex_check <- name.check(carex_phylo, carex)
carex_check
# Define a variance-covariance structure based on the model of Brownian evolution
# Or should I use Ornstein-Uhlenbeck?
carex_bm <- corBrownian(1, carex_phylo)
carex_pagel <- corPagel(1, carex_phylo)
# Also try with other models!
# Try a phylo regression
model_LINE <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_I.LINE,
                  data = carex, correlation = carex_pagel)
summary(model_LINE) # p 0.1532
model_LTR <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_I.LTR,
                 data = carex, correlation = carex_pagel)
summary(model_LTR) # p 0.0238
model_Ty1Copia <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_I.LTR.Ty1_copia,
                      data = carex, correlation = carex_pagel)
summary(model_Ty1Copia) # p 0.0663
model_Ty3gypsy <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_I.LTR.Ty3_gypsy,
                      data = carex, correlation = carex_pagel)
summary(model_Ty3gypsy) # p 0.0237
model_TIR <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_II.Subclass_1.TIR,
                 data = carex, correlation = carex_pagel)
summary(model_TIR) # p 0.3614
model_Helitron <- gls(Chromosome_2N_mean ~ All.repeat.mobile_element.Class_II.Subclass_2.Helitron,
                      data = carex, correlation = carex_pagel)
summary(model_Helitron) # p 0.6811
model_rDNA <- gls(Chromosome_2N_mean ~ All.repeat.rDNA,
                  data = carex, correlation = carex_pagel)
summary(model_rDNA) # p 0.1467
model_satellite <- gls(Chromosome_2N_mean ~ All.repeat.satellite,
                       data = carex, correlation = carex_pagel)
summary(model_satellite) # p 0.4999

## Scatterplot of repeats against chromosome number
carex_LINE_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_I.LINE)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of LINE") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_LTR_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_I.LTR)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of LTR") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_Ty1Copia_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_I.LTR.Ty1_copia)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Ty1/Copia") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_Ty3Gypsy_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_I.LTR.Ty3_gypsy)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Ty3/Gypsy") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_TIR_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_II.Subclass_1.TIR)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of TIR") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_Helitron_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.mobile_element.Class_II.Subclass_2.Helitron)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of Helitron") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_rDNA_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.rDNA)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of rDNA") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
carex_satDNA_plot <- ggplot(carex, aes(x = Chromosome_2N_mean, y = All.repeat.satellite)) +
  geom_point(size = 0.6) +
  xlab("Chromosome number (2N)") +
  ylab("Abundance of satDNA") +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))

pdf("output/carex_PGLS.pdf", width = 7.5)
carex_scatterplot <- egg::ggarrange(carex_LINE_plot, carex_LTR_plot,
                                  carex_Ty1Copia_plot, carex_Ty3Gypsy_plot, carex_TIR_plot,
                                  carex_Helitron_plot, carex_rDNA_plot, carex_satDNA_plot,
                                  label.args = list(gp = grid::gpar(fontface = "plain")),
                                  ncol = 4, nrow = 2)
dev.off()

# Correction for multiple testing
p_carex <- c(0.1532, 0.0238, 0.0663, 0.0237, 0.3614, 0.6811, 0.1467, 0.4999)
p_carex_adjust <- p.adjust(p_carex, method = "fdr")
p_carex_adjust

# Other variables in the model
model_LINE <- gls(All.repeat.mobile_element.Class_I.LINE ~
                    Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                  data = carex, correlation = carex_pagel)
summary(model_LINE)
model_LTR <- gls(All.repeat.mobile_element.Class_I.LTR ~
                   Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                 data = carex, correlation = carex_pagel)
summary(model_LTR)
model_Ty1Copia <- gls(All.repeat.mobile_element.Class_I.LTR.Ty1_copia ~
                        Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                      data = carex, correlation = carex_pagel)
summary(model_Ty1Copia)
model_Ty3gypsy <- gls(All.repeat.mobile_element.Class_I.LTR.Ty3_gypsy ~
                        Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                      data = carex, correlation = carex_pagel)
summary(model_Ty3gypsy)
model_TIR <- gls(All.repeat.mobile_element.Class_II.Subclass_1.TIR ~
                   Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                 data = carex, correlation = carex_pagel)
summary(model_TIR)
model_Helitron <- gls(All.repeat.mobile_element.Class_II.Subclass_2.Helitron ~
                        Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                      data = carex, correlation = carex_pagel)
summary(model_Helitron)
model_rDNA <- gls(All.repeat.rDNA ~
                    Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                  data = carex, correlation = carex_pagel)
summary(model_rDNA)
model_satellite <- gls(All.repeat.satellite ~
                         Chromosome_2N_mean + Nb_karyotypes + Proba_high,
                       data = carex, correlation = carex_pagel)
summary(model_satellite)

# Correction for multiple testing
p_carex_2N <- c(0.260, 0.031, 0.272, 0.029, 0.280, 0.969, 0.307, 0.394)
p_carex_adjust_2N <- p.adjust(p_carex_2N, method = "fdr")
p_carex_adjust_2N

p_carex_Nb <- c(0.182, 0.074, 0.618, 0.079, 0.314, 0.658, 0.473, 0.252)
p_carex_adjust_Nb <- p.adjust(p_carex_Nb, method = "fdr")
p_carex_adjust_Nb

p_carex_chromevo <- c(0.0002, 0.014, 0.043, 0.028, 0.248, 0.005, 0.998, 0.174)
p_carex_adjust_chromevo <- p.adjust(p_carex_chromevo, method = "fdr")
p_carex_adjust_chromevo

## Association between repeats and genome size
model_genomesize <- gls(Genome_size ~ Total, data = carex,
                        correlation = carex_pagel)
summary(model_genomesize)

## Testing for a phylogenetic signal
carex_nexus <- phylobase::readNewick("input/carex_nondated_phylo.nex")
carex_data <- cbind(carex[,7:17])[,-1]
carex_phylo4d <- phylo4d(carex_nexus, carex_data)
carex_signal <- phyloSignal(carex_phylo4d, methods = c("all"),
                            reps = 9999, W = NULL)
carex_signal
write.table(carex_signal, file = "output/Carex_phylo_signal.txt", sep = "\t")
## Correction for multiple testing
carex_signal_p <- c(0, 0.003, 0.058, 0.007, 0.026, 0.134, 1, 0.384)
carex_signal_p_adj <- p.adjust(carex_signal_p, method = "fdr")
carex_signal_p_adj

## Plot phylogeny in front of clustering according to repeats
carex_repeats_full <- carex_repeats_full[order(row.names(carex_repeats_full)), ]
carex_variables <- carex_variables[order(row.names(carex_variables)), ]
rownames(carex_repeats_full) <- carex_variables[,1]

carex_phylo$tip.label <- c('nigra', 'magellanica', 'extensa', 'laevigata',
                           'helodes', 'lepidocarpa', 'sempervirens',
                           'caryophyllea', 'lucennoiberica', 'furva',
                           'echinata', 'leporina', 'maritima', 'capitata')
pdf("output/carex_phylo.pdf", width = 6.65354)
ggtree(carex_phylo, branch.length = "none") +
  geom_tiplab(fontface = 3) +
  xlim(-50, 50) + ylim(-16,16)
dev.off()

carex_phylo_dendro <- as.dendrogram.phylo(carex_phylo)
carex_dendro <- as.dendrogram(hclust(vegdist(carex_repeats_full,
                                             method = "bray")))
labels(carex_dendro) <- c('sempervirens', 'furva', 'lucennoiberica', 'echinata',
                          'maritima', 'leporina', 'laevigata', 'helodes',
                          'lepidocarpa', 'magellanica', 'nigra', 'caryophyllea',
                          'capitata', 'extensa')

pdf("output/carex_tanglegram.pdf", width = 6.65354, height = 3.3)
carex_tangle <- dendlist(carex_phylo_dendro, carex_dendro) %>%
  dendextend::untangle(method = "step2side") %>%
  set("labels_cex", 1.4) %>%
  tanglegram(margin_inner = 9, highlight_branches_lwd = F,
             common_subtrees_color_lines = F,
             lwd = 2, axes = F, rank_branches = T,
             highlight_distinct_edges = F,
             main_left = expression(paste("(b) ", italic("Carex"))))
dev.off()
