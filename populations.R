rm(list = ls())

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ape) # fct PCoA
library(egg) # to combine plots

cbPalette <- c("#E69F00", "#0072B2", "#009E73", "#56B4E9",
               "#D55E00", "#CC79A7", "grey33")

#### PCoA Nivalis ####
# Load dataset
niv_clusters <- read.delim("input/niv_clusters.txt")
niv_clusters <- subset(niv_clusters, select = -c(SAJ_3)) # Drop sibling
niv_clusters <- subset(niv_clusters, Annotation != "All/organelle/mitochondria")
niv_clusters <- t(niv_clusters)
niv_clusters <- as.data.frame(niv_clusters)
niv_clusters <- niv_clusters[-11,] # Remove Annotation row
for (i in 1:ncol(niv_clusters)) {
  niv_clusters[,i] <- as.numeric(niv_clusters[,i])
}

# PCoA
dist_niv_clusters <- vegdist(niv_clusters,  method = "bray")
niv_pcoa_clusters <- pcoa(dist_niv_clusters)

# Plot
niv_coord_pcoa <- as.data.frame(niv_pcoa_clusters$vectors[,1:2])
niv_coord_pcoa$pop <- c('GRI', 'GRI', 'GRI', 'GRO', 'GRO',
                        'GRO', 'SAJ', 'SAJ', 'SCH', 'SCH')
## To get square grid
break_niv <- function(limits) {
  seq(floor(limits[1]), ceiling(limits[2]), 0.025)
}
niv_coord_pcoa_plot <- ggplot(niv_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(niv_pcoa_clusters$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(niv_pcoa_clusters$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = rownames(niv_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.05, color = "darkgrey") +
  geom_hline(yintercept = 0, linewidth = 0.05, color = "darkgrey") +
  coord_fixed() +
  scale_x_continuous(breaks = break_niv, limits = c(-0.08, 0.07)) +
  scale_y_continuous(breaks = break_niv, limits = c(-0.04, 0.045)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(c) ", italic("Erebia nivalis")))) +
  theme_bw() +
  theme(text = element_text(size = 8), legend.position = "bottom",
        legend.justification = 'left', panel.grid = element_blank(),
        plot.title = element_text(margin = margin(5,0,0,0)),
        plot.margin = grid::unit(c(0,5,0,0), "mm"),
        legend.margin = margin(0,-19,0,0),
        axis.title.y = element_text(margin = margin(0,-0.5,0,2))) +
  theme(panel.border = element_rect(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        legend.box.spacing = grid::unit(c(0,0), "mm"),
        legend.key.width = grid::unit(0.5, "mm")) +
  guides(colour = guide_legend(nrow = 1))

#### PCoA Tyndarus ####
# Load dataset
tyn_clusters <- read.delim("input/tyn_clusters.txt")
tyn_clusters <- subset(tyn_clusters, select = -c(GRW_3)) # Drop sibling
tyn_clusters <- subset(tyn_clusters, Annotation != "All/organelle/mitochondria")
tyn_clusters <- t(tyn_clusters)
tyn_clusters <- as.data.frame(tyn_clusters)
tyn_clusters <- tyn_clusters[-21,] # Remove Annotation row
for (i in 1:ncol(tyn_clusters)) {
  tyn_clusters[,i] <- as.numeric(tyn_clusters[,i])
}

# PCoA
dist_tyn_clusters <- vegdist(tyn_clusters,  method = "bray")
tyn_pcoa_clusters <- pcoa(dist_tyn_clusters)

# Plot
tyn_coord_pcoa <- as.data.frame(tyn_pcoa_clusters$vectors[,1:2])
tyn_coord_pcoa$pop <- c('ARO', 'ARO', 'ARO', 'COM', 'COM', 'COM',
                        'GOT', 'GOT', 'GOT', 'GRF', 'GRF', 'GRF',
                        'GRW', 'GRW', 'LAU', 'LAU', 'LAU',
                        'SAN', 'SAN', 'SAN')
## To get square grid
break_tyn <- function(limits) {
  seq(floor(limits[1]), ceiling(limits[2]), 0.025)
}
tyn_coord_pcoa_plot <- ggplot(tyn_coord_pcoa,
                              aes(x=Axis.1, y=Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(tyn_pcoa_clusters$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(tyn_pcoa_clusters$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = rownames(tyn_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 5, force = 0.01,
                  max.time = 10, box.padding = 0.1) +
  geom_vline(xintercept = 0, linewidth = 0.05, color = "darkgrey") +
  geom_hline(yintercept = 0, linewidth = 0.05, color = "darkgrey") +
  coord_fixed() +
  scale_x_continuous(breaks = break_tyn, limits = c(-0.075, 0.075)) +
  scale_y_continuous(breaks = break_tyn, limits = c(-0.055, 0.03)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(b) ", italic("Erebia tyndarus")))) +
  theme_bw() +
  theme(text = element_text(size = 8), legend.position = "bottom",
        legend.justification = 'left', panel.grid = element_blank(),
        plot.title = element_text(margin = margin(5,0,0,0)),
        plot.margin = grid::unit(c(0,5,0,0), "mm"),
        legend.margin = margin(0,-19,0,0),
        axis.title.y = element_text(margin = margin(0,-0.5,0,2))) +
  theme(panel.border = element_rect(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        legend.box.spacing = grid::unit(c(0,0), "mm"),
        legend.key.width = grid::unit(0.5, "mm")) +
  guides(colour = guide_legend(nrow = 1))

#### PCoA Cassioides ####
# Load dataset
cas_clusters <- read.delim("input/cas_clusters.txt")
cas_clusters <- subset(cas_clusters, Annotation != "All/organelle/mitochondria")
cas_clusters <- t(cas_clusters)
cas_clusters <- as.data.frame(cas_clusters)
cas_clusters <- cas_clusters[-22,] # Remove Annotation row
for (i in 1:ncol(cas_clusters)) {
  cas_clusters[,i] <- as.numeric(cas_clusters[,i])
}

# PCoA
dist_cas_clusters <- vegdist(cas_clusters, method = "bray")
cas_pcoa_clusters <- pcoa(dist_cas_clusters)

# Plot
cas_coord_pcoa <- as.data.frame(cas_pcoa_clusters$vectors[,1:2])
cas_coord_pcoa$pop <- c('GRF', 'GRF', 'GRF', 'GRW', 'GRW', 'GRW',
                        'KAN', 'KAN', 'KAN', 'PDM', 'PDM', 'PDM',
                        'ROU', 'ROU', 'ROU', 'SCH', 'SCH', 'SCH',
                        'STO', 'STO', 'STO')
## To get square grid
break_cas <- function(limits) {
  seq(floor(limits[1]), ceiling(limits[2]), 0.025)
}
cas_coord_pcoa_plot <- ggplot(cas_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(cas_pcoa_clusters$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(cas_pcoa_clusters$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = rownames(cas_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  coord_fixed() +
  scale_x_continuous(breaks = break_cas, limits = c(-0.07, 0.08)) +
  scale_y_continuous(breaks = break_cas, limits = c(-0.05, 0.035)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(a) ", italic("Erebia cassioides")))) +
  geom_vline(xintercept = 0, linewidth = 0.05, color = "darkgrey") +
  geom_hline(yintercept = 0, linewidth = 0.05, color = "darkgrey") +
  theme_bw() +
  theme(text = element_text(size = 8), legend.position = "bottom",
        legend.justification = 'left', panel.grid = element_blank(),
        plot.title = element_text(margin = margin(5,0,0,0)),
        plot.margin = grid::unit(c(0,5,0,0), "mm"),
        legend.margin = margin(0,-19,0,0),
        axis.title.y = element_text(margin = margin(0,-0.5,0,2))) +
  theme(panel.border = element_rect(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        legend.box.spacing = grid::unit(c(0,0), "mm"),
        legend.key.width = grid::unit(0.5, "mm")) +
  guides(colour = guide_legend(nrow = 1))

#### PCoA Pronoe ####
# Load dataset
pro_clusters <- read.delim("input/pro_clusters.txt")
pro_clusters <- subset(pro_clusters, Annotation != "All/organelle/mitochondria")
pro_clusters <- t(pro_clusters)
pro_clusters <- as.data.frame(pro_clusters)
pro_clusters <- pro_clusters[-12,] # Remove Annotation row
for (i in 1:ncol(pro_clusters)) {
  pro_clusters[,i] <- as.numeric(pro_clusters[,i])
}

# PCoA
dist_pro_clusters <- vegdist(pro_clusters,  method = "bray")
pro_pcoa_clusters <- pcoa(dist_pro_clusters)

# Plot
pro_coord_pcoa <- as.data.frame(pro_pcoa_clusters$vectors[,1:2])
pro_coord_pcoa$pop <- c('psathura', 'psathura', 'psathura','psathura',
                        'psathura', 'vergy', 'vergy', 'vergy', 'vergy',
                        'vergy', 'vergy')
## To get square grid
break_pro <- function(limits) {
  seq(floor(limits[1]), ceiling(limits[2]), 0.05)
}
pro_coord_pcoa_plot <- ggplot(pro_coord_pcoa,
                              aes(x = Axis.1, y  =Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(pro_pcoa_clusters$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(pro_pcoa_clusters$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Subspecies") +
  geom_text_repel(segment.size = 0.1, label = rownames(pro_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.05, color = "darkgrey") +
  geom_hline(yintercept = 0, linewidth = 0.05, color = "darkgrey") +
  coord_fixed() +
  scale_x_continuous(breaks = break_pro, limits = c(-0.1, 0.2)) +
  scale_y_continuous(breaks = break_pro, limits = c(-0.10, 0.07)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(d) ", italic("Erebia pronoe")))) +
  theme_bw() +
  theme(text = element_text(size = 8), legend.position = "bottom",
        legend.justification = 'left', panel.grid = element_blank(),
        plot.title = element_text(margin = margin(5,0,0,0)),
        plot.margin = grid::unit(c(0,5,0,0), "mm"),
        legend.margin = margin(0,-19,0,0),
        axis.title.y = element_text(margin = margin(0,-0.5,0,2))) +
  theme(panel.border = element_rect(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        legend.box.spacing = grid::unit(c(0,0), "mm"),
        legend.key.width = grid::unit(0.5, "mm")) +
  guides(colour = guide_legend(nrow = 1))

#### Final plot PCoA ####
pdf("output/Fig_PCoA.pdf", width = 6.65354)
egg::ggarrange(cas_coord_pcoa_plot, tyn_coord_pcoa_plot,
               niv_coord_pcoa_plot, pro_coord_pcoa_plot,
               label.args = list(gp = grid::gpar(fontface = "plain")),
               ncol = 2, nrow = 2)
dev.off()
