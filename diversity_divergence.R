rm(list = ls())

library(tidyverse)
library(ggrepel)
library(reshape2) # for function melt
library(vegan)
library(nlme)
library(egg) # to combine plots
library(scales) # to set nb of decimals
library(adespatial) # for function dbmem
library(ape) # for function pcoa

#### Cassioides ####
## Load repeat dataset summed by annotation
cas_clusters <- read.delim("input/cas_clusters.txt")
cas_clusters$Annotation[cas_clusters$Annotation == 'All'] <- 'Unindentified repeat'
cas_clusters$Annotation[cas_clusters$Annotation == 'All/repeat'] <- 'Unindentified repeat'
cas_clusters$Annotation[cas_clusters$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
cas_clusters$Annotation[cas_clusters$Annotation == 'All/repeat/mobile_element/Class_I'] <- 'Unindentified repeat'
cas_clusters$Annotation[grep("All/repeat/rDNA/*", cas_clusters$Annotation)] <- 'All/repeat/rDNA'
cas_clusters <- aggregate(. ~ Annotation, data = cas_clusters, FUN = sum)
cas_clusters <- t(cas_clusters)
colnames(cas_clusters) <- cas_clusters[1,]
cas_clusters <- cas_clusters[-1,]
cas_clusters <- cas_clusters[,-13] # Remove unidentified repeats
cas_clusters <- cas_clusters[,-1] # Remove organelle clusters
cas_clusters <- as.data.frame(cas_clusters)
for (i in 1:ncol(cas_clusters)) {
  cas_clusters[,i] <- as.numeric(cas_clusters[,i])
}
cas_clusters <- cas_clusters[order(row.names(cas_clusters)),]

# Calculate diversity indices
cas_simpson <- diversity(cas_clusters, index = "simpson")
cas_div <- as.data.frame(cas_simpson)
cas_div$indiv <- rownames(cas_div)

# Link with genetic diversity theta, based on individuals (with mlrho)
cas_theta <- read.delim("input/cas_theta.txt")
cas <- merge(cas_div, cas_theta, by = "indiv")
cas$pop <- c('grf', 'grf', 'grf', 'grw', 'grw', 'grw', 'kan', 'kan', 'kan',
             'pdm', 'pdm', 'pdm', 'rou', 'rou', 'rou', 'sch', 'sch', 'sch',
             'sto', 'sto', 'sto')
cas$pop <- as.factor(cas$pop)

# correlation between simpson diversity and theta
cas_cor <- cor.test(cas$cas_simpson, cas$theta_avg, method = "spearman")
cas_cor

# Load genetic distance dataset
cas_eucl_dist <- read.delim("input/cas_eucl_dist.txt", row.names = 1)
cas_eucl_dist <- as.matrix(cas_eucl_dist)
cas_eucl_dist <- cas_eucl_dist[order(rownames(cas_eucl_dist)),
                               order(colnames(cas_eucl_dist))]

# Do a redundancy analysis to compare the distance in repeats and SNP
# First, make a pcoa based on the euclidean distance matrix
cas_pcoa <- pcoa(cas_eucl_dist)
cas_X <- cas_pcoa$vectors
# Make distance matrix based on repeats
# First, load full dataset
cas_clusters <- read.delim("input/cas_clusters.txt")
cas_clusters <- subset(cas_clusters, Annotation != "All/organelle/mitochondria")
cas_clusters <- cas_clusters[,-22]
cas_clusters <- t(cas_clusters)
for (i in 1:ncol(cas_clusters)) {
  cas_clusters[,i] <- as.numeric(cas_clusters[,i])
}
cas_rep_dist <- vegdist(cas_clusters, method = "bray")
cas_rep_dist <- as.matrix(cas_rep_dist)
cas_rep_dist <- cas_rep_dist[order(rownames(cas_rep_dist)),
                             order(colnames(cas_rep_dist))]
cas_rep_dist <- as.dist(cas_rep_dist)
# Make a Moran Eigenvector's Maps based on the repeat distance matrix
cas_dbmem <- dbmem(cas_rep_dist)
# Run RDA
cas_rda <- rda(cas_X ~ ., data = cas_dbmem)
# Get significance and Rsquared
cas_Rsquared <- RsquareAdj(cas_rda)
cas_Rsquared
cas_signif <- anova.cca(cas_rda, permutations = 9999)
cas_signif

# Permanova
cas_permanova <- adonis2(cas_rep_dist ~ cas$pop,
                         data = cas, permutations = 9999)
cas_permanova

# Prepare for plotting
cas_rep_dist <- as.matrix(cas_rep_dist)
cas_rep_dist[upper.tri(cas_rep_dist, diag = TRUE)] <- NA
cas_rep_dist <- melt(cas_rep_dist, value.name = "rep_dist", na.rm = TRUE)
cas_eucl_dist <- as.matrix(cas_eucl_dist)
cas_eucl_dist[upper.tri(cas_eucl_dist, diag = TRUE)] <- NA
cas_eucl_dist <- melt(cas_eucl_dist, value.name = "eucl_dist", na.rm = TRUE)
cas_dist <- merge(cas_rep_dist, cas_eucl_dist, by = c("Var1","Var2"))


#### Tyndarus ####
## Load repeat dataset summed by annotation
tyn_clusters <- read.delim("input/tyn_clusters.txt")
tyn_clusters <- subset(tyn_clusters, select = -c(GRW_3)) # Remove sibling
tyn_clusters$Annotation[tyn_clusters$Annotation == 'All'] <- 'Unindentified repeat'
tyn_clusters$Annotation[tyn_clusters$Annotation == 'All/repeat'] <- 'Unindentified repeat'
tyn_clusters$Annotation[tyn_clusters$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
tyn_clusters$Annotation[tyn_clusters$Annotation == 'All/repeat/mobile_element/Class_I'] <- 'Unindentified repeat'
tyn_clusters$Annotation[grep("All/repeat/rDNA/*", tyn_clusters$Annotation)] <- 'All/repeat/rDNA'
tyn_clusters <- aggregate(. ~ Annotation, data = tyn_clusters, FUN = sum)
tyn_clusters <- t(tyn_clusters)
colnames(tyn_clusters) <- tyn_clusters[1,]
tyn_clusters <- tyn_clusters[-1,]
tyn_clusters <- tyn_clusters[,-12] # Remove unidentified repeats
tyn_clusters <- tyn_clusters[,-1] # Remove organelle clusters
tyn_clusters <- as.data.frame(tyn_clusters)
for (i in 1:ncol(tyn_clusters)) {
  tyn_clusters[,i] <- as.numeric(tyn_clusters[,i])
}
tyn_clusters <- tyn_clusters[order(row.names(tyn_clusters)),]

# Calculate diversity indices
tyn_simpson <- diversity(tyn_clusters, index = "simpson")
tyn_div <- as.data.frame(tyn_simpson)
tyn_div$indiv <- rownames(tyn_div)

# Link with genetic diversity theta, based on individuals (with mlrho)
tyn_theta <- read.delim("input/tyn_theta.txt")
tyn <- merge(tyn_div, tyn_theta, by = "indiv")
tyn$pop <- c('aro', 'aro', 'aro', 'com', 'com', 'com', 'got', 'got', 'got',
             'grf', 'grf', 'grf', 'grw', 'grw', 'lau', 'lau', 'lau',
             'san', 'san', 'san')
tyn$pop <- as.factor(tyn$pop)

# correlation between simpson diversity and theta
tyn_cor <- cor.test(tyn$tyn_simpson, tyn$theta_avg, method = "spearman")
tyn_cor

# Load genetic distance dataset
tyn_eucl_dist <- read.delim("input/tyn_eucl_dist.txt", row.names = 1)
tyn_eucl_dist <- tyn_eucl_dist[-c(12), -c(12)] # Drop sibling
tyn_eucl_dist <- as.matrix(tyn_eucl_dist)
tyn_eucl_dist <- tyn_eucl_dist[order(rownames(tyn_eucl_dist)),
                               order(colnames(tyn_eucl_dist))]

# Do a redundancy analysis to compare the distance in repeats and SNP
# First, make a pcoa based on the euclidean distance matrix
tyn_pcoa <- pcoa(tyn_eucl_dist)
tyn_X <- tyn_pcoa$vectors
# Make distance matrix based on repeats
# First, load full dataset
tyn_clusters <- read.delim("input/tyn_clusters.txt")
tyn_clusters <- subset(tyn_clusters, select = -c(GRW_3)) # Remove sibling
tyn_clusters <- subset(tyn_clusters, Annotation != "All/organelle/mitochondria")
tyn_clusters <- tyn_clusters[,-21]
tyn_clusters <- t(tyn_clusters)
for (i in 1:ncol(tyn_clusters)) {
  tyn_clusters[,i] <- as.numeric(tyn_clusters[,i])
}
tyn_rep_dist <- vegdist(tyn_clusters, method = "bray")
tyn_rep_dist <- as.matrix(tyn_rep_dist)
tyn_rep_dist <- tyn_rep_dist[order(rownames(tyn_rep_dist)),
                             order(colnames(tyn_rep_dist))]
tyn_rep_dist <- as.dist(tyn_rep_dist)
# Make a Moran Eigenvector's Maps based on the repeat distance matrix
tyn_dbmem <- dbmem(tyn_rep_dist)
# Run RDA
tyn_rda <- rda(tyn_X ~ ., data = tyn_dbmem)
# Get significance and Rsquared
tyn_Rsquared <- RsquareAdj(tyn_rda)
tyn_Rsquared
tyn_signif <- anova.cca(tyn_rda, permutations = 9999)
tyn_signif

# Permanova
tyn_permanova <- adonis2(tyn_rep_dist ~ tyn$pop,
                         data = tyn, permutations = 9999)
tyn_permanova

# Prepare for plotting
tyn_rep_dist <- as.matrix(tyn_rep_dist)
tyn_rep_dist[upper.tri(tyn_rep_dist, diag = TRUE)] <- NA
tyn_rep_dist <- melt(tyn_rep_dist, value.name = "rep_dist", na.rm = TRUE)
tyn_eucl_dist <- as.matrix(tyn_eucl_dist)
tyn_eucl_dist[upper.tri(tyn_eucl_dist, diag = TRUE)] <- NA
tyn_eucl_dist <- melt(tyn_eucl_dist, value.name = "eucl_dist", na.rm = TRUE)
tyn_dist <- merge(tyn_rep_dist, tyn_eucl_dist, by = c("Var1","Var2"))


#### Nivalis ####
## Load repeat dataset summed by annotation
niv_clusters <- read.delim("input/niv_clusters.txt")
niv_clusters <- subset(niv_clusters, select = -c(SAJ_3)) # Remove sibling
niv_clusters$Annotation[niv_clusters$Annotation == 'All'] <- 'Unindentified repeat'
niv_clusters$Annotation[niv_clusters$Annotation == 'All/repeat'] <- 'Unindentified repeat'
niv_clusters$Annotation[niv_clusters$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
niv_clusters$Annotation[niv_clusters$Annotation == 'All/repeat/mobile_element/Class_I'] <- 'Unindentified repeat'
niv_clusters$Annotation[grep("All/repeat/rDNA/*", niv_clusters$Annotation)] <- 'All/repeat/rDNA'
niv_clusters <- aggregate(. ~ Annotation, data = niv_clusters, FUN = sum)
niv_clusters <- t(niv_clusters)
colnames(niv_clusters) <- niv_clusters[1,]
niv_clusters <- niv_clusters[-1,]
niv_clusters <- niv_clusters[,-13] # Remove unidentified repeats
niv_clusters <- niv_clusters[,-1] # Remove organelle clusters
niv_clusters <- as.data.frame(niv_clusters)
for (i in 1:ncol(niv_clusters)) {
  niv_clusters[,i] <- as.numeric(niv_clusters[,i])
}
niv_clusters <- niv_clusters[order(row.names(niv_clusters)),]

# Calculate diversity indices
niv_simpson <- diversity(niv_clusters, index = "simpson")
niv_div <- as.data.frame(niv_simpson)
niv_div$indiv <- rownames(niv_div)

# Link with genetic diversity theta, based on individuals (with mlrho)
niv_theta <- read.delim("input/niv_theta.txt")
niv <- merge(niv_div, niv_theta, by = "indiv")
niv$pop <- c('gri', 'gri', 'gri', 'gro', 'gro', 'gro',
             'saj', 'saj', 'sch', 'sch')
niv$pop <- as.factor(niv$pop)

# correlation between simpson diversity and theta
niv_cor <- cor.test(niv$niv_simpson, niv$theta_avg, method = "spearman")
niv_cor

# Load genetic distance dataset
niv_eucl_dist <- read.delim("input/niv_eucl_dist.txt", row.names = 1)
niv_eucl_dist <- niv_eucl_dist[-c(4), -c(4)] # Drop sibling
niv_eucl_dist <- as.matrix(niv_eucl_dist)
niv_eucl_dist <- niv_eucl_dist[order(rownames(niv_eucl_dist)),
                               order(colnames(niv_eucl_dist))]

# Do a redundancy analysis to compare the distance in repeats and SNP
# First, make a pcoa based on the euclidean distance matrix
niv_pcoa <- pcoa(niv_eucl_dist)
niv_X <- niv_pcoa$vectors
# Make distance matrix based on repeats
# First, load full dataset
niv_clusters <- read.delim("input/niv_clusters.txt")
niv_clusters <- subset(niv_clusters, select = -c(SAJ_3)) # Remove sibling
niv_clusters <- subset(niv_clusters, Annotation != "All/organelle/mitochondria")
niv_clusters <- niv_clusters[,-11]
niv_clusters <- t(niv_clusters)
for (i in 1:ncol(niv_clusters)) {
  niv_clusters[,i] <- as.numeric(niv_clusters[,i])
}
niv_rep_dist <- vegdist(niv_clusters, method = "bray")
niv_rep_dist <- as.matrix(niv_rep_dist)
niv_rep_dist <- niv_rep_dist[order(rownames(niv_rep_dist)),
                             order(colnames(niv_rep_dist))]
niv_rep_dist <- as.dist(niv_rep_dist)
# Make a Moran Eigenvector's Maps based on the repeat distance matrix
niv_dbmem <- dbmem(niv_rep_dist)
# Run RDA
niv_rda <- rda(niv_X ~ ., data = niv_dbmem)
# Get significance and Rsquared
niv_Rsquared <- RsquareAdj(niv_rda)
niv_Rsquared
niv_signif <- anova.cca(niv_rda, permutations = 9999)
niv_signif

# Permanova
niv_permanova <- adonis2(niv_rep_dist ~ niv$pop,
                         data = niv, permutations = 9999)
niv_permanova

# Prepare for plotting
niv_rep_dist <- as.matrix(niv_rep_dist)
niv_rep_dist[upper.tri(niv_rep_dist, diag = TRUE)] <- NA
niv_rep_dist <- melt(niv_rep_dist, value.name = "rep_dist", na.rm = TRUE)
niv_eucl_dist <- as.matrix(niv_eucl_dist)
niv_eucl_dist[upper.tri(niv_eucl_dist, diag = TRUE)] <- NA
niv_eucl_dist <- melt(niv_eucl_dist, value.name = "eucl_dist", na.rm = TRUE)
niv_dist <- merge(niv_rep_dist, niv_eucl_dist, by = c("Var1","Var2"))

#### Pronoe ####
## Load repeat dataset summed by annotation
pro_clusters <- read.delim("input/pro_clusters.txt")
pro_clusters$Annotation[pro_clusters$Annotation == 'All'] <- 'Unindentified repeat'
pro_clusters$Annotation[pro_clusters$Annotation == 'All/repeat'] <- 'Unindentified repeat'
pro_clusters$Annotation[pro_clusters$Annotation == 'All/repeat/mobile_element'] <- 'Unindentified repeat'
pro_clusters$Annotation[pro_clusters$Annotation == 'All/repeat/mobile_element/Class_I'] <- 'Unindentified repeat'
pro_clusters$Annotation[grep("All/repeat/rDNA/*", pro_clusters$Annotation)] <- 'All/repeat/rDNA'
pro_clusters <- aggregate(. ~ Annotation, data = pro_clusters, FUN = sum)
pro_clusters <- t(pro_clusters)
colnames(pro_clusters) <- pro_clusters[1,]
pro_clusters <- pro_clusters[-1,]
pro_clusters <- pro_clusters[,-11] # Remove unidentified repeats
pro_clusters <- pro_clusters[,-1] # Remove organelle clusters
pro_clusters <- as.data.frame(pro_clusters)
for (i in 1:ncol(pro_clusters)) {
  pro_clusters[,i] <- as.numeric(pro_clusters[,i])
}
pro_clusters <- pro_clusters[order(row.names(pro_clusters)),]

# Calculate diversity indices
pro_simpson <- diversity(pro_clusters, index = "simpson")
pro_div <- as.data.frame(pro_simpson)
pro_div$indiv <- rownames(pro_div)

# Link with genetic diversity theta, based on individuals (with mlrho)
pro_theta <- read.delim("input/pro_theta.txt")
pro <- merge(pro_div, pro_theta, by = "indiv")
pro$pop <- c('psa', 'psa', 'psa','psa', 'psa',
             'ver', 'ver', 'ver', 'ver', 'ver', 'ver')
pro$pop <- as.factor(pro$pop)

# correlation between simpson diversity and theta
pro_cor <- cor.test(pro$pro_simpson, pro$theta_avg, method = "spearman")
pro_cor

# Load genetic distance dataset
pro_eucl_dist <- read.delim("input/pro_eucl_dist.txt", row.names = 1)
pro_eucl_dist <- as.matrix(pro_eucl_dist)
pro_eucl_dist <- pro_eucl_dist[order(rownames(pro_eucl_dist)),
                               order(colnames(pro_eucl_dist))]

# Do a redundancy analysis to compare the distance in repeats and SNP
# First, make a pcoa based on the euclidean distance matrix
pro_pcoa <- pcoa(pro_eucl_dist)
pro_X <- pro_pcoa$vectors
# Make distance matrix based on repeats
# First, load full dataset
pro_clusters <- read.delim("input/pro_clusters.txt")
pro_clusters <- subset(pro_clusters, Annotation != "All/organelle/mitochondria")
pro_clusters <- pro_clusters[,-12]
pro_clusters <- t(pro_clusters)
for (i in 1:ncol(pro_clusters)) {
  pro_clusters[,i] <- as.numeric(pro_clusters[,i])
}
pro_rep_dist <- vegdist(pro_clusters, method = "bray")
pro_rep_dist <- as.matrix(pro_rep_dist)
pro_rep_dist <- pro_rep_dist[order(rownames(pro_rep_dist)),
                             order(colnames(pro_rep_dist))]
pro_rep_dist <- as.dist(pro_rep_dist)
# Make a Moran Eigenvector's Maps based on the repeat distance matrix
pro_dbmem <- dbmem(pro_rep_dist)
# Run RDA
pro_rda <- rda(pro_X ~ ., data = pro_dbmem)
# Get significance and Rsquared
pro_Rsquared <- RsquareAdj(pro_rda)
pro_Rsquared
pro_signif <- anova.cca(pro_rda, permutations = 9999)
pro_signif

# Permanova
pro_permanova <- adonis2(pro_rep_dist ~ pro$pop,
                         data = pro, permutations = 9999)
pro_permanova

# Prepare for plotting
pro_rep_dist <- as.matrix(pro_rep_dist)
pro_rep_dist[upper.tri(pro_rep_dist, diag = TRUE)] <- NA
pro_rep_dist <- melt(pro_rep_dist, value.name = "rep_dist", na.rm = TRUE)
pro_eucl_dist <- as.matrix(pro_eucl_dist)
pro_eucl_dist[upper.tri(pro_eucl_dist, diag = TRUE)] <- NA
pro_eucl_dist <- melt(pro_eucl_dist, value.name = "eucl_dist", na.rm = TRUE)
pro_dist <- merge(pro_rep_dist, pro_eucl_dist, by = c("Var1","Var2"))

#### Final plot ####
cas_div_plot <- ggplot(cas, aes(x = cas_simpson, y = theta_avg)) +
  geom_point(size = 0.6) +
  xlab("Simpson diversity") +
  ylab(expression(paste(theta, " estimate"))) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = expression(paste("(a) ", italic("E. cassioides")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
cas_dist_plot <- ggplot(cas_dist, aes(x = rep_dist, y = eucl_dist)) +
  geom_point(size = 0.6) +
  geom_smooth(method = 'lm', color = "black", linewidth = 0.1, se = T) +
  xlab("Distance in repeats") +
  ylab("Individual based\ngenetic distance") +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  labs(title = expression(paste("(b) ", italic("E. cassioides")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
tyn_div_plot <- ggplot(tyn, aes(x = tyn_simpson, y = theta_avg)) +
  geom_point(size = 0.6) +
  xlab("Simpson diversity") +
  ylab(expression(paste(theta, " estimate"))) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = expression(paste(italic("E. tyndarus")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
tyn_dist_plot <- ggplot(tyn_dist, aes(x = rep_dist, y = eucl_dist)) +
  geom_point(size = 0.6) +
  geom_smooth(method = 'lm', color = "black", linewidth = 0.1, se = T) +
  xlab("Distance in repeats") +
  ylab("Individual based\ngenetic distance") +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  labs(title = expression(paste(italic("E. tyndarus")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
niv_div_plot <- ggplot(niv, aes(x = niv_simpson, y = theta_avg)) +
  geom_point(size = 0.6) +
  xlab("Simpson diversity") +
  ylab(expression(paste(theta, " estimate"))) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = expression(paste(italic("E. nivalis")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
niv_dist_plot <- ggplot(niv_dist, aes(x = rep_dist, y = eucl_dist)) +
  geom_point(size = 0.6) +
  geom_smooth(method = 'lm', color = "black", linewidth = 0.1, se = T) +
  xlab("Distance in repeats") +
  ylab("Individual based\ngenetic distance") +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  labs(title = expression(paste(italic("E. nivalis")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
pro_div_plot <- ggplot(pro, aes(x = pro_simpson, y = theta_avg)) +
  geom_point(size = 0.6) +
  geom_smooth(method = 'lm', color = "black", linewidth = 0.1, se = T) +
  xlab("Simpson diversity") +
  ylab(expression(paste(theta, " estimate"))) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = expression(paste(italic("E. pronoe")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))
pro_dist_plot <- ggplot(pro_dist, aes(x = rep_dist, y = eucl_dist)) +
  geom_point(size = 0.6) +
  geom_smooth(method = 'lm', color = "black", linewidth = 0.1, se = T) +
  xlab("Distance in repeats") +
  ylab("Individual based\ngenetic distance") +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  labs(title = expression(paste(italic("E. pronoe")))) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 8),
        panel.background = element_blank()) +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))

pdf("output/Fig_div_dist.pdf", width = 6.65354)
egg::ggarrange(cas_div_plot, tyn_div_plot, niv_div_plot, pro_div_plot,
               cas_dist_plot, tyn_dist_plot, niv_dist_plot, pro_dist_plot,
          label.args = list(gp = grid::gpar(fontface = "plain")),
          ncol = 4, nrow = 2)
dev.off()


#### Plot whole-genome PCoA ####
cbPalette <- c("#E69F00", "#0072B2", "#009E73", "#56B4E9",
               "#D55E00", "#CC79A7", "grey33")

cas_coord_pcoa <- as.data.frame(cas_pcoa$vectors[,1:2])
cas_coord_pcoa$pop <- c('GRF', 'GRF', 'GRF', 'GRW', 'GRW', 'GRW',
                        'KAN', 'KAN', 'KAN', 'PDM', 'PDM', 'PDM',
                        'ROU', 'ROU', 'ROU', 'SCH', 'SCH', 'SCH',
                        'STO', 'STO', 'STO')
tyn_coord_pcoa <- as.data.frame(tyn_pcoa$vectors[,1:2])
tyn_coord_pcoa$pop <- c('ARO', 'ARO', 'ARO', 'COM', 'COM', 'COM',
                        'GOT', 'GOT', 'GOT', 'GRF', 'GRF', 'GRF',
                        'GRW', 'GRW', 'LAU', 'LAU', 'LAU',
                        'SAN', 'SAN', 'SAN')
niv_coord_pcoa <- as.data.frame(niv_pcoa$vectors[,1:2])
niv_coord_pcoa$pop <- c('GRI', 'GRI', 'GRI', 'GRO', 'GRO',
                        'GRO', 'SAJ', 'SAJ', 'SCH', 'SCH')
pro_coord_pcoa <- as.data.frame(pro_pcoa$vectors[,1:2])
pro_coord_pcoa$pop <- c('psathura', 'psathura', 'psathura','psathura',
                        'psathura', 'vergy', 'vergy', 'vergy', 'vergy',
                        'vergy', 'vergy')

## To get square grid
cas_coord_pcoa_plot <- ggplot(cas_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(cas_pcoa$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(cas_pcoa$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = row.names(cas_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  coord_fixed() +
  scale_x_continuous(limits = c(-430, 430)) +
  scale_y_continuous(limits = c(-430, 430)) +
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

tyn_coord_pcoa_plot <- ggplot(tyn_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(tyn_pcoa$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(tyn_pcoa$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = row.names(tyn_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.1) +
  coord_fixed() +
  scale_x_continuous(limits = c(-500, 500)) +
  scale_y_continuous(limits = c(-500, 500)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(b) ", italic("Erebia tyndarus")))) +
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

niv_coord_pcoa_plot <- ggplot(niv_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(niv_pcoa$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(niv_pcoa$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Populations") +
  geom_text_repel(segment.size = 0.1, label = row.names(niv_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  coord_fixed() +
  scale_x_continuous(limits = c(-520, 520)) +
  scale_y_continuous(limits = c(-520, 520)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(c) ", italic("Erebia nivalis")))) +
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

pro_coord_pcoa_plot <- ggplot(pro_coord_pcoa,
                              aes(x = Axis.1, y = Axis.2, color = pop)) +
  geom_point(size = 0.8) +
  xlab(paste0("PCo1 (", round(pro_pcoa$values$Relative_eig[1]*100,
                              digits = 1), "%)")) +
  ylab(paste0("PCo2 (", round(pro_pcoa$values$Relative_eig[2]*100,
                              digits = 1), "%)")) +
  labs(color = "Subspecies") +
  geom_text_repel(segment.size = 0.1, label = row.names(pro_coord_pcoa),
                  show.legend = F, size = 2, force_pull = 2, force = 0.01,
                  max.time = 5, box.padding = 0.2) +
  coord_fixed() +
  scale_x_continuous(limits = c(-485, 485)) +
  scale_y_continuous(limits = c(-485, 485)) +
  scale_color_manual(values = cbPalette) +
  labs(title = expression(paste("(d) ", italic("Erebia pronoe")))) +
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


#### Final plot PCoA
pdf("output/Fig_PCoA_gen.pdf", width = 6.65354)
egg::ggarrange(cas_coord_pcoa_plot, tyn_coord_pcoa_plot,
               niv_coord_pcoa_plot, pro_coord_pcoa_plot,
               label.args = list(gp = grid::gpar(fontface = "plain")),
               ncol = 2, nrow = 2)
dev.off()
