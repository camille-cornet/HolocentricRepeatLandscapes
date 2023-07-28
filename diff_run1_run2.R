rm(list = ls())

#### Comparison Erebia comparative run1 and run2 in terms of genome proportions
ere_run1_prop <- read.delim("input/Genome_proportion_ere_run1.txt")
ere_run2_prop <- read.delim("input/Genome_proportion_ere_run2.txt")

ere_prop_merged <- merge(ere_run1_prop, ere_run2_prop,
                         by = "Annotation", all = T)
ere_prop_merged$diff <-
  abs(ere_prop_merged$Proportion....x - ere_prop_merged$Proportion....y)
ere_prop_merged$diff

sum(ere_prop_merged$diff, na.rm = T)

write.table(ere_prop_merged, file = "output/ere_diff_run1_run2.txt", sep = "\t")
