library(stringr)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

# dists <- read.csv(file.choose(), header = F, na.strings=c("","NA"), stringsAsFactors=FALSE)
# 
# dists$original_name <- str_split_fixed(dists$V2, "-v", 2)[,1]
# 
# dists_repeats <- subset(dists, V1 == original_name & V1 != V2 & V2 != "MN908947" & V2 != "MN908947" &
#                           !grepl("-v", V1))
# 
# dist_freqs <- data.frame(table(dists_repeats$V3))
# colnames(dist_freqs) <- c("No.SNPs", "num_repeat_pairs")
# 
# write.csv(dist_freqs, "repeat_snp_dists_frequencies.csv", row.names = F, quote = F)
# write.csv(dists_repeats %>% select(-original_name), "repeats_snp_dists.csv",
#           row.names = F, quote = F)

# analyze the frequency of SNP dists among all repeat groups
# as well as Ct values correlation to SNP dists 

only_repeat_partners <- read.csv(file.choose(), header = T, na.strings=c("","NA"), stringsAsFactors=FALSE)

wgs_with_ct <- read.csv(file.choose(), header = T, na.strings=c("","NA"), stringsAsFactors=FALSE)

only_repeat_partners <- merge(only_repeat_partners, wgs_with_ct, by.x = "V1", by.y = "WGS_Id",
                              all.x = T) %>% plyr::rename(c('Egene_Ct_rerun'='Ct_original'))

only_repeat_partners <- merge(only_repeat_partners, wgs_with_ct, by.x = "V2", by.y = "WGS_Id",
                              all.x = T) %>% plyr::rename(c('Egene_Ct_rerun'='Ct_repeat'))

grouped_by_snps <- only_repeat_partners %>% group_by(V3) %>% summarise(counts = n(), mean_original = mean(Ct_original,
                                                                                            na.rm = T),
                                                                       mean_repeat = mean(Ct_repeat,
                                                                                          na.rm = T),
                                                                       dev_original = sd(Ct_original,
                                                                                         na.rm = T),
                                                                       dev_repeat = sd(Ct_repeat,
                                                                                       na.rm = T))

grouped_by_snps$percent <- 100*(grouped_by_snps$counts/sum(grouped_by_snps$counts))

ggplot(grouped_by_snps, aes(x = as.numeric(V3),
                            y = percent,
                            )) + geom_bar(stat = "identity",
                                          fill = "blue") +
  xlim(c(-1, 100)) + xlab("SNP distance between repeats") +
  ylab("Repeat percentage (%)") + ylim(c(0, 110))

# make into tidy format before plotting
tidy_snps <- grouped_by_snps %>% select(V3, mean_original, mean_repeat) %>%
  pivot_longer(c("mean_original", "mean_repeat"), names_to = "category", values_to = "Ct")

# reverse the pivot
# pivot_wider(tidy_snps, names_from = category, values_from = Ct)

ggplot(tidy_snps, aes(x = as.factor(V3), y = Ct, fill = category)) + geom_bar(stat = "identity",
                                                                   position = "dodge") +
  xlab("SNP distance between repeats")

# ggplot(tidy_snps, aes(x = as.factor(V3), y = Ct, fill = category)) + geom_boxplot()
