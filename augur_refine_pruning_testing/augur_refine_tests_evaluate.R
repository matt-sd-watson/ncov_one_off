library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(ggpubr)
library(corrplot)

testing <- read.csv("/home/mwatson/COVID-19/one_off/augur_refine_tests/augur_refine_test_results.csv",
                    header = TRUE, na.strings=c("","NA"),
                    stringsAsFactors=FALSE)

# lineages <- read.csv(file.choose(), header = TRUE, na.strings=c("","NA"),
#                      stringsAsFactors=FALSE)

testing$wgs_id <- paste("PHLON",
                        str_split_fixed(testing$sample, "-", 4)[,3],
                        "-SARS", str_split_fixed(testing$sample, "-", 4)[,4],
                        sep = "")

nextstrain_available <- read.table(
  file.choose(),
  header = T,
  sep = ',',
  fill = TRUE,
  quote = ""
)


merged <- merge(testing, subset(nextstrain_available, select = c(WGS_Id, Pango.Lineage)),
                by.x = "wgs_id", by.y = "WGS_Id",
                all.x = T) %>% distinct()

new_high <- subset(merged, subset_number >= 20 & clock_value >= 4)

old_high <- subset(merged, subset_number <= 20 & clock_value >= 4)

table(new_high$Pango.Lineage)

table(old_high$Pango.Lineage)

high_numbers <- length(unique(merged$subset_number[merged$subset_number >= 20]))
not_just_10 <- length(unique(merged$subset_number[merged$subset_number > 20 & merged$subset_number < 200]))

numbers_over_clock_4 <- subset(merged, clock_value >= 2) %>% 
  mutate(cat = ifelse(subset_number <= 20, "Oldest 20", "New/Recent")) %>% group_by(cat, clock_value) %>%
  dplyr::summarize(size = n()) %>% distinct() %>% mutate(num = ifelse(cat == "Oldest 20", 20,
                                                                ifelse(cat == "New/Recent" &
                                                                  clock_value >= 2 & clock_value <= 9, 
                                                                not_just_10, high_numbers))) %>%
  mutate(avg = size /num)

# 
# group_counts <-  subset(merged, clock_value >= 2) %>% mutate(cat = ifelse(subset_number <= 20,
#                                                                           "Oldest 20", "New/Recent")) %>%
#   group_by(cat, clock_value) %>% dplyr::summarize(counts = n())
# 
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# df2 <- data_summary(group_counts, varname="counts", 
#                     groupnames=c("cat", "clock_value"))
# 
# std_by_cat<- subset(merged, clock_value >= 2) %>% mutate(cat = ifelse(subset_number <= 20,
#                                                                              "Oldest 20", "New/Recent")) %>%
#   group_by(cat, clock_value, subset_number) %>% summarize(counts = n()) %>%
#   group_by(cat, clock_value) %>% summarise(dev = sd(counts))
# 
# numbers_over_clock_4 <- merge(numbers_over_clock_4, std_by_cat, by=c("cat", "clock_value"),
#                               all = T)

ggplot(numbers_over_clock_4, aes(x = clock_value, y = avg, fill = cat)) + geom_bar(stat = "identity",
                                                                                   position = "dodge") +
  scale_x_continuous(breaks= seq(2, 10, 1)) +
  ggtitle("Average number of pruned sequences per clock filter, all lineages") +
  theme(panel.background = element_blank()) + scale_fill_manual(values = c("light blue", "grey30")) +
  labs(fill='Time Category') + ylab("Average Number") + xlab("TimeTree Refine Clock Value")
  # geom_errorbar(aes(ymin = avg, ymax = avg+dev), width = 0.2,
  #               position=position_dodge(.9))
  
# subset(merged, subset_number <= 20 & clock_value == 10)
# nrow(subset(merged, subset_number >= 20 & clock_value == 10))

# only B.1.617.2

numbers_over_clock_4_delta <- subset(merged, clock_value >= 2 & Pango.Lineage == "B.1.617.2 (Delta)") %>% 
  mutate(cat = ifelse(subset_number <= 20, "Oldest 20", "New/Recent")) %>% group_by(cat, clock_value) %>%
  dplyr::summarize(size = n()) %>% distinct() %>% mutate(num = ifelse(cat == "Oldest 20", 20,
                                                                          ifelse(cat == "New/Recent" &
                                                                                   clock_value >= 2 & clock_value <= 9, 
                                                                                 not_just_10, high_numbers))) %>%
  mutate(avg = size /num)

ggplot(numbers_over_clock_4_delta, aes(x = clock_value, y = avg, fill = cat)) + geom_bar(stat = "identity",
                                                                                   position = "dodge") +
  scale_x_continuous(breaks= seq(2, 10, 1)) +
  ggtitle("Average number of pruned sequences that are B.1.617.2 (Delta)") +
  theme(panel.background = element_blank())  + scale_fill_manual(values = c("light blue", "grey30")) +
  labs(fill='Time Category') + ylab("Average Number") + xlab("TimeTree Refine Clock Value")

merged_2 <- merge(numbers_over_clock_4_delta, numbers_over_clock_4, by = c("cat", "clock_value"),
                  all.x = T) %>% mutate(avg_delta = 100*(size.x/size.y))

ggplot(merged_2, aes(x = clock_value, y = avg_delta, fill = cat)) + geom_bar(stat = "identity",
                                                                            position = "dodge") +
  scale_x_continuous(breaks= seq(2, 10, 1)) +
  ggtitle("Average Percentage of pruned sequences that are B.1.617.2 (Delta)") +
  theme(panel.background = element_blank()) + scale_fill_manual(values = c("light blue", "grey30")) +
  labs(fill='Time Category') + ylab("Average Percentage (%)") + xlab("TimeTree Refine Clock Value") 


only_delta_high <- subset(merged, clock_value >= 10 & Pango.Lineage == "B.1.617.2 (Delta)")

most_pruned_delta <- as.data.frame(table(only_delta_high$wgs_id)) %>% arrange(-Freq)

only_delta_high_unique <- only_delta_high[!duplicated(only_delta_high$wgs_id),]

setwd("/home/mwatson/")
write.table(only_delta_high_unique$wgs_id, "only_delta_high.txt", row.names = F,
            quote = F)

# ## Nextclade outputs of the delta variants that have been filtered at clock value of 10
# 
# nextclade_delta_high <- read.csv(file.choose(), header = TRUE, na.strings=c("","NA"),
#                                  stringsAsFactors=FALSE, sep=";")
# 
# frame_missing_regions <- data.frame(table(unlist(strsplit(as.vector(nextclade_delta_high$missing),
#                       split = ","))))
# 
# # look for samples that have missing regions starting around 21717, for amplicon 71_left in the v3
# # primer set
# 
# frame_missing_regions_amp_71 <- frame_missing_regions[grepl("21717", frame_missing_regions$Var1),]
# sum(frame_missing_regions_amp_71$Freq)
# 
# table(nextclade_delta_high$aaDeletions)

### get all the delta variants and compare pruned vs non-pruned


lineages_delta <- subset(nextstrain_available, Pango.Lineage == "B.1.617.2 (Delta)")

nextclade <- list.files("/NetDrive/Projects/COVID-19/Other/master_fasta/", pattern = "*nextclade_all*")
nextclade_data <-read.table(paste("/NetDrive/Projects/COVID-19/Other/master_fasta/", nextclade,
                                  sep=""), header = T, sep = '\t', fill = TRUE, quote = "") 

nextclade_delta_all <- nextclade_data[nextclade_data$seqName %in% lineages_delta$WGS_Id & nextclade_data$clade == "21A (Delta)",]
nextclade_delta_all$pruned <- ifelse(nextclade_delta_all$seqName %in% only_delta_high_unique$wgs_id, "yes",
                                     "no")

prune_status_table <- data.frame(table(nextclade_delta_all$pruned))

for (i in unique(nextclade_delta_all$pruned)) {
  subset_cat <- subset(nextclade_delta_all, pruned == i)
  print(i)
  print(nrow(subset_cat[grepl("21717", subset_cat$missing) & grepl("21990", subset_cat$missing),])/nrow(subset_cat))
}


delta_pruned_nextclade <- subset(nextclade_delta_all, pruned == "yes" & totalAminoacidDeletions <= 5)
delta_no_prune_nextclade <- subset(nextclade_delta_all, pruned == "no" &
                                     seqName %in% nextstrain_available$WGS_Id &
                                     totalAminoacidDeletions <= 5)

missing_delta_pruned <- unlist(strsplit(as.vector(delta_pruned_nextclade$missing),
                split = ","))

missing_delta_not_pruned <- unlist(strsplit(as.vector(delta_no_prune_nextclade$missing),
                                            split = ","))

table_delta_pruned <- data.frame(table(missing_delta_pruned))
table_delta_not_pruned <- data.frame(table(missing_delta_not_pruned))

table_delta_pruned$perc <- 100*(table_delta_pruned$Freq / nrow(delta_pruned_nextclade))
table_delta_not_pruned$perc <- 100*(table_delta_not_pruned$Freq / nrow(delta_no_prune_nextclade))

# sum(table_delta_pruned[grepl("21717", table_delta_pruned$missing_delta_pruned),]$perc)
# sum(table_delta_not_pruned[grepl("21717", table_delta_not_pruned$missing_delta_not_pruned),]$perc)

table_delta_pruned$miss_start <- as.numeric(str_split_fixed(table_delta_pruned$missing_delta_pruned, "-", 2)[,1])
table_delta_pruned$miss_end <- as.numeric(str_split_fixed(table_delta_pruned$missing_delta_pruned, "-", 2)[,2])

table_delta_not_pruned$miss_start <- as.numeric(str_split_fixed(table_delta_not_pruned$missing_delta_not_pruned, "-", 2)[,1])
table_delta_not_pruned$miss_end <- as.numeric(str_split_fixed(table_delta_not_pruned$missing_delta_not_pruned, "-", 2)[,2])


table_pruned_amp_71 <- subset(table_delta_pruned, miss_start == 21717 & miss_end >= 22034)
sum(table_pruned_amp_71$perc)


table_not_pruned_amp_71 <- subset(table_delta_not_pruned, miss_start == 21717 & miss_end >= 22034)
sum(table_not_pruned_amp_71$perc)

t.test(delta_pruned_nextclade$totalSubstitutions, delta_no_prune_nextclade$totalSubstitutions)

t.test(delta_pruned_nextclade$totalAminoacidDeletions, delta_no_prune_nextclade$totalAminoacidDeletions)

# test for normality-= these data are not normally distributed

shapiro.test(delta_pruned_nextclade$totalAminoacidDeletions)
shapiro.test(delta_no_prune_nextclade$totalAminoacidDeletions)

# wilcoxon rank sum test for non paired non parametric


# for (i in seq(1, 100, 1)) {
#   random_no_prune <- sample_n(delta_no_prune_nextclade, 200)
#   
#   sub_counts_prune <- data.frame(counts = delta_pruned_nextclade$totalAminoacidDeletions,
#                                  cat = "pruned")
#   sub_counts_not_prune <- data.frame(counts = random_no_prune$totalAminoacidDeletions,
#                                      cat = "not pruned")
#   
#   
#   del_counts_all <- rbind(sub_counts_prune, sub_counts_not_prune)
#   
#   ggboxplot(del_counts_all, x = "cat", y = "counts", 
#             color = "cat", palette = c("#00AFBB", "#E7B800"),
#             ylab = "Deletions", xlab = "Prune Status")
#   
#   
#   res <- wilcox.test(counts ~ cat, data = del_counts_all,
#               alternative = "less", exact = F)
#   print(res$p.value)
# }

random_no_prune <- sample_n(delta_no_prune_nextclade, nrow(delta_no_prune_nextclade))

sub_counts_prune <- data.frame(counts = delta_pruned_nextclade$totalAminoacidDeletions,
                               cat = "pruned")
sub_counts_not_prune <- data.frame(counts = random_no_prune$totalAminoacidDeletions,
                               cat = "not pruned")


del_counts_all <- rbind(sub_counts_prune, sub_counts_not_prune)

ggboxplot(del_counts_all, x = "cat", y = "counts", 
          color = "cat", palette = c("#00AFBB", "#E7B800"),
          ylab = "Deletions", xlab = "Prune Status")


wilcox.test(counts ~ cat, data = del_counts_all,
            alternative = "less", exact = F)

missing_delta_pruned <- unlist(strsplit(as.vector(delta_pruned_nextclade$aaDeletions),
                                        split = ","))

missing_delta_not_pruned <- unlist(strsplit(as.vector(delta_no_prune_nextclade$aaDeletions),
                                            split = ","))

table_delta_pruned <- data.frame(table(missing_delta_pruned))
table_delta_not_pruned <- data.frame(table(missing_delta_not_pruned))

table_delta_pruned$perc <- 100*(table_delta_pruned$Freq / nrow(delta_pruned_nextclade))
table_delta_not_pruned$perc <- 100*(table_delta_not_pruned$Freq / nrow(delta_no_prune_nextclade))

table_delta_pruned$cat <- "pruned"
table_delta_not_pruned$cat <- "not pruned"
colnames(table_delta_pruned) <- c("del", "freq", "Percent", "category")
colnames(table_delta_not_pruned) <- colnames(table_delta_pruned)

rbind_del <- rbind(table_delta_pruned, table_delta_not_pruned) %>% filter(Percent > 1)
rbind_del$n_cat <- ifelse(rbind_del$category == "pruned", paste(rbind_del$category, " (",
                                                                nrow(only_delta_high_unique), ")",
                                                                sep=""),
                          paste(rbind_del$category, " (",
                                nrow(delta_no_prune_nextclade), ")",
                                sep=""))

ggplot(rbind_del, aes(x = del, y = Percent, fill = n_cat)) + geom_bar(stat = "identity",
                                                                         position = "dodge") +
  theme(panel.background = element_blank()) +
  labs(fill='Prune Status') + ylab("Percentage (%)") + xlab("AA Deletion") +
  ggtitle("Percentage of Nextstrain-eligible Delta sequences with AA Deletions")
  # geom_text(aes(label = n_lab), hjust = 0.5, size = 3)


## chi squared test

random_no_prune <- sample_n(delta_no_prune_nextclade, nrow(only_delta_high_unique))

pruned_with_del <- nrow(delta_pruned_nextclade[grepl("S:E156-|S:F157-", delta_pruned_nextclade$aaDeletions),])
pruned_no_del <- nrow(delta_pruned_nextclade) - pruned_with_del
pruned_counts <- c(pruned_with_del, pruned_no_del)

no_prune_with_del <- nrow(random_no_prune[grepl("S:E156-|S:F157-", random_no_prune$aaDeletions),])
no_prune_no_del <- nrow(random_no_prune) - no_prune_with_del
no_prune_counts <- c(no_prune_with_del, no_prune_no_del)

frame_for_chi_del <- data.frame(pruned = pruned_counts,
                                not_pruned = no_prune_counts)
row.names(frame_for_chi_del) <- c("deletion", "no deletion")

chi_test <- chisq.test(frame_for_chi_del)
chi_test
chi_test$observed
chi_test$expected

corrplot(chi_test$residuals, is.cor = FALSE)

## ch square contributions

contrib <- 100*chi_test$residuals^2/chi_test$statistic
corrplot(contrib, is.cor = FALSE)

