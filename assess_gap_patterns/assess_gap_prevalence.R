library(stringr)
library(dplyr)
library(magrittr)
library(tidyr)
library(splitstackshape)
library(ggplot2)
library(ggrepel)

find_nextclade <- list.files("/NetDrive/Projects/COVID-19/Other/master_fasta/", pattern = "nextclade_all*")
nextclade_all <- read.csv(paste("/NetDrive/Projects/COVID-19/Other/master_fasta/", find_nextclade, sep=""), header = TRUE, na.strings=c("","NA"), stringsAsFactors=FALSE,
                          sep="\t")

# nextclade_all <- nextclade_all %>% filter(qc.overallStatus == "good")

gap_freqs <- data.frame(table(unlist(strsplit(as.vector(nextclade_all$missing),
                split = ","))))

gap_freqs$distance <- as.numeric(str_split_fixed(gap_freqs$Var1, "-", 2)[,2]) -
  as.numeric(str_split_fixed(gap_freqs$Var1, "-", 2)[,1]) + 1


tidy_gaps <- cSplit(nextclade_all %>% select(seqName, missing), splitCols = "missing", ",") %>%
                        # gather from the first missing column to the last one automatically
  gather("name", "gap", missing_001:.[,ncol(.)]) %>% select(-name) %>% filter(! is.na(gap)) %>%
  sample_n(nrow(.))

tidy_gaps$gap[is.na(tidy_gaps$gap)] <- "None"

# tidy_gaps$gap[is.na(tidy_gaps$gap)] <- "None"

tidy_gaps$distance <- as.numeric(as.numeric(str_split_fixed(tidy_gaps$gap, "-", 2)[,2]) -
  as.numeric(str_split_fixed(tidy_gaps$gap, "-", 2)[,1]) + 1)

tidy_gaps$distance[is.na(tidy_gaps$distance)] <- 1


as.data.frame(table(tidy_gaps$gap)) %>% arrange(-Freq)

dist_freqs <- tidy_gaps %>% group_by(distance) %>% summarise(count = n()) %>% filter(!is.na(distance))

ggplot(dist_freqs %>% filter(distance > 200), aes(x = distance, y = count)) + geom_point()

# samples flagged during NML submission as being the shortest
very_short <- c("PHLON21-SARS28095", "PHLON21-SARS28178", "PHLON21-SARES28111",
                "PHLON21-SARS28177")

tidy_gap_flagged <- tidy_gaps %>% filter(seqName %in% very_short)

# read a CSV file containing all of the WGS Ids that have been Gisaid accepted
all_gisaid_accepted <- read.csv(file.choose(), header = TRUE, na.strings=c("","NA"), stringsAsFactors=FALSE)

colnames(all_gisaid_accepted) <- "name"

all_gisaid_accepted$WGS_Id <- paste("PHLON", str_split_fixed(all_gisaid_accepted$name, "-", 5)[,4],
                                    "-SARS",
                                    str_split_fixed(str_split_fixed(all_gisaid_accepted$name, "-", 5)[,5],
                                    "/", 2)[,1], sep="")

all_gisaid_accepted <- all_gisaid_accepted %>% filter(!grepl("/2020-SARS", WGS_Id))

tidy_gaps <- tidy_gaps %>% mutate(category = ifelse(seqName %in% all_gisaid_accepted$WGS_Id, "in_gisaid",
                                                    "not_in_gisaid")) %>%
  filter(!grepl(
    "NTC|TWIST|pos|Neg|ntc|POS|twist|ON-PHL|SCTSF|Pos|blank|negative|MN908947|Twist", seqName)) %>%
  # cap the most recent sequence ID to the most recent one in Gisaid for comparison
  filter(seqName <= max(all_gisaid_accepted$WGS_Id))

d <- density(subset(tidy_gaps, category == "in_gisaid")$distance)
plot(d, main="Density of gap sizes, Gisaid accepted sequences")

# confirm non normality for the gap sizes for 
library(nortest)
ad.test(subset(tidy_gaps, category == "in_gisaid")$distance)$p.value

totals_gisaid <- tidy_gaps %>% distinct(seqName, category) %>% group_by(category) %>% summarise(totals = n())

tallies <- as.data.frame(table(tidy_gaps$distance, tidy_gaps$category))

length(unique(tidy_gaps$seqName))

nextclade_all$category <- ifelse(nextclade_all$seqName %in% all_gisaid_accepted$WGS_Id, "in_gisaid",
                                 "not_in_gisaid")

tidy_gaps %>% filter(distance != "None") %>% group_by(category) %>% summarise(mean_gap_length = mean(as.numeric(distance)))

frame_gap_freqs <- as.data.frame(tidy_gaps %>% filter(distance != "None") %>%
                                   group_by(category, gap, distance) %>%
                                   summarise(counts = n()) %>% arrange(-counts) %>%
                                   # ignore the common masking and prime ends of the sequence                             
  merge(totals_gisaid, by = "category") %>%
  mutate(percent = 100*(counts/totals)) %>% mutate(distance = as.numeric(distance)) %>%
    filter(!grepl("1-54|29837-29903", gap)))

ggplot(frame_gap_freqs %>% filter(distance <= 10000), aes(x = distance, y = percent, color = category)) +
  geom_point() +
  geom_text_repel(data = frame_gap_freqs %>% filter(distance <= 10000 & percent >= 5),
                  aes(label = paste(gap, " (", distance, ")", sep="")),
            size = 3.5, color = "black") +
  facet_wrap(~category) +
  ylim(c(0, max(frame_gap_freqs$percent))) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") + xlab("gap distance")

tidy_gaps %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(seqName, category) %>% summarise(num_gaps = n()) %>% group_by(category) %>%
  summarise(mean = mean(num_gaps))

tidy_gaps %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(seqName, category) %>%
  summarise(avg_gap_sam = mean(as.numeric(distance))) %>% group_by(category) %>%
  summarise(mean = mean(avg_gap_sam))

### gap size median

tidy_gaps %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(category) %>%
  summarise(median = median(as.numeric(distance)))

tidy_gaps %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(seqName, category) %>%
  summarise(avg_gap_sam = mean(as.numeric(distance))) %>% group_by(category) %>%
  summarise(dev = sd(avg_gap_sam))


library(karyoploteR)
library(BSgenome)

cov2 <- data.frame(chr = c("ORF1a", "ORF1b", "Spike", "ORF3a",
                           "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10",
                           "5_UTR", "3_UTR"),
                   start = c(266, 13469, 21563, 25393, 26245, 26523, 27202,
                             27394, 27894, 28274, 29558, 1, 29675),
                   end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387,
                           27759, 28259, 29533, 29674, 266, 29903))

single_genome <- data.frame(chr = "ncov", start = 1, end = 29903)

tidy_gaps_gisaid <- tidy_gaps %>% filter(category == "in_gisaid") %>% sample_n(nrow(.))

regions_frame <- data.frame(chr = "ncov",
                            start = as.numeric(str_split_fixed(tidy_gaps_gisaid$gap, "-", 2)[,1]),
                            end = as.numeric(str_split_fixed(tidy_gaps_gisaid$gap, "-", 2)[,2]))

regions_frame <- regions_frame %>% mutate(end = ifelse(is.na(end), start, end))
# regions_frame$width <- regions_frame$end - regions_frame$start + 1
row.names(regions_frame) <- NULL

kp <- plotKaryotype(genome = single_genome, main = "Gap mapping for PHO sequences in Gisaid",
                    cex = 1.5)

data_for_plot <- paste(regions_frame$chr, ":", regions_frame$start, "-",
                       regions_frame$end, sep="")

# ensure that the trimmed ends are not included or for sequences that cannot
# identify gaps
data_for_plot <- data_for_plot[!grepl("NA|1-54|29837-29903", data_for_plot)]

kpPlotRegions(kp, data=data_for_plot, layer.margin = 0.01, r0=0, r1=0.4)

kpPlotCoverage(kp, data = data_for_plot, r0=-0.3, r1=-0.1)


# do not show markers for the largest genes, add them as rectangles
markers <- data.frame(chr = "ncov", pos = cov2$start, labels = cov2$chr) %>%
  filter(! labels %in% c("ORF1a", "ORF1b", "Spike"))

kpPlotMarkers(kp, chr=markers$chr, x=markers$pos, labels=markers$labels)

genome_regions <- paste("ncov", ":", cov2$start, "-",
                                         cov2$end, sep="")

library(RColorBrewer)
n <- length(unique(markers$labels))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


kpPlotRegions(kp, data=genome_regions, col=col_vector, r0=0.45, r1=0.6)

kpRect(kp, chr="ncov", x0 = (266 + 100), x1 = (13468 - 100), y0 = 0.2, y1 = 0.8, r0 = 0.8, r1 = 1)
kpText(kp, chr="ncov", x = (266 + 13468)/2, col="red", y = 0.5, r0 = 0.8, r1 = 1, labels="ORF1a")
kpRect(kp, chr="ncov", x0 = (13469 + 100), x1 = (21555 - 100), y0 = 0.2, y1 = 0.8, r0 = 0.8, r1 = 1)
kpText(kp, chr="ncov", x = (13469 + 21555)/2, col="red", y = 0.5, r0 = 0.8, r1 = 1, labels="ORF1b")
kpRect(kp, chr="ncov", x0 = (21563 + 100), x1 = (25384 - 300), y0 = 0.2, y1 = 0.8, r0 = 0.8, r1 = 1)
kpText(kp, chr="ncov", x = (21563 + 25384)/2, col="red", y = 0.5, r0 = 0.8, r1 = 1, labels="Spike")

kpText(kp, chr="ncov", x = 8700, col="black", y = 0.5, r0 = 0.8, r1 = 1.5,
       labels="Note: Trimmed gaps from the UTRs, 1-54 & 29837-29903 are not included", cex = 0.9)

# manually annotate most prevalent gap locations given by previous plot

kpPlotMarkers(kp, chr="ncov", x=21991, labels="21991 to 21993", cex = 0.7, y = 0.55,
              text.orientation = "horizontal")
kpRect(kp, chr="ncov", x0 = 3508, x1 = 3795, y0 = 0.2, y1 = 0.8, r0 = 0.6, r1 = 0.8,
       col = "black")
kpText(kp, chr="ncov", x = 4950, col="black", y = 0.5, r0 = 0.7, r1 = 0.7,
       labels="3508 to 3795", cex = 0.8)


kpRect(kp, chr="ncov", x0 = 19276, x1 = 19570, y0 = 0.2, y1 = 0.8, r0 = 0.6, r1 = 0.8,
       col = "black")
kpText(kp, chr="ncov", x = 20500, col="black", y = 0.5, r0 = 0.7, r1 = 0.7,
       labels=paste("19276 to ", "19570", sep="\n"), cex = 0.8)


kpRect(kp, chr="ncov", x0 = 22325, x1 = 22542, y0 = 0.2, y1 = 0.8, r0 = 0.6, r1 = 0.8,
       col = "black")
kpText(kp, chr="ncov", x = 23900, col="black", y = 0.5, r0 = 0.7, r1 = 0.7,
       labels="22325 to 22542", cex = 0.8)

kpText(kp, chr="ncov", x = 1250, col="black", y = 0.5, r0 = -0.15, r1 = -0.15,
       labels="per base density", cex = 0.8)

kpText(kp, chr="ncov", x = 1300, col="black", y = 0.5, r0 = 0.3, r1 = 0.3,
       labels="gap contigs", cex = 0.8)


### compare the PHO sequence lengths 

pho_lengths <- read.csv(file.choose(), header = TRUE, na.strings=c("","NA"),
                        stringsAsFactors=FALSE)

pho_lengths <- pho_lengths %>% filter(!grepl(
  "NTC|TWIST|pos|Neg|ntc|POS|twist|ON-PHL|SCTSF|Pos|blank|negative|MN908947|Twist", sample) &
    length > 0)


tidy_gaps_with_length <- merge(tidy_gaps_gisaid, pho_lengths, by.x = "seqName",
                               by.y = "sample")

tidy_gaps_with_length <- tidy_gaps_with_length %>% mutate(cat = ifelse(length <= 29885, "29825_to_29885",
                                    ifelse(length > 29885 & length <= 29896, "29885_to_29896",
                                           "over_29896"))) %>% filter(!grepl("NA",
                                                                             gap))
# ggplot(tidy_gaps_with_length, aes(x = cat, y = length, fill = cat)) + geom_boxplot()


gisaid_length_table <- tidy_gaps_with_length[!duplicated(tidy_gaps_with_length$seqName),] %>%
  group_by(cat) %>% summarise(counts = n())
colnames(gisaid_length_table) <- c("cat", "totals")

frame_freqs_length <- as.data.frame(tidy_gaps_with_length %>% filter(distance != "None") %>%
                                   group_by(cat, gap, distance) %>%
                                   summarise(counts = n()) %>% arrange(-counts) %>%
                                   # ignore the common masking and prime ends of the sequence                             
                                   filter(!grepl("1-54|29837-29903", gap))) %>%
  merge(gisaid_length_table, by = "cat") %>%
  mutate(percent = 100*(counts/totals)) %>% mutate(distance = as.numeric(distance))

ggplot(frame_freqs_length, aes(x = distance, y = percent, color = cat)) +
  geom_point() +
  geom_text_repel(data = frame_freqs_length %>% filter(percent >= 5),
                  aes(label = paste(gap, " (", distance, ")", sep="")),
                  size = 4, color = "black") +
  facet_wrap(~cat) +
  ylim(c(0, max(frame_freqs_length$percent) + 2.5)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=10),
                     strip.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "none") + xlab("Gap Distance (bp)") + ylab("Sequence Percentage (%)") +
  ggtitle("Gap prevalence by sequence length, Gisaid accepted")

for (i in unique(tidy_gaps_with_length$cat)) {
  subset_length <- subset(tidy_gaps_with_length, cat == i) %>% filter(!grepl("NA|1-54|29837-29903",
                                                                  gap))
  length_frame <- data.frame(chr = "ncov",
                              start = as.numeric(str_split_fixed(subset_length$gap, "-", 2)[,1]),
                              end = as.numeric(str_split_fixed(subset_length$gap, "-", 2)[,2]))
  
  length_frame <- length_frame %>% mutate(end = ifelse(is.na(end), start, end))
  # regions_frame$width <- regions_frame$end - regions_frame$start + 1
  row.names(length_frame) <- NULL
  
  data_for_plot <- paste(length_frame$chr, ":", length_frame$start, "-",
                         length_frame$end, sep="")
  
  plot_name <- paste("data_for_plot_", i, sep="")
  assign(plot_name, data_for_plot)
  
}

kp <- plotKaryotype(genome = single_genome, main = "Gap mapping for PHO sequences in Gisaid by length",
                    cex = 1.5)

kpPlotRegions(kp, data=data_for_plot_29825_to_29885, r0=0, r1=0.3,
              col = "red")
kpPlotRegions(kp, data=data_for_plot_29885_to_29896, r0=0.3, r1=0.6,
              col = "dark green")
kpPlotRegions(kp, data=data_for_plot_over_29896, r0=0.6, r1=0.9,
              col = "blue")

kpPlotMarkers(kp, chr=markers$chr, x=markers$pos, labels=markers$labels)

kpRect(kp, chr="ncov", x0 = (266 + 100), x1 = (13468 - 100), y0 = 0.2, y1 = 0.8, r0 = 0.95, r1 = 1.15)
kpText(kp, chr="ncov", x = (266 + 13468)/2, col="red", y = 0.5, r0 = 0.95, r1 = 1.15, labels="ORF1a")
kpRect(kp, chr="ncov", x0 = (13469 + 100), x1 = (21555 - 100), y0 = 0.2, y1 = 0.8, r0 = 0.95, r1 = 1.15)
kpText(kp, chr="ncov", x = (13469 + 21555)/2, col="red", y = 0.5, r0 = 0.95, r1 = 1.15, labels="ORF1b")
kpRect(kp, chr="ncov", x0 = (21563 + 100), x1 = (25384 - 300), y0 = 0.2, y1 = 0.8, r0 = 0.95, r1 = 1.15)
kpText(kp, chr="ncov", x = (21563 + 25384)/2, col="red", y = 0.5, r0 = 0.95, r1 = 1.15, labels="Spike")

kpText(kp, chr="ncov", x = 8700, col="black", y = 0.5, r0 = 1.2, r1 = 1.4,
       labels="Note: Trimmed gaps from the UTRs, 1-54 & 29837-29903 are not included", cex = 0.9)

# manually annotate most prevalent gap locations given by previous plot

kpPlotMarkers(kp, chr="ncov", x=21991, labels="21991 to 21993", cex = 0.7, y = 0.85,
              text.orientation = "horizontal")
# kpRect(kp, chr="ncov", x0 = 3508, x1 = 3795, y0 = 0.2, y1 = 0.8, r0 = 0.9, r1 = 1,
#        col = "black")
kpText(kp, chr="ncov", x = 5200, col="black", y = 0.5, r0 = 0.9, r1 = 0.9,
       labels="3508 to 3795", cex = 0.8)


# kpRect(kp, chr="ncov", x0 = 19276, x1 = 19570, y0 = 0.2, y1 = 0.8, r0 = 0.8, r1 = 0.9,
#        col = "black")
# kpText(kp, chr="ncov", x = 20500, col="black", y = 0.5, r0 = 0.5, r1 = 0.5,
#        labels=paste("19276 to ", "19570", sep="\n"), cex = 0.8)


# kpRect(kp, chr="ncov", x0 = 22325, x1 = 22542, y0 = 0.2, y1 = 0.8, r0 = 0.6, r1 = 0.8,
#        col = "black")
kpText(kp, chr="ncov", x = 23900, col="black", y = 0.5, r0 = 0.2, r1 = 0.2,
       labels="22325 to 22542", cex = 0.8)

kpText(kp, chr="ncov", x = -1300, col="black", y = 0.5, r0 = 0.1, r1 = 0.1,
       labels="29825 to 29885", cex = 0.8)

kpText(kp, chr="ncov", x = -1300, col="black", y = 0.5, r0 = 0.4, r1 = 0.4,
       labels="29885 to 29896", cex = 0.8)

kpText(kp, chr="ncov", x = -1100, col="black", y = 0.5, r0 = 0.7, r1 = 0.7,
       labels="over 29896", cex = 0.8)


#### add lineage assignments to lengths

lineages <- read.csv(file.choose(), header = TRUE, na.strings=c("","NA"),
                     stringsAsFactors=FALSE)

tidy_gaps_with_lineages <- merge(tidy_gaps_with_length, lineages %>%
                                   select(taxon, lineage), by.x = "seqName",
                                 by.y = "taxon")

all_gisaid_lin_counts <- merge(all_gisaid_accepted %>% select(WGS_Id), lineages %>%
                                 select(taxon, lineage), by.x = "WGS_Id",
                               by.y = "taxon") %>% group_by(lineage) %>%
  summarise(counts = n())



gaps_lin_freqs <- tidy_gaps_with_lineages[!duplicated(tidy_gaps_with_lineages[1]),] %>%
                  group_by(cat, lineage) %>% summarise(counts = n()) %>%
  filter(counts > 200) %>% mutate(lineage = ifelse(grepl("AY", lineage),
                                                   "B.1.617.2", lineage))

length(unique(gaps_lin_freqs$lineage)) 
ggplot(gaps_lin_freqs, aes(x = lineage, y = counts, fill = cat)) + geom_bar(stat = "identity",
                                                                                position = "stack") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12)) +
  xlab("Pango lineage") +labs(fill = "Sequence length") +
  coord_flip()


tidy_gaps_with_lineages[!duplicated(tidy_gaps_with_lineages$seqName),] %>% group_by(lineage) %>% summarise(avg = mean(length)) %>%
  filter(lineage %in% c("B.1.1.7", "P.1", "B.1.617.2", "P.1", "B.1.351"))


tidy_gaps_with_length %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(seqName, cat) %>%
  summarise(num_gaps = n()) %>% group_by(cat) %>%
  summarise(mean = mean(num_gaps))

tidy_gaps_with_length %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(seqName, cat) %>%
  summarise(avg_gap_sam = mean(as.numeric(distance))) %>% group_by(cat) %>%
  summarise(mean = mean(avg_gap_sam))

### gap size median

tidy_gaps_with_length %>% filter(distance != "None", !grepl("1-54|29837-29903", gap)) %>% group_by(cat) %>%
  summarise(avg_gap_sam = median(as.numeric(distance)))



