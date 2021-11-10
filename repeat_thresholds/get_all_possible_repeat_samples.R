library(stringr)

# from a list of all possible PHO WGS Ids, get the repeat names
repeat_names <- read.csv(file.choose(), header = F, na.strings=c("","NA"),
                         stringsAsFactors=FALSE)

repeats <- data.frame(repeat_name = repeat_names[grepl("-v", repeat_names$V1) &
                                                   !grepl("twist|Neg|NTC", repeat_names$V1),])

repeats$original_name <- str_split_fixed(repeats$repeat_name, "-v", 2)[,1]

repeat_partners <- data.frame(all_names = unique(c(as.character(repeats$original_name),
                                            as.character(repeats$repeat_name))))

write.table(repeat_partners, "all_repeat_names.txt", row.names = F, quote = F)
