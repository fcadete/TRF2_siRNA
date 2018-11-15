
library(tidyverse)

pdf("repeat_count_analysis.pdf", width = 10)

repeats <- data.frame()

for (file in list.files("repeat_counts/", pattern = "filtered")) {
  
  read_in <- read_delim(paste0("repeat_counts/", file), " ", col_names = c("number", "index"))
  read_in$sample <- file
  read_in_counts <- read_in %>% count(sample, number)
  repeats <- rbind(repeats, read_in_counts)
  rm(read_in, read_in_counts)
  
}

repeats <- repeats %>%
            separate(sample,
                     into = c("shRNA", "timePoint",
                              "replicate", "end",
                              "filtered", "paired",
                              "repeatType", "repeat",
                              "count")) %>%
            select(shRNA, timePoint,
                   replicate, end,
                   repeatType, number, n)


p <- repeats %>%
    group_by(shRNA, timePoint, repeatType, number) %>%
    summarise(n = sum(n)) %>%
  ggplot(aes(x = number, y = n, fill = repeatType)) +
    geom_col(position = "dodge") +
    facet_grid(shRNA ~ timePoint) +
    scale_y_log10()
print(p)



repeats <- data.frame()

for (file in list.files("repeat_counts/", pattern = "fq")) {
  
  read_in <- read_delim(paste0("repeat_counts/", file), " ", col_names = c("number", "index"))
  read_in$sample <- file
  read_in_counts <- read_in %>% count(sample, number)
  repeats <- rbind(repeats, read_in_counts)
  rm(read_in, read_in_counts)
  
}

repeats <- repeats %>%
  separate(sample,
           into = c("shRNA", "timePoint",
                    "replicate", "end",
                    "fq", "repeatType", "repeat",
                    "count")) %>%
  select(shRNA, timePoint,
         replicate, end,
         repeatType, number, n)


p <- repeats %>%
    group_by(shRNA, timePoint, repeatType, number) %>%
    summarise(n = sum(n)) %>%
  ggplot(aes(x = number, y = n, fill = repeatType)) +
    geom_col(position = "dodge") +
    facet_grid(shRNA ~ timePoint) +
    scale_y_log10()
print(p)

p <- repeats %>%
    group_by(shRNA, timePoint, repeatType, number) %>%
    summarise(n = sum(n)) %>%
  ggplot(aes(x = number, y = n, fill = shRNA)) +
    geom_col(position = "dodge") +
    facet_grid(repeatType ~ timePoint) +
    scale_y_log10()
print(p)

dev.off()