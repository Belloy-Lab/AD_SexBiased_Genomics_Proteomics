#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(rrvgo)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(tm)
  library(wordcloud)
  library(stringr)
  library(RColorBrewer)
  library(simplifyEnrichment)
})

# Define command-line options
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--female_input", type = "character", help = "Female cluster file (.csv)"),
  make_option("--male_input", type = "character", help = "Male cluster file (.csv)")
)
opt <- parse_args(OptionParser(option_list = option_list))
out_dir <- opt$work_dir

# Function to clean GO descriptions with custom stopwords and preserve hyphenated terms
clean_text <- function(text) {
  text <- gsub("amyloid-beta", "amyloid_β", text, ignore.case = TRUE)
  text <- gsub("-", "_", text)
  corpus <- Corpus(VectorSource(text))
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub("[^[:alnum:]_]", " ", x))) # nolint: line_length_linter.
  custom_stopwords <- c("regulation", "positive", "negative", "production", "process", "activity",
                        "eif", "alpha", "via", "eif2", "class", "across", "resting") # nolint
  corpus <- tm_map(corpus, removeWords, c("and", "the", "of", "in", "to", custom_stopwords))
  corpus <- tm_map(corpus, stripWhitespace)
  corpus <- tm_map(corpus, content_transformer(function(x) gsub("_", "-", x)))
  sapply(corpus, as.character)
}

female_go_sim = fread(file.path(out_dir, opt$female_input)) %>% 
  mutate(log10_p = -log10(p.adjust),
         cluster_num = paste0("Cluster ", Cluster))

female_go_sim$`Cluster Name` <- str_wrap(female_go_sim$`Cluster Name`, width = 20)

fp1a = ggplot(female_go_sim, aes(x = "", y = log10_p)) +
  geom_boxplot(color = "#700000", fill = NA, outlier.shape = NA) +
  geom_point(color = "gray", size = 2, position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(. ~ cluster_num, scales = "free_x", space = "fixed") +
  coord_cartesian(ylim = c(1, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")) +
  labs(x = "Cluster Group", y = "-log10(p.adjust)")

ggsave(filename = file.path(out_dir, "Female_raw_cluster_pathways_box.jpg"), plot = fp1a, device = "jpg", height = 8, width = 14, units = "in", dpi = 320)

# Create Plot - violin
fp1b = ggplot(female_go_sim, aes(x = "", y = log10_p)) +
  geom_violin(color = "#700000", fill = NA, outlier.shape = NA) +
  geom_point(color = "gray", size = 2, position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(. ~ cluster_num, scales = "free_x", space = "fixed") +
  coord_cartesian(ylim = c(1, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")) +
  labs(x = "Cluster Group", y = "-log10(p.adjust)")

ggsave(filename = file.path(out_dir, "Female_raw_cluster_pathways_violin.jpg"), plot = fp1b, device = "jpg", height = 8, width = 14, units = "in", dpi = 320)

## Create Word Cloud for Each Cluster
cluster_groups <- split(female_go_sim, female_go_sim$`Cluster Name`)

for (cluster in names(cluster_groups)) {
  cluster_data <- cluster_groups[[cluster]]
  descriptions <- cluster_data$Description
  weights <- -log10(cluster_data$`Female p-value`)
  if (length(descriptions) == 0 || all(is.na(weights))) next
  cleaned_text <- clean_text(descriptions)
  tdm <- TermDocumentMatrix(Corpus(VectorSource(cleaned_text)))
  matrix <- as.matrix(tdm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  word_freq <- word_freq[word_freq > 0]
  if (length(word_freq) == 0) next
  terms <- names(word_freq)
  term_weights <- sapply(terms, function(term) {
    indices <- grep(term, cleaned_text, ignore.case = TRUE)
    if (length(indices) > 0) sum(weights[indices], na.rm = TRUE) else 0
  })
  output_file <- file.path(out_dir, paste0("female_raw_sim_", str_replace_all(cluster, "[^[:alnum:]]", "_"), "_pval_color.png"))
  png(output_file, width = 1600, height = 1200, res = 300)
  wordcloud(words = names(word_freq), freq = term_weights, min.freq = 0, max.words = 50, scale = c(2.5, 0.3), rot.per = 0.1, random.order = FALSE, colors = brewer.pal(8, "Dark2"))
  dev.off()
}

# Create a word cloud for each cluster - Using p-value weights - black
for (cluster in names(cluster_groups)) {
  cluster_data <- cluster_groups[[cluster]]
  descriptions <- cluster_data$Description
  weights <- -log10(cluster_data$`Female p-value`)
  if (length(descriptions) == 0 || all(is.na(weights))) next
  cleaned_text <- clean_text(descriptions)
  tdm <- TermDocumentMatrix(Corpus(VectorSource(cleaned_text)))
  matrix <- as.matrix(tdm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  word_freq <- word_freq[word_freq > 0]
  if (length(word_freq) == 0) next
  terms <- names(word_freq)
  term_weights <- sapply(terms, function(term) {
    indices <- grep(term, cleaned_text, ignore.case = TRUE)
    if (length(indices) > 0) sum(weights[indices], na.rm = TRUE) else 0
  })
  output_file <- file.path(out_dir, paste0("female_raw_sim_", str_replace_all(cluster, "[^[:alnum:]]", "_"), "_pval_bw.png"))
  png(output_file, width = 1600, height = 1200, res = 300)
  wordcloud(words = names(word_freq), freq = term_weights, min.freq = 0, max.words = 50, scale = c(2.5, 0.3), rot.per = 0.1, random.order = FALSE, colors = "black")
  dev.off()
}

rm(list = setdiff(ls(), c("opt", "out_dir", "clean_text")))

## Visualization for Raw Male - Clustering Approach 1
male_go_sim = fread(file.path(out_dir, opt$male_input)) %>% 
  mutate(log10_p = -log10(p.adjust),
         cluster_num = paste0("Cluster ", Cluster))

male_go_sim$`Cluster Name` <- str_wrap(male_go_sim$`Cluster Name`, width = 20)

mp1a = ggplot(male_go_sim, aes(x = "", y = log10_p)) +
  geom_boxplot(color = "#3780b5", fill = NA, outlier.shape = NA) +
  geom_point(color = "gray", size = 2, position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(. ~ cluster_num, scales = "free_x", space = "fixed") +
  coord_cartesian(ylim = c(1, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")) +
  labs(x = "Cluster Group", y = "-log10(p.adjust)")

ggsave(filename = file.path(out_dir, "Male_raw_cluster_pathways_box.jpg"), plot = mp1a, device = "jpg", height = 8, width = 12, units = "in", dpi = 320)

# Create Plot - violin
mp1b = ggplot(male_go_sim, aes(x = "", y = log10_p)) +
  geom_violin(color = "#3780b5", fill = NA, outlier.shape = NA) +
  geom_point(color = "gray", size = 2, position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(. ~ cluster_num, scales = "free_x", space = "fixed") +
  coord_cartesian(ylim = c(1, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")) +
  labs(x = "Cluster Group", y = "-log10(p.adjust)")

ggsave(filename = file.path(out_dir, "Male_raw_cluster_pathways_violin.jpg"), plot = mp1b, device = "jpg", height = 8, width = 12, units = "in", dpi = 320)

## Create Word Cloud for Each Cluster
cluster_groups <- split(male_go_sim, male_go_sim$`Cluster Name`)

# Create a word cloud for each cluster - Using p-value weights - colored
for (cluster in names(cluster_groups)) {
  cluster_data <- cluster_groups[[cluster]]
  descriptions <- cluster_data$Description
  weights <- -log10(cluster_data$`Male p-value`)
  if (length(descriptions) == 0 || all(is.na(weights))) {
    cat("No valid data to plot for cluster:", cluster, "\n")
    next
  }
  cleaned_text <- clean_text(descriptions)
  tdm <- TermDocumentMatrix(Corpus(VectorSource(cleaned_text)))
  matrix <- as.matrix(tdm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  word_freq <- word_freq[word_freq > 0]
  if (length(word_freq) == 0) {
    cat("No words to plot for cluster:", cluster, "\n")
    next
  }
  terms <- names(word_freq)
  term_weights <- sapply(terms, function(term) {
     indices <- grep(term, cleaned_text, ignore.case = TRUE)
    if (length(indices) > 0) sum(weights[indices], na.rm = TRUE) else 0
  })
  output_file <- file.path(out_dir, paste0("male_raw_sim_", str_replace_all(cluster, "[^[:alnum:]]", "_"), "_pval_color.png"))
  png(output_file, width = 1600, height = 1200, res = 300)
  cat("Generating word cloud for cluster:", cluster, "with", length(terms), "words\n")
  wordcloud(words = names(word_freq),freq = term_weights,min.freq = 0,max.words = 50,scale = c(2.5, 0.3),rot.per = 0.1,random.order = FALSE,colors = brewer.pal(8, "Dark2"))
  dev.off()
}

# Create a word cloud for each cluster - Using p-value weights - black
for (cluster in names(cluster_groups)) {
  cluster_data <- cluster_groups[[cluster]]
  descriptions <- cluster_data$Description
  weights <- -log10(cluster_data$`Male p-value`)
  if (length(descriptions) == 0 || all(is.na(weights))) {
    cat("No valid data to plot for cluster:", cluster, "\n")
    next
  }
  cleaned_text <- clean_text(descriptions)
  tdm <- TermDocumentMatrix(Corpus(VectorSource(cleaned_text)))
  matrix <- as.matrix(tdm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  word_freq <- word_freq[word_freq > 0]
  if (length(word_freq) == 0) {
    cat("No words to plot for cluster:", cluster, "\n")
    next
  }
  terms <- names(word_freq)
  term_weights <- sapply(terms, function(term) {
    indices <- grep(term, cleaned_text, ignore.case = TRUE)
    if (length(indices) > 0) sum(weights[indices], na.rm = TRUE) else 0
  })
  output_file <- file.path(out_dir, paste0("male_raw_sim_", str_replace_all(cluster, "[^[:alnum:]]", "_"), "_pval_bw.png"))
  png(output_file, width = 1600, height = 1200, res = 300)
  cat("Generating word cloud for cluster:", cluster, "with", length(terms), "words\n")
  wordcloud(words = names(word_freq), freq = term_weights, min.freq = 0, max.words = 50, scale = c(2.5, 0.3), rot.per = 0.1, random.order = FALSE, colors = "black")
  dev.off()
}