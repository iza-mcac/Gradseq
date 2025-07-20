---
title: "eCLIP"
author: "Iza"
date: "2025-05-13"
output: html_document
---

```{r}
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyverse)


metadata_2 <- read.delim("eCLIP_encode_data/experiment_report_2025_5_13_22h_42m.tsv", 
                      skip = 1, header = TRUE, 
                      na.strings = c("", "NA"), 
                      quote = "")

# Get list of BED files
setwd("eCLIP_encode_data/")
bed_files <- list.files(pattern = "\\.bed\\.gz$")
file_ids <- gsub("\\.bed\\.gz$", "", bed_files)

```

Process the metadata

```{r}
processed_metadata <- metadata_2 %>%
  # Select relevant columns
  select(Accession, Target = Target.gene.symbol, 
         CellLine = Biosample.term.name, Files) %>%
  # Extract all ENCFF IDs from the Files column
  mutate(ENCFF_IDs = str_extract_all(Files, "ENCFF[0-9A-Z]+")) %>%
  # Create one row per ENCFF ID
  unnest(ENCFF_IDs) %>%
  # Filter to only include files we actually have
  filter(ENCFF_IDs %in% file_ids) %>%
  # Add the full filename
  mutate(filename = paste0(ENCFF_IDs, ".bed.gz"))

write.csv(processed_metadata, "file_metadata_mapping.csv", row.names = FALSE)
```

Import .bed and add metadata
```{r}
library(GenomicRanges)
encode_bed_cols <- c("chrom", "start", "end", "name", "score", "strand", 
                    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")


bed_files <- list.files("eCLIP_encode_data", pattern = "\\.bed\\.gz$", full.names = TRUE)

all_granges <- list()

for (file in bed_files) {
  # Read the BED file
  bed_data <- read.table(file, sep = "\t", header = FALSE,
                        col.names = encode_bed_cols[1:count.fields(file)[1]],
                        fill = TRUE)
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = bed_data$chrom,
    ranges = IRanges(start = bed_data$start + 1, end = bed_data$end), # BED is 0-based
    strand = bed_data$strand,
    name = bed_data$name,
    score = bed_data$score
  )
  
  # Add metadata column with filename
  mcols(gr)$source_file <- basename(file)
  
  # Store in list
  all_granges[[basename(file)]] <- gr
}

# Combine all GRanges objects
combined_gr <- do.call(c, unname(all_granges))

# Convert to dataframe if needed
combined_df <- as.data.frame(combined_gr)

write_rds(combined_df, "combined_df.rds")
```
Processign I already did no need to re-do
```{r}
combined_df<-readRDS("combined_df.rds")


combined_df %>% 
  mutate(num_parts = lengths(strsplit(name, "_"))) %>% 
  count(num_parts)

combined_df_separated <- combined_df %>%
  separate(name, into = c("part1", "part2", "part3"),
           sep = "_", remove = TRUE)


combined_df_separated <- combined_df %>%
  separate(name, into = c("part1", "part2", "part3"), 
           sep = "_", remove = TRUE, fill = "left")

library(data.table)

# Convert to data.table (modifies by reference, low memory overhead)
setDT(combined_df)

# Split the column efficiently
combined_df[, c("part1", "part2", "part3", "part4") := tstrsplit(name, "_", fixed = TRUE, fill = "")]

combined_df_whatiwant<-combined_df %>% dplyr::filter(part2 %in% c("HepG2","K562"))

```

Now the stuff really
```{r}
combined_df<-readRDS("combined_df.rds")

geenlist<-readRDS("Processed_res/Isabela_processed_res/")

gau_true_change_actino_logFC<-readRDS("Processed_res/Isabela_processed_res/gau_true_change_actino_logFC.rds")

head(combined_df)

```

Import anno
```{r}
# Your table of interest
introns_of_interest <- gau_true_change_actino_logFC %>%
  select(seqnames, start, end, strand, gene_name, gene_id, log2FC)

# Load GFF file (adjust path as needed)
gff <- rtracklayer::import("Ref/gencode.v38.annotation.gff3.gz")
gene_gff <- as.data.frame(gff) %>% filter(type == "gene")


library(GenomicRanges)
library(dplyr)
library(rtracklayer)


gff <- rtracklayer::import("Ref/gencode.v38.annotation.gff3.gz")

exons_df <- as.data.frame(gff) %>%
  filter(type == "exon") %>%
  mutate(
    transcript_id = ifelse(
      !is.na(transcript_id), transcript_id,
      ifelse(
        grepl("transcript_id=", attributes),
        sub(".*transcript_id=([^;]+);.*", "\\1", attributes),
        NA
      )
    )
  ) %>%
  select(seqnames, start, end, strand, gene_id, gene_name, transcript_id, exon_number) %>%
  arrange(gene_id, transcript_id, as.numeric(exon_number))

calculate_introns <- function(exon_group) {
  exon_group <- exon_group %>%
    mutate(exon_number = as.numeric(exon_number)) %>%
    arrange(exon_number)
  
  if (nrow(exon_group) < 2) return(tibble())
  
  intron_starts <- exon_group$end[-nrow(exon_group)] + 1
  intron_ends <- exon_group$start[-1] - 1
  valid <- intron_starts < intron_ends
  
  if (sum(valid) == 0) return(tibble())
  
  n <- sum(valid)
  
  tibble(
    seqnames = rep(first(exon_group$seqnames), n),
    start = intron_starts[valid],
    end = intron_ends[valid],
    strand = rep(first(exon_group$strand), n),
    gene_id = rep(first(exon_group$gene_id), n),
    gene_name = rep(first(exon_group$gene_name), n),
    transcript_id = rep(first(exon_group$transcript_id), n),
    intron_number = seq_len(n)
  )
}


all_introns_df <- exons_df %>%
  group_by(transcript_id) %>%
  group_modify(~ calculate_introns(.x), .progress = TRUE) %>%
  ungroup()


same_genes_introns <- all_introns_df %>%
  filter(gene_id %in% unique(gau_true_change_actino_logFC$gene_id))

your_introns_gr <- makeGRangesFromDataFrame(
  gau_true_change_actino_logFC,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

other_introns_gr <- makeGRangesFromDataFrame(
  same_genes_introns %>%
    anti_join(as.data.frame(your_introns_gr), 
  by = c("seqnames", "start", "end", "strand")),
  keep.extra.columns = TRUE
)

protein_binding_gr <- makeGRangesFromDataFrame(
  combined_df,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

```

Intersect the results and maybe stats

```{r}
hits_your_introns <- findOverlaps(protein_binding_gr, your_introns_gr)
hits_other_introns <- findOverlaps(protein_binding_gr, other_introns_gr)

summary_df <- data.frame(
  Category = c("Only GETs", "Only other introns", "Both", "Neither"),
  Count = c(
    length(setdiff(unique(queryHits(hits_your_introns)), unique(queryHits(hits_other_introns)))),
    length(setdiff(unique(queryHits(hits_other_introns)), unique(queryHits(hits_your_introns)))),
    length(intersect(unique(queryHits(hits_your_introns)), unique(queryHits(hits_other_introns)))),
    length(protein_binding_gr) - length(union(unique(queryHits(hits_your_introns)), unique(queryHits(hits_other_introns))))
  )
)
print(summary_df)

```

Fisher exact
```{r}
# Count overlaps
your_introns_df <- as.data.frame(your_introns_gr)
your_introns_df$binding_count <- countOverlaps(your_introns_gr, protein_binding_gr)

other_introns_df <- as.data.frame(other_introns_gr)
other_introns_df$binding_count <- countOverlaps(other_introns_gr, protein_binding_gr)


library(ggplot2)

your_introns_df_GEtrihiconf<-your_introns_df %>% dplyr::filter(gauRI_high_conf == "yes" & status == "up")

ggplot() +
  geom_density(data = your_introns_df_GEtrihiconf, aes(x = binding_count), fill = "grey", alpha = 0.7) +
  geom_density(data = other_introns_df, aes(x = binding_count), fill = "#923346", alpha = 0.7) +
  theme_minimal()

cowplot::ggsave2("density_gauRIhiconf_vs_other.pdf", height = 4, width = 7)

```

Are RIs enrihce din biding?
```{r}
# How many RIs have ≥1 binding
k <- sum(your_introns_df$binding_count > 0)

# Total RIs
m <- nrow(your_introns_df)

# Background: other introns
n <- nrow(other_introns_df)

# How many other introns have ≥1 binding
x <- sum(other_introns_df$binding_count > 0)

# Hypergeometric test: probability of ≥k binding RIs by chance
phyper(q = k - 1, m = m + n, n = n, k = m, lower.tail = FALSE)


```

Exclusive events
```{r}
# Hits
hits_your <- findOverlaps(protein_binding_gr, your_introns_gr)
hits_other <- findOverlaps(protein_binding_gr, other_introns_gr)

# Which peaks are exclusive to your introns
your_peak_ids <- unique(queryHits(hits_your))
other_peak_ids <- unique(queryHits(hits_other))

exclusive_to_your <- setdiff(your_peak_ids, other_peak_ids)
exclusive_to_other <- setdiff(other_peak_ids, your_peak_ids)
shared <- intersect(your_peak_ids, other_peak_ids)

length(exclusive_to_your)  # count
length(shared)
length(exclusive_to_other)

```

Random forest
```{r}
library(randomForest)

your_introns_df$type <- "RI"
other_introns_df$type <- "Other"


missing_in_other <- setdiff(names(your_introns_df), names(other_introns_df))
missing_in_your  <- setdiff(names(other_introns_df), names(your_introns_df))


for (col in missing_in_other) {
  other_introns_df[[col]] <- NA
}
for (col in missing_in_your) {
  your_introns_df[[col]] <- NA
}


other_introns_df <- other_introns_df[, names(your_introns_df)]


your_introns_df$type <- "RI"
other_introns_df$type <- "Other"


full_df <- rbind(your_introns_df, other_introns_df)
full_df$length <- full_df$end - full_df$start + 1

full_df$binding <- as.factor(full_df$binding_count > 0)


model_df <- full_df %>%
  select(binding, length, strand, gene_name) %>%
  mutate(across(c(strand, gene_name), as.factor)) %>%
  na.omit()

set.seed(42)
rf_model <- randomForest(binding ~ ., data = model_df, importance = TRUE)


print(rf_model)
varImpPlot(rf_model)

```

```{r}
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(pheatmap)

# Your introns
hits_your <- findOverlaps(protein_binding_gr, your_introns_gr)
df_your <- as.data.frame(hits_your)
df_your$intron_type <- "RI"
df_your$intron_id <- paste0("RI_", df_your$subjectHits)

# Other introns
hits_other <- findOverlaps(protein_binding_gr, other_introns_gr)
df_other <- as.data.frame(hits_other)
df_other$intron_type <- "Other"
df_other$intron_id <- paste0("Other_", df_other$subjectHits)

# Combine
df_all <- bind_rows(df_your, df_other)

# Add protein info and score
df_all$protein <- protein_binding_gr$part1[df_all$queryHits]
df_all$score <- protein_binding_gr$score[df_all$queryHits]

```

protein thing
```{r}
# Pivot to wide format
heatmap_mat <- df_all %>%
  group_by(protein, intron_id) %>%
  summarise(total_score = sum(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = intron_id, values_from = total_score, values_fill = 0)

# Convert to matrix (rownames = protein)
heatmap_matrix <- as.matrix(heatmap_mat[,-1])
rownames(heatmap_matrix) <- heatmap_mat$protein


```

Ploting the rf pretty
```{r}
library(ggplot2)
library(randomForest)
library(dplyr)
library(tibble)

# Extract importance measures
imp_df <- as.data.frame(importance(rf_model)) %>% 
  tibble::rownames_to_column("RBP") %>%
  arrange(desc(MeanDecreaseAccuracy))

# Select top 10 RBPs (by MeanDecreaseAccuracy)
top_rbps <- imp_df %>% 
  slice_max(MeanDecreaseGini, n = 15) %>% 
  pull(RBP)

# Create plot
ggplot(imp_df, aes(x = MeanDecreaseGini, y = MeanDecreaseAccuracy)) +
  geom_point(color = ifelse(imp_df$RBP %in% top_rbps, "#C1325A", "grey70"),
             size = ifelse(imp_df$RBP %in% top_rbps, 3, 2),
             alpha = 0.8) +
  geom_text_repel(data = subset(imp_df, RBP %in% top_rbps),
                  aes(label = RBP),
                  color = "#C1325A",
                  size = 3.5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(title = "Random Forest Variable Importance",
       x = "Mean Decrease in Gini Index",
       y = "Mean Decrease in Accuracy") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(
    breaks = 1:14,
    limits = c(1, 14)
  )

cowplot::ggsave2("Rf_results_markedtop15decreaseGini.pdf", height = 6, width = 8)

```

