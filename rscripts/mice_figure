---
title: "mice_new_figures"
author: "Iza"
date: "2025-04-26"
output: html_document
---

import salmon results
```{r}

gsingle_m_filter2_mice<-gsingle_m_filter2 %>% filter(grad_fraction_max %in% c("3", "4"))

mean(gsingle_m_filter2_mice$tx_length)

gsingle_m_filter_human<-gsinglefinal %>% filter(grad_fraction_max %in% c("3", "4"))

median(gsingle_m_filter_human$tx_length)

gsingle_final2023_noRI<-gsinglefinal %>% filter(transcript_type != "retained_intron")

length(unique(gsingle_final2023_noRI$gene_name))

load("01_txs_filtered.RData")

load("Processed_res/mouse/gaussian_results.RData")

extract_total_mice<-txs_tpm_mean_mESC %>% filter(sample_type %in% c("total_RNA", "cell_extract"))

extract_total_mice_meantpm1<-extract_total_mice %>% filter(mean_tpm > 1)

length(unique(extract_total_mice_meantpm1$transcript_name))

library(tidyverse)
library(isoformic)
library(readr)
library(tximport)
library(DESeq2)
library(tximeta)

library(tibble)

# Your existing data
vector_samples <- c("JR2988","JR2989","JR2990","JR2991","JR2992","JR2993",
                    "JR2994","JR2995","JR2996","JR2997","JR2998","JR2999",
                    "JR3000","JR3001","JR3002","JR3003","JR3004","JR3005",
                    "JR3006","JR3007")

sampledata <- c("mR1F2","mR1F3","mR1F4","mR1F5","mR1F6","mR1F7","mR1F8",
                "mR1F9","mR1F10","mR1F11","mR2F2","mR2F3","mR2F4","mR2F5",
                "mR2F6","mR2F7","mR2F8","mR2F9","mR2F10","mR2F11")

# Create coldata table
coldata <- tibble(
  names = vector_samples,  # JR2988, JR2989, etc.
  cond = sampledata,    # mR1F2, mR1F3, etc.
  files = file.path(
    "/Users/imamede/My Drive/Projetos/RinnLab/terminus_res_mice/terminus_mice",
    paste0(vector_samples, "_quant"),
    "quant.sf"
  )
)

library(tximport)
library(readr)
library(purrr)

library(readr)
library(purrr)

# 1. Ler todos os arquivos e pegar só a coluna de IDs
all_names <- map(cleaned_quant, ~ {
  q <- read_tsv(.x, show_col_types = FALSE)
  q$Name
})

# 2. Descobrir quais IDs são comuns a todos
common_names <- reduce(all_names, intersect)

cat("Número de IDs comuns:", length(common_names), "\n")

# 3. Filtrar os arquivos para manter só esses IDs
cleaned_quant_common <- map_chr(cleaned_quant, ~ {
  q <- read_tsv(.x, show_col_types = FALSE)
  q <- q[q$Name %in% common_names, ]  # Filtra só os IDs comuns
  temp_file <- tempfile(fileext = ".sf")
  write_tsv(q, temp_file)
  temp_file
})

# 4. Agora faz o tximport!
library(tximport)

txi <- tximport(
  cleaned_quant_common,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

counts <- txi$counts

# Filtrar: manter apenas linhas onde o valor máximo em alguma coluna seja >=10
filtered_counts <- counts[rowSums(counts >= 10) > 0, ]

colnames(filtered_counts) <- colnames(txi$counts)

filtered_counts<-as.data.frame(filtered_counts)

colnames(filtered_counts) <-c("mR1F2","mR1F3","mR1F4","mR1F5","mR1F6","mR1F7","mR1F8",
                "mR1F9","mR1F10","mR1F11","mR2F2","mR2F3","mR2F4","mR2F5",
                "mR2F6","mR2F7","mR2F8","mR2F9","mR2F10","mR2F11")

```



```{r}
library(isoformic)
library(tidyverse)

#download_reference(version = "M25", file_type = "fasta", organism = "mouse")

fasta <- Biostrings::readDNAStringSet("data-raw/gencode.vM25.transcripts.fa.gz")

# Then extract transcript to gene mapping
tx2gene <- data.frame(
  tx_id = names(fasta),
  gene_id = sub(".*gene_id=([^;]+).*", "\\1", names(fasta))
)

df_separated <- tx2gene %>%
  separate_wider_delim(
    tx_id,
    delim = "|",
    names_sep = "_"
  )

df_separated <- df_separated %>%
  separate_wider_delim(
    gene_id,
    delim = "|",
    names_sep = "_"
  )

df_separated_names<-df_separated %>% dplyr::select(1,2,5,6,8)

names(df_separated_names) <-c("transcript_id", "gene_id", "transcript_name","gene_name", "transcript_type")

df_separated_names %>% select(gene_name, transcript_name)

gsingle_m_filter2 %>% left_join()

rlog_gauss_single_fits

```

do the barplot for all expressed
```{r}
gsingle_m_filter2 <- gsingle_m_filter2 %>%
  mutate(transcript_id = sub("\\..*", "", transcript_id))

df_separated_names<-df_separated_names %>%
  mutate(transcript_id = sub("\\..*", "", transcript_id))

gsingle_m_filter2_dic<-gsingle_m_filter2 %>% left_join(df_separated_names, by = "transcript_id")

df_toplot_BB<-df_separated_names %>% filter(transcript_id %in% rownames(filtered_counts))

library(dplyr)
library(ggplot2)

# 1. Group gene types AND reorder factor levels
gsingle_m_filter2_dic <- gsingle_m_filter2_dic %>%
  mutate(
    gene_type_grouped = case_when(
      gene_type %in% c("lincRNA", "protein_coding") ~ gene_type,
      TRUE ~ "other"
    ),
    # Force "other" to be the first level (will appear at the bottom)
    gene_type_grouped = factor(gene_type_grouped, 
                              levels = c("lincRNA", "protein_coding","other"))
  )

# 2. Plot with "other" at the bottom
gsingle_m_filter2_dic %>%
  mutate(All = "All") %>%
  ggplot() +
  geom_bar(
    aes(y = All, fill = gene_type_grouped),
    position = "fill",  # Stacks bars from bottom to top in factor level order
    #color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "other" = "gray70",          # Bottom layer
      "lincRNA" = "darkred",       # Middle layer
      "protein_coding" = "darkblue"  # Top layer
    ),
    name = "Gene Type"
  ) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = NULL)

cowplot::ggsave2("BB_mice_new_colors_GETs.pdf", height = 4, width = 2.7)
```

length versus peak fraction graph
```{r}

gsingle_m_filter2$grad_fraction_max<-as.character(gsingle_m_filter2$grad_fraction_max)

p<-gsingle_m_filter2 %>%
  #filter(gene_type != "lncRNA") %>%
  ggplot(aes(y=log10(tx_length), x=as.numeric(grad_fraction_max))) +
  geom_point(alpha=0.7, colour = "black", fill = "grey", stroke = 0.5, shape = 21) +
  labs(x="Gaussian Fraction Max",
       y="log10 length in bases") +
  theme_minimal() 

p_marg <- ggExtra::ggMarginal(p, data = gsinglefinalwlength, type="densigram")

p_marg

cowplot::ggsave2("Densigram_length_all_transcripts.pdf", height = 6, width = 6)
```

Now length comparing all the other transcripts

pvalue: 2.2e-16
log2fc: -0.430811

```{r}
tx_fasta <- Biostrings::readDNAStringSet("data-raw/gencode.vM25.transcripts.fa.gz")

# Extract lengths and match to your filtered_counts
transcript_lengths <- data.frame(
  transcript_id = names(tx_fasta),
  length = width(tx_fasta),
  stringsAsFactors = FALSE
) 

transcript_lengths <- transcript_lengths %>%
  mutate(transcript_id = sub("\\..*", "", transcript_id)) %>%
  filter(transcript_id %in% rownames(filtered_counts))


gsingl_forplot<-gsingle_m_filter2_dic %>% select(tx_length, transcript_id) %>%
  mutate(class = "GETs") 

gsingle_m_filter2_dic %>% select(tx_length)

transcript_lengths_allexpr<-transcript_lengths  %>%
  mutate(class = "All Expressed")

names(gsingl_forplot) <- c("length", "transcript_id", "class")

plotting_length<-rbind(gsingl_forplot, transcript_lengths_allexpr)

ggplot(plotting_length, aes(x = class, y = log10(length), fill = class)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside
  scale_fill_manual(values = c("GETs" = "#C1325A", "All Expressed" = "grey")) +
  labs(
    x = NULL,
    y = "Transcript Length log10 (bp)",
    #title = "Transcript Length Distribution Comparison"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",  # Remove legend (redundant with x-axis)
    plot.title = element_text(hjust = 0.5)
  )

cowplot::ggsave2("Violin_length_all_transcripts.pdf", height = 6, width = 6)

t_test_result_ot <- t.test(
  length ~ class, 
  data = plotting_length,
  var.equal = FALSE  # Default (Welch)
)


```

Violin size versus TPM per fraction
```{r}
long_counts <- filtered_counts %>%
  rownames_to_column("Transcript") %>%
  pivot_longer(
    cols = -Transcript,
    names_to = "Condition",
    values_to = "Count"
  )

long_counts <- long_counts %>%
  mutate(
    Fraction = sub("mR\\d+F(\\d+)", "\\1", Condition), 
    Replicate = sub("mR(\\d+)F\\d+", "\\1", Condition)  
  ) %>%
  mutate(Fraction = as.numeric(Fraction)) 

summary_counts <- long_counts %>%
  group_by(Transcript, Fraction) %>%
  summarize(
    Mean = mean(Count),
    Sum = sum(Count),
    Median = median(Count),
    .groups = "drop"
  )

summary_counts <- summary_counts %>%
  mutate(
    GET_Status = ifelse(Transcript %in% gsingle_m_filter2_dic$transcript_id, "GET", "All expressed")
  ) %>%
  filter(Mean > 1)

p<-summary_counts %>%
  #filter(Mean > 100) %>%
  ggplot(aes(x = as.numeric(Fraction), y = log10(Mean))) +
  geom_point(
    aes(col = GET_Status),
    position = position_jitter(width = 0.2, height = 0),
    size = 0,
    alpha = 0)+
  labs(x="Gaussian Fraction Max",
       y="log10 Mean TPM") +
  geom_violin(aes(group = interaction(Fraction, GET_Status), 
                  colour = GET_Status, fill = GET_Status), 
              alpha = 0.7)+
  theme_minimal() +
  scale_x_continuous(breaks = 2:11, n.breaks = 10, minor_breaks = NULL)+
  scale_color_manual(values = c("GET" = "#C1325A", "All expressed" = "grey"))+
  scale_fill_manual(values = c("GET" = "#C1325A", "All expressed" = "grey"))

p_marg <- ggExtra::ggMarginal(p, type = "densigram", groupColour = TRUE,
                      margins = "y")

p_marg

# Save plot
ggsave("Mice_Violin_sizevsTPM_ALL.pdf", plot = p_marg, height = 5, width = 6)
```

Now the violin with the same thing
pval:  0.6242 log2fc -0.08538392
```{r}
p<-summary_counts %>% ggplot(aes(x = GET_Status, y = log10(Sum), fill = GET_Status)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside
  scale_fill_manual(values = c("GET" = "#C1325A", "All Expressed" = "grey")) +
  labs(
    x = NULL,
    y = "Transcript Length log10 TPM",
    #title = "Transcript Length Distribution Comparison"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",  # Remove legend (redundant with x-axis)
    plot.title = element_text(hjust = 0.5)
  )
p
cowplot::ggsave2("Mice_Violin_TPM_ALL.pdf", plot = p, height = 5, width = 5)

summary_counts_getsonly<-summary_counts %>% filter(GET_Status == "GET") 

mean(summary_counts_getsonly$Sum)

t_test_result <- t.test(
  Sum ~ GET_Status, 
  data = summary_counts,
  var.equal = FALSE  # Default (Welch)
)


```

no mane but refseq selectMANE LIST for mice
43 form the mice list are the refseq select isoform 

```{r}
library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Fetch all RefSeq transcript IDs (NM_*, XM_*)
refseq_ensembl <- getBM(
  attributes = c(
    "refseq_mrna",        # RefSeq mRNA (NM_*)
    "refseq_ncrna",       # RefSeq non-coding RNA (NR_*, XR_*)
    "ensembl_transcript_id",  # ENSMUST*
    "external_gene_name"      # Gene symbol (e.g., "Pax6")
  ),
  mart = mart
)

refseq_filtered <- refseq_ensembl %>%
  filter(refseq_mrna != "" | refseq_ncrna != "")

intersect(summary_counts_getsonly$Transcript,
          refseq_filtered$ensembl_transcript_id)

refseq_filtered_toplot<-refseq_filtered %>% filter(external_gene_name %in% gsingle_m_filter2_dic$gene_name.x)

summary_counts_refseqpair<-summary_counts %>% filter(Transcript %in%refseq_filtered_toplot$ensembl_transcript_id) %>% mutate(new_class= "refseq_select")

summary_counts_GETpair<-summary_counts %>% filter(Transcript %in% gsingle_m_filter2_dic$transcript_id) %>% mutate(new_class = "GET")

to_plot_tpm<-rbind(summary_counts_refseqpair, summary_counts_GETpair)
```

Now the plot

pvalue: 2.2e-16, log2fc:2.232343
```{r}
ggplot(to_plot_tpm, aes(x = new_class, y = log10(Sum), fill = new_class)) +
  geom_violin() +
  scale_fill_manual(values = c("refseq_select" = "#8EB1C7", "GET" = "#27187E")) +
  labs( ) +
  theme_bw()
cowplot::ggsave2("Mice_Violin_GET_vs_refseqselect_tpm.pdf",height = 5, width = 5)

t_test_result <- t.test(
  Sum ~ new_class, 
  data = to_plot_tpm,
  var.equal = FALSE  # Default (Welch)
)
```

Now same for length

pval: 1.591e-13, log2fc: 0.5586474
```{r}
transcript_lengths_allexpr_refseqpair<-transcript_lengths_allexpr %>% filter(transcript_id %in%refseq_filtered_toplot$ensembl_transcript_id) %>% mutate(new_class= "refseq_select")

summary_counts_GETpair<-gsingl_forplot %>% filter(transcript_id %in% gsingle_m_filter2_dic$transcript_id & !transcript_id %in%refseq_filtered_toplot$ensembl_transcript_id) %>% mutate(new_class = "GET")

to_plot_length<-rbind(transcript_lengths_allexpr_refseqpair, summary_counts_GETpair)

ggplot(to_plot_length, aes(x = new_class, y = log10(length), fill = new_class)) +
  geom_violin() +
  scale_fill_manual(values = c("refseq_select" = "#8EB1C7", "GET" = "#27187E")) +
  labs( ) +
  theme_bw()
cowplot::ggsave2("Mice_Violin_GET_vs_refseqselect_legth.pdf",height = 5, width = 5)

t_test_result <- t.test(
  length ~ new_class, 
  data = to_plot_length,
  var.equal = FALSE  # Default (Welch)
)
```

BB all expressed transcript types
```{r}

# Define transcript types to keep (others will be grouped as "other")
keep_types <- c("lincRNA", "protein_coding", "nonsense_mediated_decay",
                "retained_intron", "processed_transcript")

gsingle_m_filter2_dic %>%
  mutate(
    All = "All",
    transcript_type_grouped = if_else(transcript_type %in% keep_types, transcript_type, "other"),
    transcript_type_grouped = factor(transcript_type_grouped, 
                                     levels = c("lincRNA", "protein_coding", 
                                                "retained_intron", "processed_transcript", 
                                                "nonsense_mediated_decay", "other"))
  ) %>%
  ggplot() +
  geom_bar(
    aes(y = All, fill = transcript_type_grouped),
    position = "fill",
    color = "grey"
  ) +
  scale_fill_manual(
    values = c(
      "nonsense_mediated_decay" = "#8FB2C7",
      "processed_transcript" = "#5C112A",
      "retained_intron" = "#C2335B",
      "other" = "gray70",          # Bottom layer
      "lincRNA" = "#282661",       # Middle layer
      "protein_coding" = "#2D2C79" 
    ),
    name = "Gene Type"
  ) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = NULL)

cowplot::ggsave2("BB_mice_new_colors_GETs_Transcripts.pdf", height = 4, width = 5)


filtered_counts_rn<-filtered_counts %>% rownames_to_column()

all_transcripts<-df_separated_names %>% filter(transcript_id %in%filtered_counts_rn$rowname)


all_transcripts %>%
    mutate(
    All = "All",
    transcript_type_grouped = if_else(transcript_type %in% keep_types, transcript_type, "other"),
    transcript_type_grouped = factor(transcript_type_grouped, 
                                     levels = c("lincRNA", "protein_coding", 
                                                "retained_intron", "processed_transcript", 
                                                "nonsense_mediated_decay", "other"))) %>%
  ggplot() +
  geom_bar(
    aes(y = All, fill = transcript_type_grouped),
    position = "fill",  # Stacks bars from bottom to top in factor level order
    color = "grey"
  ) +
  scale_fill_manual(
    values = c(
      "nonsense_mediated_decay" = "#8FB2C7",
      "processed_transcript" = "#5C112A",
      "retained_intron" = "#C2335B",
      "other" = "gray70",          # Bottom layer
      "lincRNA" = "#282661",       # Middle layer
      "protein_coding" = "#2D2C79"  # Top layer
    ),
    name = "Gene Type"
  ) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = NULL)

cowplot::ggsave2("BB_mice_new_colors_Allexpr_Transcripts.pdf", height = 4, width = 5)


```

Now some specific transcripts

```{r}

extract_total_mice

tx_to_tx<-df_separated_names %>% dplyr::select(transcript_id, transcript_name)

long_counts_wnames<-long_counts %>% left_join(tx_to_tx, by = c("Transcript"="transcript_id"))

selected_transcript <- "Akap8-203"
# Filter for that transcript
filtered_data <- long_counts_wnames %>%
  filter(transcript_name == selected_transcript)

# Calculate mean and standard deviation for each Fraction
summary_data <- filtered_data %>%
  group_by(Fraction) %>%
  summarize(
    mean_count = mean(Count, na.rm = TRUE),
    sd_count = sd(Count, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with error bars
ggplot(summary_data, aes(x = Fraction, y = mean_count)) +
  geom_line(color = "#27187E", size = 0.5) +
  geom_point(color = "#27187E", size = 1) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                width = 0.2, color = "#27187E") +
  theme_minimal() +
  labs(title = paste0(selected_transcript),
       x = "Gradient fraction",
       y = "TPM (mean)") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11))+
  theme_classic()

cowplot::ggsave2("mice_Akap8-203gradient.pdf", height = 3.5, width = 4)

```

Enrichment human and mice against GO BP
```{r}
gsingle_final2023<- readRDS("gsingle_final2023.rds")
library(enrichR)

# Run Enrichr analysis with GO Biological Process (GO BP)
enrich_result <- enrichr(gsingle_final2023$gene_name, "GO_Biological_Process_2018")

# View results
enrich_resulthuamn<-as.data.frame(enrich_result)


enrich_result <- enrichr(gsingle_final2023$gene_name, "GO_Biological_Process_2018")

# View results
enrich_resulthuamn<-as.data.frame(enrich_result)

enrich_result_mouse <- enrichr(gsingle_m$gene_name, "GO_Biological_Process_2018")

enrich_resulthuamn_filt<-enrich_resulthuamn %>% dplyr::top_n(20, GO_Biological_Process_2018.Combined.Score) %>% mutate(class = "GET_human")

enrich_resultmouse<-as.data.frame(enrich_result_mouse)

enrich_resultmouse_filt<-enrich_resultmouse %>% dplyr::top_n(20, GO_Biological_Process_2018.Combined.Score) %>% mutate(class = "GET_mouse")

result_enrichment<-rbind(enrich_resulthuamn_filt, enrich_resultmouse_filt)

library(ggplot2)
library(dplyr)

# Prepare data
result_enrichment <- result_enrichment %>%
  mutate(Term = GO_Biological_Process_2018.Term,
         Score = GO_Biological_Process_2018.Combined.Score,
         PValue = GO_Biological_Process_2018.P.value,
         logP = -log10(PValue))

# Reorder terms for better visualization
result_enrichment$Term <- factor(result_enrichment$Term,
                                 levels = rev(unique(result_enrichment$Term[order(result_enrichment$Score)])))

# Plot
ggplot(result_enrichment, aes(x = class, y = Term)) +
  geom_point(aes(size = logP, color =  Score)) +
  theme_minimal(base_size = 12) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(x = "Combined Score", y = "GO Biological Process") +
  theme_bw()+
  scale_color_gradient( low = "grey", high = "#C2335B")

cowplot::ggsave2("Dotplot_enrichmentres_GETmouseandhum.pdf", height = 5, width = 8.5)


```

Intersect mice and human
```{r}

load("orthologs.RData")


mouse_orthologs_filtered_merge
gene_name_human<-gsinglefinal 

gsingle_m_filter2

mouse_GETs_orthologs<-mouse_orthologs_filtered_merge %>% filter(mouse_gene_name %in% gsingle_m_filter2$gene_name)

GET_human_transcirptid_clean <- sub("\\..*$", "", gsinglefinal$transcript_id)

n_intersect_both<-intersect(GET_human_transcirptid_clean, mouse_GETs_orthologs$transcript_id_with_ortholog)

intersect_both_all<-mouse_GETs_orthologs %>% filter(transcript_id_with_ortholog %in%n_intersect_both) 

unique(intersect_both_all$gene_id_stable)

unique(gsingle_m_filter2$gene_name)


unique(gsinglefinal$gene_name)

intersect()

```

