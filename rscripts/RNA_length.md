---
title: "Malat1_Graph"
author: "Iza"
date: "2024-06-15"
output: html_document
---

```{r}
gsinglefinal<-readRDS("gsingle_final2023.rds")

counts_ext_new<-readRDS("counts_ext_new.rds")

gaussian_f11<-gaussian_counts_human %>% filter(Condition == "f11")
gsinglefinal 

gaussianf11_lncrna<-gaussian_f11 %>% filter(transcript_type == "lncRNA")

gaussianf11_lncrna_tb<- as.data.frame(table(gaussianf11_lncrna$parent_gene))

gsingle_dic<-gsingle %>% left_join(txtogene_v38, by = "transcript_id")

gsingle_11<-gsingle_dic %>% filter(grad_fraction_max == "11")

gsingle_10<-gsingle_dic %>% filter(grad_fraction_max == "10")

gsingle_9<-gsingle_dic %>% filter(grad_fraction_max == "9")

```

Try to make a length versus size thing
```{r}
gaussian_finalf10<-gsinglefinal %>% filter(grad_fraction_max == "10" & transcript_type == "lncRNA")

gaussian_finalf9<-gsinglefinalwlength %>% filter(grad_fraction_max == "9" & transcript_type == "lncRNA")

gsinglefinalwlength<-gsinglefinal %>% left_join(length_df, by = c("transcript_id"="tx_id"))

gaussian_final_onlylnc<-gsinglefinalwlength %>% filter(transcript_type == "lncRNA")

gaussian_final_notlnc<-gsinglefinalwlength %>% filter(transcript_type != "lncRNA")

gsinglefinalwlength_forcorr<-gsinglefinalwlength %>% dplyr::select(1, grad_fraction_max, length)

p<-gsinglefinalwlength %>%
  filter(gene_type != "lncRNA") %>%
  ggplot(aes(y=log10(length), x=as.numeric(grad_fraction_max))) +
  geom_point(alpha=0.7, colour = "black", fill = "grey", stroke = 0.5, shape = 21) +
  labs(x="Gaussian Fraction Max",
       y="log10 length in bases") +
  theme_minimal() 

p_marg <- ggMarginal(p, data = gsinglefinalwlength, type="densigram")

```


looking for noncoding on fraction 11 and 12

```{r}

txtogene_v38<-isoformic::make_tx_to_gene("/Users/imamede/My Drive/Projetos/RinnLab/Grad-seq/Ref/gencode.v38.transcripts.fa.gz", file_type = "fasta")

txtotx_v38<-txtogene_v38 %>% select(1,5,6,8)

counts_all <-readRDS("/Users/imamede/My Drive/Projetos/RinnLab/Grad-seq/Processed_res/Isabela_processed_res/Human/salmon.merged.transcript_counts.rds")

counts_all$names

counts_gauss_human<-as.data.frame(counts_all@assays@data$counts)

counts_gauss_human_rn<-rownames_to_column(counts_gauss_human)

counts_gauss_human_rn_dic<-counts_gauss_human_rn %>% left_join(txtotx_v38 , by = c("rowname"="transcript_id"))

counts_gauss_human_rn_dic_11<-counts_gauss_human_rn_dic %>% select(matches("11"), 64, 65, 66)

counts_gauss_human_rn_dic_11<-counts_gauss_human_rn_dic_11 %>% 
  mutate(ext11_sum =rowSums(select(.,contains("ext"))),
         tot11_sum =rowSums(select(.,contains("tot"))))

counts_gauss_human_rn_dic_11<-counts_gauss_human_rn_dic_11 %>%filter(transcript_type =="lncRNA")

counts_gauss_human_rn_dic_11_0.5<- counts_gauss_human_rn_dic_11[counts_gauss_human_rn_dic_11$tot11_sum < 0.5 * counts_gauss_human_rn_dic_11$ext11_sum, ] 

counts_gauss_human_rn_dic_11_0.25_high300_heat<-counts_gauss_human_rn_dic_11_0.25_high300 %>% select(9,10)

counts_gauss_human_rn_dic_11_0.25_high300_gauss<-counts_gauss_human_rn_dic_11_0.25_high300 %>% filter(transcript_name %in% gsinglefinal$transcript_name)

```

plotting malat1 isoforms

```{r}
counts_all 

counts_gauss_human_rn_dic_tot<-counts_gauss_human_rn_dic %>% select(matches("tot"), 64, 65, 66) %>% mutate(exp = "total")

counts_gauss_human_rn_dic_ext<-counts_gauss_human_rn_dic %>% select(matches("ext"), 64, 65, 66) %>% mutate(exp = "extract")

  counts_human_binded<-rbind(counts_gauss_human_rn_dic_tot, counts_gauss_human_rn_dic_ext)

```

plot the BB with all the filters and correct colors

```{r}
gsingle_dic<-gsingle %>% left_join(txtogene_v38, by = "transcript_id")

gsinglefilter3_dic<-gsingle_filter3 %>% left_join(txtogene_v38, by = "transcript_id")

# Desired order for the transcript types
desired_order <- c("lncRNA", "protein_coding", "retained_intron", 
                   "processed_transcript", "nonsense_mediated_decay",
                   "Gene", "processed_pseudogene", "transcribed_unprocessed_pseudogene")

# Convert transcript_type to a factor with the specified levels
gsinglefilter3_dic <- gsinglefilter3_dic %>%
  mutate(transcript_type = factor(transcript_type, levels = desired_order))

# Define the color names for the biotypes
tx_type_color_names <- c("#933447", "#171b48", "#33518d", "#7295c1",
                         "#620f1e", "#b3b3b3", "#b3b3b3", "#b3b3b3")

# Assign names to the color vector
names(tx_type_color_names) <- c("lncRNA", "protein_coding", "retained_intron", 
                                "processed_transcript", "nonsense_mediated_decay",
                                "Gene", "processed_pseudogene", "transcribed_unprocessed_pseudogene")

# Create the bar plot
gsinglefilter3_dic %>%
  mutate(All = "All") %>%
  ggplot() +
  geom_bar(aes(y = All, fill = transcript_type), position = "fill") +
  scale_fill_manual(values = tx_type_color_names) +
  theme_bw() +
  coord_flip()

cowplot::ggsave2("BB_gsinglefilter3.pdf", width = 4.5, height = 6)

```

Correlation
Table done by script: drafttogetlength e this oen untill line 40
```{r}
forcorr_length<-gsinglefinalwlength_forcorr %>% column_to_rownames( var = "transcript_id")

human_countsss_res <-Hmisc::rcorr(as.matrix(forcorr_length),
                                    type = "pearson")

r_value <- round(human_countsss_res$r[1, 2], 2)  # Correlation coefficient
p_value <- signif(human_countsss_res$P[1, 2], 3)

human_countsss_res$P


summary_counts <- summary_counts %>% column_to_rownames(var = "Transcript")

# Perform statistical test and calculate fold change
group1_data <- summary_counts$Mean[summary_counts$GEt_Status == "GET"]
group2_data <- summary_counts$Mean[summary_counts$GEt_Status == "notGET"]

group1_mean<- mean(group1_data)

group2_mean <- mean(group2_data)

fold_change_blastn_gau <- log2(group1_mean/group2_mean)
t_result_blastn_tracks <- t.test(Mean ~ GEt_Status, data = summary_counts)
t_result_blastn_tracks_p_gau<-t_result_blastn_tracks$p.value


```


