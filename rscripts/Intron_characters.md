---
title: "Final_figures"
author: "Iza"
date: "2024-10-12"
output: html_document
---

```{r}
genelist_forvalidated_gauss<-readRDS("Processed_res/Isabela_processed_res/validatedRI_3reps_genelist.rds")

all_intron_chars<-readRDS("actino_all_intron_characters.rds")

all_intron_chars_truegausssel<-all_intron_chars %>% mutate(all_3rep_val = case_when(gene_id %in%genelist_forvalidated_gauss$gene_id ~ "true_gauss", TRUE ~ "not_true_gauss"))

table(all_intron_chars_truegausssel$all_3rep_val)

all_intron_chars_truegausssel_dpeth <- all_intron_chars_truegausssel %>% 
  select(-9, -10, -11, -12, -13, -14, -15, -16, -17, -18, -19)

all_intron_chars_truegausssel_dpeth <-distinct(all_intron_chars_truegausssel_dpeth)

test <-as.data.frame(table(all_intron_chars_truegausssel$gene_name))
```

intorn length
```{r}

#count poly P
group1_data <- all_intron_chars_truegausssel$count_polyP[all_intron_chars_truegausssel$all_3rep_val == "true_gauss"]
group2_data <- all_intron_chars_truegausssel$count_polyP[all_intron_chars_truegausssel$all_3rep_val == "not_true_gauss"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_polyP_tracks <- log2(group1_mean/group2_mean)
t_result_polyP_tracks <- t.test(count_polyP ~ all_3rep_val, data = all_intron_chars_truegausssel)
wilcox_result_polyP_tracks <- wilcox.test(count_polyP ~ all_3rep_val, data = all_intron_chars_truegausssel)
t_result_polyP_tracks_p<-t_result_polyP_tracks$p.value

#intron size
group1_data <- all_intron_chars_truegausssel$intron_size[all_intron_chars_truegausssel$all_3rep_val == "true_gauss"]
group2_data <- all_intron_chars_truegausssel$intron_size[all_intron_chars_truegausssel$all_3rep_val == "not_true_gauss"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_intronsize <- log2(group1_mean/group2_mean)
t_result_intronsize_tracks <- t.test(intron_size ~ all_3rep_val, data = all_intron_chars_truegausssel)
wilcox_result_intronsize_tracks <- wilcox.test(intron_size ~ all_3rep_val, data = all_intron_chars_truegausssel)
t_result_intronsize_tracks_p<-t_result_intronsize_tracks$p.value

#percentage CG

group1_data <- all_intron_chars_truegausssel$percentage_CG[all_intron_chars_truegausssel$all_3rep_val == "true_gauss"]
group2_data <- all_intron_chars_truegausssel$percentage_CG[all_intron_chars_truegausssel$all_3rep_val == "not_true_gauss"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_percentCG <- log2(group1_mean/group2_mean)
t_result_percentCG_tracks <- t.test(percentage_CG ~ all_3rep_val, data = all_intron_chars_truegausssel)
wilcox_result_percentCG_tracks <- wilcox.test(percentage_CG ~ all_3rep_val, data = all_intron_chars_truegausssel)

t_result_percentCG_tracks_p<-t_result_percentCG_tracks$p.value

# Intron depth
group1_data <- all_intron_chars_truegausssel$IntronDepth[all_intron_chars_truegausssel$all_3rep_val == "true_gauss"]
group2_data <- all_intron_chars_truegausssel$IntronDepth[all_intron_chars_truegausssel$all_3rep_val == "not_true_gauss"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_intorndepth <- log2(group1_mean/group2_mean)
t_result_IntronDepth_tracks <- t.test(IntronDepth ~ all_3rep_val, data = all_intron_chars_truegausssel)
wilcox_result_IntronDepth_tracks <- wilcox.test(IntronDepth ~ all_3rep_val, data = all_intron_chars_truegausssel)
t_result_IntronDepth_tracks_p<-t_result_IntronDepth_tracks$p.value

#IRatio
group1_data <- all_intron_chars_truegausssel$IRratio[all_intron_chars_truegausssel$all_3rep_val == "true_gauss"]
group2_data <- all_intron_chars_truegausssel$IRratio[all_intron_chars_truegausssel$all_3rep_val == "not_true_gauss"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_iratio <- log2(group1_mean/group2_mean)
t_result_IRratio_tracks <- t.test(IRratio ~ all_3rep_val, data = all_intron_chars_truegausssel)
wilcox_result_IRratio_tracks <- wilcox.test(IRratio ~ all_3rep_val, data = all_intron_chars_truegausssel)
t_result_IRratio_tracks_p<-t_result_IRratio_tracks$p.value
```

Make the table

```{r}
foldchanges_all <- c(fold_change_polyP_tracks, fold_change_intronsize,
                     fold_change_percentCG, fold_change_intorndepth, fold_change_iratio) #fold_change_spliceleft, fold_change_spliceright)

pvalues_all <- c(t_result_polyP_tracks_p, t_result_intronsize_tracks_p,
                     t_result_percentCG_tracks_p, t_result_IntronDepth_tracks_p, t_result_IRratio_tracks_p) #t_result_SpliceLeft_tracks_p, t_result_SpliceRight_tracks_p)

names_vector <- c("PolypyrimidineTract", "Intron_Length", "PercentageCG", "IntronDepth", "IRratio")#, "SpliceLeft", "SpliceRight")

table_forlollipop_actino<-data.frame(foldchanges_all, pvalues_all)

colnames(table_forlollipop_actino) <- c("log2FoldChange", "pvalue(ttest)")

rownames(table_forlollipop_actino) <- names_vector

```

```{r}
#count poly P
group1_data <- all_intron_chars_truegausssel$count_polyP[all_intron_chars_truegausssel$ov_to_gau == "TRUE"]
group2_data <- all_intron_chars_truegausssel$count_polyP[all_intron_chars_truegausssel$ov_to_gau == "FALSE"]


fold_change_polyP_tracks_hi <- log2(group1_mean/group2_mean)
t_result_polyP_tracks <- t.test(count_polyP ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_polyP_tracks <- wilcox.test(count_polyP ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_polyP_tracks_p_hi<-t_result_polyP_tracks$p.value


#intron size
group1_data <- all_intron_chars_truegausssel$intron_size[all_intron_chars_truegausssel$ov_to_gau == "TRUE"]
group2_data <- all_intron_chars_truegausssel$intron_size[all_intron_chars_truegausssel$ov_to_gau == "FALSE"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_intronsize_hi <- log2(group1_mean/group2_mean)
t_result_intronsize_tracks <- t.test(intron_size ~ ov_to_gau, data = all_intron_chars_truegausssel)
wilcox_result_intronsize_tracks <- wilcox.test(intron_size ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_intronsize_tracks_p_hi<-t_result_intronsize_tracks$p.value


#percentage CG
group1_data <- all_intron_chars_truegausssel$percentage_CG[all_intron_chars_truegausssel$ov_to_gau == "TRUE"]
group2_data <- all_intron_chars_truegausssel$percentage_CG[all_intron_chars_truegausssel$ov_to_gau == "FALSE"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_percentCG_hi <- log2(group1_mean/group2_mean)
t_result_percentCG_tracks <- t.test(percentage_CG ~ ov_to_gau, data = all_intron_chars_truegausssel)
wilcox_result_percentCG_tracks <- wilcox.test(percentage_CG ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_percentCG_tracks_p_hi<-t_result_percentCG_tracks$p.value

# Intron depth
group1_data <- all_intron_chars_truegausssel$IntronDepth[all_intron_chars_truegausssel$ov_to_gau == "TRUE"]
group2_data <- all_intron_chars_truegausssel$IntronDepth[all_intron_chars_truegausssel$ov_to_gau == "FALSE"]

group1_mean <- median(group1_data)
group2_mean <- median(group2_data)

fold_change_intorndepth_hi <- log2(group1_mean/group2_mean)
t_result_IntronDepth_tracks <- t.test(IntronDepth ~ ov_to_gau, data = all_intron_chars_truegausssel)
wilcox_result_IntronDepth_tracks <- wilcox.test(IntronDepth ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_IntronDepth_tracks_p_hi<-t_result_IntronDepth_tracks$p.value

#IRatio
group1_data <- all_intron_chars_truegausssel$IRratio[all_intron_chars_truegausssel$ov_to_gau == "TRUE"]
group2_data <- all_intron_chars_truegausssel$IRratio[all_intron_chars_truegausssel$ov_to_gau == "FALSE"]

fold_change_iratio_hi <- log2(group1_mean/group2_mean)
t_result_IRratio_tracks <- t.test(IRratio ~ ov_to_gau, data = all_intron_chars_truegausssel)
wilcox_result_IRratio_tracks <- wilcox.test(IRratio ~ ov_to_gau, data = all_intron_chars_truegausssel)
t_result_IRratio_tracks_p_hi<-t_result_IRratio_tracks$p.value


```

```{r}
foldchanges_all_hi <- c(fold_change_polyP_tracks_hi, fold_change_intronsize_hi,
                     fold_change_percentCG_hi, fold_change_intorndepth_hi, fold_change_iratio_hi)

pvalues_all_hi <- c(t_result_polyP_tracks_p_hi, t_result_intronsize_tracks_p_hi,
                     t_result_percentCG_tracks_p_hi, t_result_IntronDepth_tracks_p_hi, t_result_IRratio_tracks_p_hi)

names_vector <- c("PolypyrimidineTract", "Intron_Length", "PercentageCG", "IntronDepth", "IRratio")

table_forlollipop_hi<-data.frame(foldchanges_all_hi, pvalues_all_hi)

colnames(table_forlollipop_hi) <- c("log2FoldChange", "pvalue(ttest)")

rownames(table_forlollipop_hi) <- names_vector
```

Bind both tables

```{r}
table_forlollipop<-table_forlollipop_actino %>% mutate(class = "gauss_hi_conf") %>% rownames_to_column()
table_forlollipop_hi<-table_forlollipop_hi %>% mutate(class = "gaussian") %>% rownames_to_column()

table_forlollipop_hi_final<-rbind(table_forlollipop, table_forlollipop_hi)

#table_forlollipop_hi_final<- readRDS("gaussian_logfc_info.rds")

names(table_forlollipop_hi_final) <- c("cond", "log2FC", "pvalue", "group")

#with actino info


#write_rds(table_forlollipop_hi_final, "table_forlolipop_final.rds")

all_intron_chars_now <-readRDS("all_intron_chars_toplot.rds")

write.csv(all_intron_chars_now, "all_intron_chars_now.csv")

```
Plot

```{r}
all_intron_chars_toplot %>%
  #filter(!group %in% c("all_gauss")) %>%
  ggplot(aes(group, cond))+
  geom_point(aes(size = -log(pvalue), col = log2FC))+
  theme_minimal()+
  scale_color_gradient2(low="#0e1949",mid = "white", high="#923346")
  #facet_wrap(~group)

cowplot::ggsave2("IRfinder_results_plots/Lollipop_results_new_rT.pdf", height = 4, width = 4)

```


