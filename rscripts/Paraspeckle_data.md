---
title: "Nuclear_Speckle_pd"
author: "Iza"
date: "2024-10-21"
output: html_document
---

Packages
```{r}
library(tidyverse)
library(ggplot2)
```

Our data
```{r, warning = FALSE, results = "hide", message=FALSE}
val_RIgau_3reps<-readRDS("Processed_res/Isabela_processed_res/val_RIgau_3reps.rds")

gsinglefinal<-readRDS("gsingle_final2023.rds")
```

Data Barutcu Blencowe sup tables sep by IP
```{r, message=FALSE, warning = FALSE, results = "hide"}
library(readxl)

SRSF1_targ<- read_excel("Ref/ref_blencowe/SRSF1_targets.xlsx")
SRSF7_targ<- read_excel("Ref/ref_blencowe/SRSF7_targets.xlsx")
RPS1_targ<- read_excel("Ref/ref_blencowe/RNPS1_targets.xlsx")
LMNA_targ<- read_excel("Ref/ref_blencowe/LMNA_targets.xlsx")
PML_targ<- read_excel("Ref/ref_blencowe/PML_targets.xlsx")
SP100_targ<- read_excel("Ref/ref_blencowe/SP100_targets.xlsx")
SMN2_targ<- read_excel("Ref/ref_blencowe/SMN2_targets.xlsx")
WRAP53_targ<- read_excel("Ref/ref_blencowe/WRAP53_targets.xlsx")
SAM68_targ<- read_excel("Ref/ref_blencowe/SAM58_targets.xlsx")
FBL_targ<- read_excel("Ref/ref_blencowe/FBL_targets.xlsx")
NPAT_targ<- read_excel("Ref/ref_blencowe/NPAT_targets.xlsx")

```

Test the first three (speckle markers) and us ethe IP versus Input ratio and mark the gaussians for ridge
Using the 909 genes that have gaussian transcripts

```{r}
SRSF1_targ_ip_input <- SRSF1_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(SRSF1_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(SRSF1_targ[["...11"]])) %>% mutate(factor = "SRSF1")

SRSF7_targ_ip_input <-SRSF7_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(SRSF7_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))%>% 
  dplyr::filter(!is.na(SRSF7_targ[["...11"]])) %>% mutate(factor = "SRSF7")


RPS1_targ_ip_input<-RPS1_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(RPS1_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(RPS1_targ[["...11"]])) %>% mutate(factor = "RPS1")

# combine to get a ratio of gauss to not

speckle_ip<-rbind(SRSF1_targ_ip_input, SRSF7_targ_ip_input, RPS1_targ_ip_input)

speckle_ip_pos<-speckle_ip %>% filter(`pulldown vs input log2FoldChange residuals` > 1)

speckle_ip_pos

table(speckle_ip_pos[[3]], speckle_ip_pos[[4]])


```

Now three non-Speckle stuff

```{r}
LMNA_targ_ip_input<-LMNA_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(RPS1_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(RPS1_targ[["...11"]])) %>% mutate(factor = "LMNA")

PML_targ_ip_input<-PML_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(PML_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(PML_targ[["...11"]])) %>% mutate(factor = "PML")

SP100_targ_ip_input<-SP100_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(SP100_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(SP100_targ[["...11"]])) %>% mutate(factor = "SP100")

SMN2_targ_ip_input<-SMN2_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(SMN2_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(SMN2_targ[["...11"]])) %>% mutate(factor = "SMN2")

WRAP53_targ_ip_input<-WRAP53_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(WRAP53_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(WRAP53_targ[["...11"]])) %>% mutate(factor = "WRAP53")

SAM68_targ_ip_input<-SAM68_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(SAM68_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(SAM68_targ[["...11"]])) %>% mutate(factor = "SAM68")

FBL_targ_ip_input<-FBL_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(FBL_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(FBL_targ[["...11"]])) %>% mutate(factor = "FBL")

NPAT_targ_ip_input<-NPAT_targ %>% 
  dplyr::select(11, 12) %>% 
  mutate(g_status = ifelse(NPAT_targ[["...11"]] %in% gsinglefinal$gene_name, "gaussian", "no"))  %>% 
  dplyr::filter(!is.na(NPAT_targ[["...11"]])) %>% mutate(factor = "NPAT")

```


```{r}
nonspeckle_ip<-rbind(LMNA_targ_ip_input, PML_targ_ip_input, SP100_targ_ip_input, SMN2_targ_ip_input, WRAP53_targ_ip_input, SAM68_targ_ip_input, FBL_targ_ip_input, NPAT_targ_ip_input)

nonspeckle_ip_pos<-nonspeckle_ip %>% filter(`pulldown vs input log2FoldChange residuals` > 1)

table(nonspeckle_ip_pos[[3]], nonspeckle_ip_pos[[4]])

t_test_speckle <- t.test(
  speckle_gaussian$`pulldown vs input log2FoldChange residuals`,
  speckle_non_gaussian$`pulldown vs input log2FoldChange residuals`
)

```



Ridge code

```{r, warning=FALSE, message=FALSE}
library(ggridges)

speckle_and_non<-rbind(nonspeckle_ip_pos, speckle_ip_pos)

names(speckle_and_non) <- c("gn_name", "log2fc_IP_vs_Input", "gaussian", "target")

speckle_and_non %>%
  ggplot(aes(x = log2fc_IP_vs_Input, y = gaussian)) +
  geom_density_ridges()+
  theme_bw()+
  facet_wrap(~target)


```

Check those validated gaussian introns genes (135 total)

101 are SRSF1 targets

```{r, warning=FALSE}
speckle_and_non<-speckle_and_non %>% mutate(true_g_status = ifelse(speckle_and_non[["gn_name"]] %in% val_RIgau_3reps$gene_name, "truegaussian", "no"))

teste<-speckle_and_non %>% filter(gaussian == "gaussian" & target %in% c("RPS1", "SRSF1", "SRSF7"))

teste_2<-speckle_and_non %>% filter(true_g_status == "truegaussian" & target %in% c("RPS1", "SRSF1", "SRSF7"))

teste_l<-speckle_and_non %>% filter(gaussian == "gaussian" & target %in% c("LMNA"))

teste_2_l<-speckle_and_non %>% filter(true_g_status == "truegaussian" & target %in% c("LMNA"))

table(speckle_and_non[[3]], speckle_and_non[[4]], speckle_and_non[[5]])
```




Ridge
```{r, warning=FALSE, message=FALSE}
speckle_and_non %>%
  ggplot(aes(x = log2fc_IP_vs_Input, y = true_g_status)) +
  geom_density_ridges()+
  theme_bw()+
  facet_wrap(~target)
```

Now data from the Wu and Fei

```{r}
Hepg2_NT_SON<-read_excel("Ref/ref_fei/HEpG2_NT_SON_intron.xlsx")
HELA_NT_SON<-read_excel("Ref/ref_fei/HELA_NT_SON_intron.xlsx")
Heatshock_IRfinder<-read_excel("Ref/ref_fei/Table S6_IRFinder analysis.xlsx")

```

Filter by all gaussian transcripts (1057)

```{r}
Hepg2_NT_SON<-Hepg2_NT_SON %>% 
  dplyr::select(Geneid, log2FoldChange) %>%
  mutate(g_status = ifelse(Hepg2_NT_SON[["Geneid"]] %in% gsinglefinal$gene_id, "gaussian", "no")) %>% mutate(cell = "Hepg2")


HELA_NT_SON<-HELA_NT_SON %>% 
  dplyr::select(Geneid, log2FoldChange) %>% 
  mutate(g_status = ifelse(HELA_NT_SON[["Geneid"]] %in% gsinglefinal$gene_id, "gaussian", "no")) %>% mutate(cell = "HELA")

NT_SON_intron_bind<-rbind(Hepg2_NT_SON, HELA_NT_SON)

NT_SON_intron_bind_up<-NT_SON_intron_bind %>% filter(log2FoldChange > 1)

table(NT_SON_intron_bind_up[[3]], NT_SON_intron_bind_up[[4]])

true_allgets_hepg2<-NT_SON_intron_bind_up %>% filter(cell == "Hepg2" & g_status == "gaussian")

true_allgets_hela<-NT_SON_intron_bind_up %>% filter(cell == "HELA" & g_status == "gaussian")

allgets_hepg2<-NT_SON_intron_bind_up %>% filter(cell == "Hepg2")

allgets_hela<-NT_SON_intron_bind_up %>% filter(cell == "HELA" )

length(intersect(allgets_hela$Geneid, allgets_hepg2$Geneid))

length(intersect(true_allgets_hela$Geneid, true_allgets_hepg2$Geneid))

```
300 from 1057 are enriched on their paraspeckle introns



Now only the gaussian hi confidence
```{r}
NT_SON_intron_bind_up<-NT_SON_intron_bind_up %>% mutate(true_g_status = ifelse(NT_SON_intron_bind_up[["Geneid"]] %in% val_RIgau_3reps$gene_id, "truegaussian", "no"))

table(NT_SON_intron_bind_up[[5]], NT_SON_intron_bind_up[[4]])

true_allgets_hepg2<-NT_SON_intron_bind_up %>% filter(cell == "Hepg2" & true_g_status == "truegaussian")

true_allgets_hela<-NT_SON_intron_bind_up %>% filter(cell == "HELA" & true_g_status == "truegaussian")

allgets_hepg2<-NT_SON_intron_bind_up %>% filter(cell == "Hepg2")

allgets_hela<-NT_SON_intron_bind_up %>% filter(cell == "HELA" )

length(intersect(allgets_hela$Geneid, allgets_hepg2$Geneid))

length(intersect(true_allgets_hela$Geneid, true_allgets_hepg2$Geneid))

```
60 from 135 are nuclear specle introns

venn for this fist the total gaussian
```{r}
library(eulerr)


venn_data <- euler(c(
  GETgenes = 909 -  262,      
  SpeckleEnrichedHELA = 5231- 351 - 262,
  SpeckleEnrichedHepg2 = 4995 - 308 - 262,
  "GETgenes&SpeckleEnrichedHELA" = 351 - 262 ,
  "GETgenes&SpeckleEnrichedHepg2" = 308 - 262,
  "SpeckleEnrichedHELA&SpeckleEnrichedHepg2" = 3131 - 262,
  "GETgenes&SpeckleEnrichedHELA&SpeckleEnrichedHepg2" = 262
))

venn_plot <-plot(venn_data, 
     fills = c("#923346", "grey"),  # Custom colors for the circles
     labels = TRUE,
     quantities = TRUE)
cowplot::ggsave2("Venn_SpeckleEnriched_allGETs.pdf",venn_plot,  height = 5, width = 7)

#now the validated

venn_data <- euler(c(
  GETgenes = 135 -  51,      
  SpeckleEnrichedHELA = 5231- 351 - 51,
  SpeckleEnrichedHepg2 = 4995 - 308 - 51,
  "GETgenes&SpeckleEnrichedHELA" = 351 - 51 ,
  "GETgenes&SpeckleEnrichedHepg2" = 308 - 51,
  "SpeckleEnrichedHELA&SpeckleEnrichedHepg2" = 3131 - 51,
  "GETgenes&SpeckleEnrichedHELA&SpeckleEnrichedHepg2" = 51
))

venn_plot <-plot(venn_data, 
     fills = c("#923346", "grey"),  # Custom colors for the circles
     labels = TRUE,
     quantities = TRUE)
cowplot::ggsave2("Venn_SpeckleEnriched_allGETs.pdf",venn_plot,  height = 5, width = 7)

```

