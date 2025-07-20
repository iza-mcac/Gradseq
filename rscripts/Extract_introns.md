---
title: "Extract_exon_and_intron_for_one_transcript"
author: "Iza"
date: "2023-09-19"
output: html_document
---

```{r}
library(GenomicFeatures)
library(GenomicRanges)
library(readr)
library(tidyverse)
```
 
 test using gff portion of make_gff
 
```{r}
gencode_v38_annotation <- read_delim("ref/gencode.v38.chr_patch_hapl_scaff.annotation.gff3.gz", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_names = FALSE, trim_ws = TRUE, skip = 7)
gencode_v38_annotation_tomod<-gencode_v38_annotation
gencode_v38_annotation_tomod[1:12] <-
  gencode_v38_annotation_tomod$X9 %>% str_split_fixed(';', 12)
gff_v38<-cbind(gencode_v38_annotation, gencode_v38_annotation_tomod)
names(gff_v38) <- c("chr", "HAVANA", "type", "left_pos", "right_pos", "exclude_1",
                    "strand", "exclude_2", "etc", "ID", "V1", "V2", "V3", "V4",
                    "V5", "V6", "V7", "V8", "V9", "V10", "V11")
gff_v38_filt <-gff_v38 %>% dplyr::select(1,3:5, 7, 10)
gff_v38_filt$ID <- str_sub(gff_v38_filt$ID, start = 4)
gff_v38_filt$ID <- str_replace(gff_v38_filt$ID, "exon:", "")
chr_info<-gff_v38_filt %>% dplyr::select(1)
df_split <- separate(gff_v38_filt, ID, into = c("tx_id", "exon_count"), sep = ":")
df_split_2<-cbind(df_split, chr_info)
df_split_2<-df_split_2 %>% dplyr::select(-1)
df_split_notgene<-df_split_2 %>% dplyr::filter(type !="gene")

```

merge that with our dic

```{r}
#isoformic::make_tx_to_gene()
txtotx_v38<-txtogene_v38 %>% dplyr::select("transcript_id", "transcript_name")
df_split_notgenedic<-df_split_notgene %>% left_join(txtotx_v38, by = c("tx_id"="transcript_id"))
```

try to do things

```{r}
df_split_notgenedic_plus<-df_split_notgenedic %>% filter(strand == "+")

exon_plus<-df_split_notgenedic_plus %>% filter(type == "exon")

df_split_notgenedic_minus<-df_split_notgenedic %>% filter(strand == "-")

exon_minus<-df_split_notgenedic_minus %>% filter(type == "exon")

```

get all introns for the plus strand

```{r}

calculate_introns_plus <- function(transcript_data) {
  exon_ends <- transcript_data$right_pos
  exon_starts <- c(transcript_data$left_pos[-1], NA)
  intron_start <- exon_ends + 1
  intron_end <- exon_starts - 1
  intron_start[is.na(intron_start)] <- 0  # Replace NA with 0 for the last intron
  intron_data <- data.frame(
    transcript = transcript_data$transcript[1],
    intron_start = intron_start,
    intron_end = intron_end
  )
  intron_data$intron_number <- seq_len(nrow(intron_data))
  return(intron_data)
}

intron_data <- do.call(rbind,
                       lapply(split(exon_plus, exon_plus$tx_id),
                              calculate_introns_plus))

intron_data_plus_strand <- na.omit(intron_data)

plus_tojoin<-exon_plus %>% dplyr::select(chr, transcript_name)

intron_data_plus_strand<-intron_data_plus_strand %>% left_join(distinct(plus_tojoin), c("transcript"="transcript_name"))

```

Now for the minus strand

```{r}

calculate_introns_minus <- function(transcript_data) {
  exon_starts <- transcript_data$right_pos
  exon_ends <- c(NA, transcript_data$left_pos[-length(transcript_data$left_pos)])
  intron_start <- exon_ends -1
  intron_end <- exon_starts - 1
  intron_start[is.na(intron_start)] <- 0 # Replace NA with 0 for the last intron
  intron_data <- data.frame(
    transcript = transcript_data$transcript[1],
    intron_start = intron_start,
    intron_end = intron_end
  )
  intron_data$intron_number <- seq_len(nrow(intron_data))
  return(intron_data)
}

intron_data <- do.call(rbind,
                       lapply(split(exon_minus, exon_minus$tx_id),
                              calculate_introns_minus))

intron_data_minus_strand <- na.omit(intron_data)

minus_tojoin<-exon_minus %>% dplyr::select(chr, transcript_name)

intron_data_minus_strand<-intron_data_minus_strand %>% left_join(distinct(minus_tojoin), c("transcript"="transcript_name"))

intron_data_minus_strand <- subset(intron_data_minus_strand, intron_start != 0)

```

I am a genius it worked

now a function to get the fasta
 
```{r}
library(GRanges)
library(IRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
# Load the human genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Create a GRanges object to specify the region of interest

start_position <- 8949404
end_position <- 8949487
chromosome <- "chr12"


# Extract the sequence
sequence <- getSeq(Hsapiens, chromosome, start_position, end_position)

# Print or manipulate the sequence
cat(as.character(sequence))

```

Filter all Gaussian Transcripts to get a final table

```{r}
gaussian_trans_table_teste_4_dic<-gaussian_trans_table_teste_4 %>% left_join(txtotx_v38, by =c("transcript_id"="transcript_id"))

all_gaussian_minus<-intron_data_minus_strand %>% filter(transcript %in% gaussian_trans_table_teste_4_dic$transcript_name)

all_gaussian_plus<-intron_data_plus_strand %>% filter(transcript %in% gaussian_trans_table_teste_4_dic$transcript_name)


```

make a code that does that

```{r}
my_transcript <-exon_plus %>% filter(transcript_name == "ORC6-205")
transcript_name <- "ORC6-205"
#my_transcript <-exon_minus %>% filter(transcript_name == "MRPS18B-209")
my_exon<- my_transcript %>% filter(exon_count == 1)
start_position <- my_exon$left_pos
end_position <- my_exon$right_pos
chromosome <- my_exon$chr
sequence <- getSeq(Hsapiens, chromosome, start_position, end_position)
sequences <- as.data.frame(as.character(sequence))
sequences<-sequences %>% mutate(transcript_name = transcript_name,
                                exon_number = 1)

names(sequences) <- c("sequence", "transcript_name", "ret_intron_n")

fasta_seqs <- sequences

```

Two subtrack two sequences when it is a partial intron inclusion

```{r}
my_transcript <- exon_minus %>% dplyr::filter(transcript_name == "BSCL2-233")
transcript_name <- "BSCL2-233"
#my_transcript <-exon_minus %>% filter(transcript_name == "MRPS18B-209")
my_exon<- my_transcript %>% filter(exon_count ==7)
start_position <- my_exon$left_pos
end_position <- my_exon$right_pos
chromosome <- my_exon$chr
sequence <- getSeq(Hsapiens, chromosome, start_position, end_position)
sequence
#
my_transcript <- exon_minus %>% filter(transcript_name =="BSCL2-229")
my_exon<- my_transcript %>% filter(exon_count == 8)
start_position <- my_exon$left_pos
end_position <- my_exon$right_pos
chromosome <- my_exon$chr
sequence_sub <- getSeq(Hsapiens, chromosome, start_position, end_position)
sequence_sub

my_transcript <-exon_minus %>% filter(transcript_name == "BSCL2-229")
my_exon<- my_transcript %>% filter(exon_count == 9)
start_position <- my_exon$left_pos
end_position <- my_exon$right_pos
chromosome <- my_exon$chr
sequence_sub_2 <- getSeq(Hsapiens, chromosome, start_position, end_position)
sequence_sub_2


# Subtract sequence2 from sequence1
subtracted_sequence <- gsub(sequence_sub, "", sequence)
#subtracted_sequence <- gsub(sequence_sub_2, "", subtracted_sequence)
#subtracted_sequence <- gsub(sequence_sub_3, "", subtracted_sequence)

#nchar(subtracted_sequence)

rm(fasta_seqs_subtrack)
rm(fasta_seqs)
sequences <- as.data.frame(as.character(subtracted_sequence))
#sequences <- as.data.frame(as.character(sequence))
sequences_SUB<-sequences %>% mutate(transcript_name = transcript_name, exon_number = 7)
names(sequences_SUB) <- c("sequence", "transcript_name", "ret_intron_n")
fasta_seqs_subtrack <- sequences_SUB


```

Percentage of sequence before versus after the intron retention

Best mode: get the main isoform and sum its exons before and after the retention on the RI isoform

```{r}
my_transcript <-exon_minus %>% filter(transcript_name == "BSCL2-229")
max_exon_count <- max(as.numeric(my_transcript$exon_count))
exon_sequences <- list()

for (exon_number in 1:max_exon_count) {
  my_exon <- my_transcript %>% filter(exon_count == exon_number)
  if (nrow(my_exon) > 0) {
    start_position <- my_exon$left_pos[1] 
    end_position <- my_exon$right_pos[1]    
    chromosome <- my_exon$chr[1]  
    exon_sequence <- getSeq(Hsapiens, chromosome, start_position, end_position)
    exon_name <- paste("sequence_exon", exon_number, sep = "")
    exon_sequences[exon_name] <- exon_sequence
  }
}

exon_sequences["sequence_exon1"]


exon_not_present<-paste(#as.character(exon_sequences$sequence_exon1),
                      #as.character(exon_sequences$sequence_exon2)#,
as.character(exon_sequences$sequence_exon3),
#as.character(exon_sequences$sequence_exon4)#,
#as.character(exon_sequences$sequence_exon5)#,
#as.character(exon_sequences$sequence_exon6),
#as.character(exon_sequences$sequence_exon7),
as.character(exon_sequences$sequence_exon8),
as.character(exon_sequences$sequence_exon9),
as.character(exon_sequences$sequence_exon10),
as.character(exon_sequences$sequence_exon11)#,
#as.character(exon_sequences$sequence_exon12)#,
#as.character(exon_sequences$sequence_exon13),
#as.character(exon_sequences$sequence_exon14),
#as.character(exon_sequences$sequence_exon15),
#as.character(exon_sequences$sequence_exon16),
#as.character(exon_sequences$sequence_exon17),
#as.character(exon_sequences$sequence_exon18),
#as.character(exon_sequences$sequence_exon19),
#as.character(exon_sequences$sequence_exon20),
#as.character(exon_sequences$sequence_exon21),
#as.character(exon_sequences$sequence_exon22),
#as.character(exon_sequences$sequence_exon23),
#as.character(exon_sequences$sequence_exon24),
#as.character(exon_sequences$sequence_exon25),
#as.character(exon_sequences$sequence_exon26),
#as.character(exon_sequences$sequence_exon27)#,
#as.character(exon_sequences$sequence_exon29)
)

exon_present <- paste(#0
  as.character(exon_sequences$sequence_exon1),
  as.character(exon_sequences$sequence_exon2),
  #as.character(exon_sequences$sequence_exon3),
  as.character(exon_sequences$sequence_exon4),
  as.character(exon_sequences$sequence_exon5),
  as.character(exon_sequences$sequence_exon6),
  as.character(exon_sequences$sequence_exon7)#,
  #as.character(exon_sequences$sequence_exon8),
  #as.character(exon_sequences$sequence_exon9),
  #as.character(exon_sequences$sequence_exon10),
  #as.character(exon_sequences$sequence_exon11)#,
  #as.character(exon_sequences$sequence_exon12)
                    )

nchar(exon_present)
nchar(exon_not_present)

fasta_seqs_subtrack <-fasta_seqs_subtrack %>% mutate(n_nucleotides_present = nchar(exon_present),
                               n_nucleotides_notpresent = nchar(exon_not_present))

#fasta_seqs_subtrack<-fasta_seqs_subtrack %>% mutate(n_nucleotides_present = "lncRNA",              n_nucleotides_notpresent = "lncRNA")

write.csv(fasta_seqs_subtrack, "manual_hell/BSCL2-233_table.csv")


```







