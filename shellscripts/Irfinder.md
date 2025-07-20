# Annotated IRFinder Workflow Script (English)

## 1. **Set up Singularity environment**
Build the IRFinder container image using Singularity:
```bash
singularity build -F ./IRFinder docker://cloxd/irfinder:2.0
```

## 2. **Download Reference Files**
Create the reference directories:
```bash
mkdir REF
mkdir REF/gencode38
```
Download GENCODE release 38 transcriptome and annotation files:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
```
Unzip both files and rename:
```bash
gunzip gencode.v38.transcripts.fa.gz
gunzip gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
mv gencode.v38.transcripts.fa genome.fa
mv gencode.v38.chr_patch_hapl_scaff.annotation.gtf transcripts.gtf
```

## 3. **Build IRFinder Reference (Old method)**
```bash
singularity exec -B /home/izamamede/isoswitch/IRFinder/ IRFinder IRFinder -m BuildRefProcess -r REF/gencode38
```

## 4. **Build STAR index manually**
```bash
STAR --runThreadN 22 --runMode genomeGenerate \
--genomeDir /path/to/STAR_index \
--genomeFastaFiles genome.fa \
--sjdbGTFfile transcripts.gtf \
--genomeChrBinNbits 8
```

## 5. **Modify GTF for transcript types**
Change lncRNA and retained_intron to processed_transcript:
```bash
awk -F'\t' 'BEGIN{OFS=FS} {if ($3 == "exon") {gsub(/transcript_type \"lncRNA\"/, "transcript_type \"processed_transcript\""); gsub(/transcript_type \"retained_intron\"/, "transcript_type \"processed_transcript\"")} print}' transcripts.gtf > modified_file.gtf
```

## 6. **STAR genomeParameters.txt example**
Adjust STAR parameters for genome index generation.

## 7. **Build IRFinder reference from STAR index**
```bash
singularity run IRFinder -m BuildRefFromSTARRef -r REF/gencode38_new -x REF/STAR_index -M REF/gencode38/Mapability/MapabilityExclusion.bed.gz
```

## 8. **Run IRFinder (fastq input)**
```bash
IRFinder -r REF/gencode38 -t 12 -d res/sample_irfinder sample_1.fq sample_2.fq
```
Using Singularity:
```bash
singularity run IRFinder FastQ -r REF/fullgencode38_new -d align/sample_irfinder -t 20 sample_1.fq sample_2.fq
```

## 9. **Run IRFinder in BAM mode**
```bash
singularity run IRFinder BAM -r REF/fullgencode38_new -d align/sample_irfinder align/sample_irfinder/Sorted.bam
```

## 10. **Gunzip all human fastq files**
```bash
for ACCESSION in `cat index_human.txt`; do
  echo "Processing $ACCESSION";
  gunzip -c rawdata/${ACCESSION}_2.fq.gz > gz_rawdata/${ACCESSION}_2.fq;
done
```

## 11. **Loop to run IRFinder on all samples (FASTQ mode)**
```bash
for ACCESSION in `cat index_human_3.txt`; do
  echo "Processing $ACCESSION";
  singularity run IRFinder FastQ -r REF/small_gencode38 -d align/${ACCESSION}_irfinder -t 20 gz_rawdata/${ACCESSION}_1.fq gz_rawdata/${ACCESSION}_2.fq;
done
```

## 12. **Align everything with STAR**
```bash
for ACCESSION in `cat index_human_3.txt`; do
  echo "Processing $ACCESSION";
  mkdir star_res/${ACCESSION}
  STAR --runThreadN 18 --genomeDir STAR_index/ \
  --readFilesIn gz_rawdata/${ACCESSION}_1.fq gz_rawdata/${ACCESSION}_2.fq \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix star_res/${ACCESSION}/${ACCESSION}_;
done
```

## 13. **Run IRFinder on STAR-aligned BAMs**
```bash
for ACCESSION in `cat index_irf.txt`; do
  singularity run IRFinder -m BAM -r REF/fullgencode38_new -d align/full_${ACCESSION}_irfinder all_bam/${ACCESSION}_Aligned.out.bam;
done
```

## 14. **Run IRFinder on stability dataset**
```bash
for ACCESSION in `cat index.txt`; do
  singularity run IRFinder -m BAM -r REF/fullgencode38_new -d align_stability/full_${ACCESSION}_irfinder res_stability/${ACCESSION}_Aligned.out.bam;
done
```

## 15. **Fix for IRFinder source code**
Ensure `#include <cmath>` exists in `includefine.h` to avoid compile errors.

## 16. **Example SLURM job script (for STAR index build)**
```bash
#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --partition=hpc
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=isaabelaa@gmail.com
#SBATCH --cpus-per-task=22
#SBATCH --output=logs%j.log

STAR --runThreadN 22 --runMode genomeGenerate \
--genomeDir /scratch/STAR_index/ \
--genomeFastaFiles GRCh38.p13.genome.fa \
--sjdbGTFfile gencode.v38.annotation.gtf
```

## 17. **Decompress FASTQ before STAR alignment (if needed)**
```bash
for ACCESSION in `cat index_extract`; do
  gunzip -c fastq/${ACCESSION}_2.fq.gz > fastq/gunzip/${ACCESSION}_2.fq;
done
```

---
