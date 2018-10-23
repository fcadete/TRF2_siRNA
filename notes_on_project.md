---
title: "Notes on the ongoing analysis of TRF2 siRNA data"
output: html_document
author: Filipe Tavares-Cadete
---

## Data

### Summary

This data was received from NovoGene by hard drive on the 17th of October, 2018. This is the result of a RNA-seq experiment done in the CMAzzalin lab at Instituto de Medicina Molecular, Lisboa, Portugal. There are two short-hairpin RNA conditions (control and TRF2 shRNA) and three time points after inducing short-hairpin expression with Dox (0h, 48h, 96h), with three replicates per combination of condition and time point, for a total of 18 samples. Each sample was sequenced by NovoGene using Illumina paired-end sequencing, yielding two files per sample, for a total of 36 fq.gz files in the `raw_data` folder.

### People involved

 + Experimentalists: Bruno Adriano de Sousa Silva, Ana Margarida da Silva Figueira
 + Principal investigator: Claus M. Azzalin
 + Bioinformatician: Filipe Tavares-Cadete

### File description

Filename | siRNA | TimePoint | Replicate | Pair
--- | --- | --- | --- | ---
C_d48h_1_1.fq.gz | control | 48h | 1 | 1
C_d48h_1_2.fq.gz | control | 48h | 1 | 2
C_d48h_2_1.fq.gz | control | 48h | 2 | 1
C_d48h_2_2.fq.gz | control | 48h | 2 | 2
C_d48h_3_1.fq.gz | control | 48h | 3 | 1
C_d48h_3_2.fq.gz | control | 48h | 3 | 2
C_d96h_1_1.fq.gz | control | 96h | 1 | 1
C_d96h_1_2.fq.gz | control | 96h | 1 | 2
C_d96h_2_1.fq.gz | control | 96h | 2 | 1
C_d96h_2_2.fq.gz | control | 96h | 2 | 2
C_d96h_3_1.fq.gz | control | 96h | 3 | 1
C_d96h_3_2.fq.gz | control | 96h | 3 | 2
C_NT_1_1.fq.gz | control | 0h | 1 | 1
C_NT_1_2.fq.gz | control | 0h | 1 | 2
C_NT_2_1.fq.gz | control | 0h | 2 | 1
C_NT_2_2.fq.gz | control | 0h | 2 | 2
C_NT_3_1.fq.gz | control | 0h | 3 | 1
C_NT_3_2.fq.gz | control | 0h | 3 | 2
T_d48h_1_1.fq.gz | TRF2 | 48h | 1 | 1
T_d48h_1_2.fq.gz | TRF2 | 48h | 1 | 2
T_d48h_2_1.fq.gz | TRF2 | 48h | 2 | 1
T_d48h_2_2.fq.gz | TRF2 | 48h | 2 | 2
T_d48h_3_1.fq.gz | TRF2 | 48h | 3 | 1
T_d48h_3_2.fq.gz | TRF2 | 48h | 3 | 2
T_d96h_1_1.fq.gz | TRF2 | 96h | 1 | 1
T_d96h_1_2.fq.gz | TRF2 | 96h | 1 | 2
T_d96h_2_1.fq.gz | TRF2 | 96h | 2 | 1
T_d96h_2_2.fq.gz | TRF2 | 96h | 2 | 2
T_d96h_3_1.fq.gz | TRF2 | 96h | 3 | 1
T_d96h_3_2.fq.gz | TRF2 | 96h | 3 | 2
T_NT_1_1.fq.gz | TRF2 | 0h | 1 | 1
T_NT_1_2.fq.gz | TRF2 | 0h | 1 | 2
T_NT_2_1.fq.gz | TRF2 | 0h | 2 | 1
T_NT_2_2.fq.gz | TRF2 | 0h | 2 | 2
T_NT_3_1.fq.gz | TRF2 | 0h | 3 | 1
T_NT_3_2.fq.gz | TRF2 | 0h | 3 | 2

## FastQC reports

### File by file

I ran FastQC for every file in the dataset. Here I will summarise the results, sample by sample (so combining the two pairs per sample).

+ `C_d48h_1`: Minimum per base sequence quality drops bellow 30 from 75bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 43.74% - 45.41% reads remaining if deduplicated.
+ `C_d48h_2`: Minimum per base sequence quality drops bellow 30 from 50bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 50.6% - 53.68% reads remaining if deduplicated.
+ `C_d48h_3`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-end and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 53.47% - 59.57% reads remaining if deduplicated.
+ `C_d96h_1`: Minimum per base sequence quality drops bellow 30 from 75bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 45.5% - 48.23% reads remaining if deduplicated.
+ `C_d96h_2`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-end and 75bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 48.51% - 50.13% reads remaining if deduplicated.
+ `C_d96h_3`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-end and 75bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 51.56% - 53.2% reads remaining if deduplicated.
+ `C_NT_1`: Minimum per base sequence quality drops bellow 30 from 60bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 54.02% - 56.68% reads remaining if deduplicated.
+ `C_NT_2`: Minimum per base sequence quality drops bellow 30 from 50bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 45.65% - 47.79% reads remaining if deduplicated.
+ `C_NT_3`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-end and 75bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 45.12% - 48.1% reads remaining if deduplicated.
+ `T_d48h_1`: Minimum per base sequence quality drops bellow 30 from 60bp on both pair-ends. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 48.24% - 50.91% reads remaining if deduplicated.
+ `T_d48h_2`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-ends and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 46.22% - 47.97% reads remaining if deduplicated.
+ `T_d48h_3`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-ends and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 47.91% - 50.68% reads remaining if deduplicated.
+ `T_d96h_1`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-ends and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 46.42% - 48.51% reads remaining if deduplicated.
+ `T_d96h_2`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-ends and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 50.07% - 55.52% reads remaining if deduplicated.
+ `T_d96h_3`: Minimum per base sequence quality drops bellow 30 from 50bp on first pair-ends and 70bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 47.79% - 49.69% reads remaining if deduplicated.
+ `T_NT_1`: Minimum per base sequence quality drops bellow 30 from 60bp on first pair-ends and 90bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 46.4% - 48.72% reads remaining if deduplicated.
+ `T_NT_2`: Minimum per base sequence quality drops bellow 30 from 50bp on first pair-ends and 75bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 53.68% - 56.28% reads remaining if deduplicated.
+ `T_NT_3`: Minimum per base sequence quality drops bellow 30 from 50bp on first pair-ends and 75bp on second pair-end. Usual RNA-seq random hexamer distortion in per base sequence content in the first bases. Shifted and non-normal per sequence GC content. 49.39% - 51.8% reads remaining if deduplicated.

### General remarks

There is a drop in sequence quality but not worrisome. Very few positions below a Phred score of 20 (1% chance of error). What is a bit concerning is the GC content and duplication levels. There is little adapter contamination, so my guess is that we have a lot of ribossomal or mithocondrial RNA.

### Mitigating steps

+ Run Trimmomatic to remove low quality bases and adapters (if present);
+ Pay attention to rRNA and mtRNA hits after mapping.

## Data after Trimmomatic

### General remarks

Ran Trimmomatic on the raw data. The code is in the file `scripts/trimmomatic_filtering.sbatch`. The following Trimmomatic parameters were used:

```
ILLUMINACLIP:"/NGStools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
ILLUMINACLIP:"/NGStools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa":2:30:10 \
SLIDINGWINDOW:5:20 \
MINLEN:50
```

These will filter out parts of reads that contain adapters or a sliding window of five nucleotides with a Phred score under 20. If after filtering any read has a length below 50bp the entire read is excluded.

### Table with percentages of filtered reads

Sample | Total reads | Both surviving | Forward only surviving | Reverse only surviving | Dropped
--- | --- | --- | --- | --- | ---
C_NT_1 | 56007528 | 49544618 (88.46%) | 3637447 (6.49%) | 1680749 (3.00%) | 1144714 (2.04%)
C_NT_2 | 52972895 | 45968209 (86.78%) | 3985717 (7.52%) | 1646321 (3.11%) | 1372648 (2.59%)
C_NT_3 | 54933944 | 48747199 (88.74%) | 3489874 (6.35%) | 1654793 (3.01%) | 1042078 (1.90%)
C_d48h_1 | 52250236 | 46039091 (88.11%) | 3447020 (6.60%) | 1645480 (3.15%) | 1118645 (2.14%)
C_d48h_2 | 51392522 | 44873544 (87.32%) | 3360284 (6.54%) | 1952072 (3.80%) | 1206622 (2.35%)
C_d48h_3 | 52616241 | 46686730 (88.73%) | 3470310 (6.60%) | 1453862 (2.76%) | 1005339 (1.91%)
C_d96h_1 | 53902982 | 47842941 (88.76%) | 3487130 (6.47%) | 1507478 (2.80%) | 1065433 (1.98%)
C_d96h_2 | 52892121 | 46900394 (88.67%) | 3271019 (6.18%) | 1745919 (3.30%) | 974789 (1.84%)
C_d96h_3 | 52085349 | 46506559 (89.29%) | 2901899 (5.57%) | 1751393 (3.36%) | 925498 (1.78%)
T_NT_1 | 53064408 | 46541131 (87.71%) | 3547404 (6.69%) | 1781607 (3.36%) | 1194266 (2.25%)
T_NT_2 | 51913373 | 45762932 (88.15%) | 3355360 (6.46%) | 1798006 (3.46%) | 997075 (1.92%)
T_NT_3 | 51969224 | 45697535 (87.93%) | 3188372 (6.14%) | 2012843 (3.87%) | 1070474 (2.06%)
T_d48h_1 | 51405408 | 45222885 (87.97%) | 3381514 (6.58%) | 1683474 (3.27%) | 1117535 (2.17%)
T_d48h_2 | 54218977 | 48362456 (89.20%) | 3160842 (5.83%) | 1711996 (3.16%) | 983683 (1.81%)
T_d48h_3 | 51360243 | 45890571 (89.35%) | 2861711 (5.57%) | 1707410 (3.32%) | 900551 (1.75%)
T_d96h_1 | 52468176 | 46784592 (89.17%) | 2861698 (5.45%) | 1838412 (3.50%) | 983474 (1.87%)
T_d96h_2 | 55515089 | 49220034 (88.66%) | 3452668 (6.22%) | 1741113 (3.14%) | 1101274 (1.98%)
T_d96h_3 | 54431968 | 48130521 (88.42%) | 3240831 (5.95%) | 1970999 (3.62%) | 1089617 (2.00%)

### Results

Trimmomatic did not significantly reduce the number of reads we had in each sample. The lowest percentage of kept full pairs was 88.78% for C_NT_2, with all other samples hovering between 87 and 89%.

After filtering I ran FastQC on the filtered data. Per base sequence quality greatly improved, as expected. There were no significant changes in the per sequence GC content and number of unique sequences. I'm still convinced that they come from rRNA and mtRNA sequences.

### Next steps

+ Align the reads to a rRNA database for further filtering;
+ I'm also running TopHat with the Trimmomatic-filtered dataset to see if we get meaningful results there.


## Mapping to rRNA

### Procedure

Download all rRNA sequences from SILVA that belong to _Homo sapiens_. From there I created a Bowtie2 index and mapped the filtered reads to this index. Took care to save the unmapped reads to a file. Once this is done I will be able to quantify the levels of rRNA in our samples and rerun the FastQC to see if the peaks in per read GC content disappear.

### Results

There was very little mapping to the rRNA sequences using Bowtie2 -- less than 1% mapping there for all samples. The peaks in GC-content did not disappear.

### Next steps

+ TopHat is still running, but the first samples to finish show great mapping. Maybe I can filter after mapping to the genome/transcriptome?
+ On specialised forums bioinformaticians say SortMeRNA has better results than Bowtie2 when mapping to rRNA. Trying to run that but running into weird error messages.

## SortMeRNA

### Procedure

Used SortMeRNA with the SILVA and RFAM databases to identify rRNA reads in our samples. Afterwards ran FastQC on the non-rRNA files.

### Results

Each sample had around 1.5% of rRNA. This was not enough to remove the peaks in GC content. I'm thinking it also includes tRNAs.

### Next steps

+ TopHat still running.
+ Started Bowtie2 to run against a fasta file extracted from UCSC that has rRNA and tRNA.



