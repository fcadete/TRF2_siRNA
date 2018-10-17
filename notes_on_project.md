---
title: "Notes on the ongoing analysis of TRF2 siRNA data"
output: html_document
author: Filipe Tavares-Cadete
---

# Data

## Summary

This data was received from NovoGene by hard drive on the 17th of October, 2018. This is the result of a RNA-seq experiment done in the CMAzzalin lab at Instituto de Medicina Molecular, Lisboa, Portugal. There are two short-hairpin RNA conditions (control and TRF2 shRNA) and three time points after inducing short-hairpin expression with Dox (0h, 48h, 96h), with three replicates per combination of condition and time point, for a total of 18 samples. Each sample was sequenced by NovoGene using Illumina paired-end sequencing, yielding two files per sample, for a total of 36 fq.gz files in the `raw_data` folder.

## People involved

 + Experimentalists: Bruno Adriano de Sousa Silva, Ana Margarida da Silva Figueira
 + Principal investigator: Claus M. Azzalin
 + Bioinformatician: Filipe Tavares-Cadete

## File description

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

# FastQC reports

## File by file

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

## General remarks

There is a drop in sequence quality but not worrisome. Very few positions below a Phred score of 20 (1% chance of error). What is a bit concerning is the GC content and duplication levels. There is little adapter contamination, so my guess is that we have a lot of ribossomal or mithocondrial RNA.

## Mitigating steps

+ Run Trimmomatic to remove low quality bases and adapters (if present);
+ Pay attention to rRNA and mtRNA hits after mapping.


