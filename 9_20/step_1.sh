#!/bin/bash

# Step1: Obtain the data

## step 1.4 sort and merge
for FEATURE in exons genes cCREs #Loops through the three features and sorts them, then merges duplicates
do
    bedtools sort -i ${FEATURE}.bed > ${FEATURE}_sort.bed
    bedtools merge -i  ${FEATURE}_sort.bed > ${FEATURE}_chr1.bed
done
## step 1.5 get introns

bedtools subtract -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed

## step 1.6
bedtools subtract -a genome_chr1.bed -b exons_chr1.bed cCREs_chr1.bed introns_chr1.bed > other_chr1.bed
