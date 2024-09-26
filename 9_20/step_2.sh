#!/bin/bash

# Step2: Count feature SNPs and determine enrichment

# Create the results file with a header
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

for MAF in 0.1 0.2 0.3 0.4 0.5 # Loop through each possible MAF value
do 
    MAF_file=chr1_snps_${MAF}.bed # file name for the SNP MAF file

    # Find the SNP coverage of the whole chromosome
    bedtools coverage -a genome_chr1.bed -b ${MAF_file} > coverage_${MAF}.txt 

    # Sum SNPs from coverage
    num_SNP=$(awk '{s+=$4}END{print s}' coverage_${MAF}.txt)
    # Sum total bases from coverage
    num_bases=$(awk '{s+=$6}END{print s}' coverage_${MAF}.txt)

    # Calculate the background density
    density=$(bc -l -e ${num_SNP}/${num_bases})

    for feature in exons introns cCREs other # Loop through each feature name
    do
        # Use the feature value to get the file name for the feature file
        feature_file=${feature}_chr1.bed

        # Find the SNP coverage of the current feature
        bedtools coverage -a ${feature_file} -b ${MAF_file} > ${MAF}_${feature}.txt

        # Find the SNP coverage & total bases of the current feature
        feature_num_SNP=$(awk '{s+=$4}END{print s}' ${MAF}_${feature}.txt)
        feature_num_bases=$(awk '{s+=$6}END{print s}' ${MAF}_${feature}.txt)
        echo ${feature_num_bases}
        # Calculate density
        feature_density=$(bc -l -e ${feature_num_SNP}/${feature_num_bases})
        # Calculate enrichment
        enrichment=$(bc -l -e ${feature_density}/${density})

        # Save result to results file
        echo -e "${MAF}\t${feature}\t${enrichment}"  >> snp_counts.txt
        
    done
done