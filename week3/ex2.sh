#!/user/bin/env bash

# step 2.1 -------------------------------------------------------------------

wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz

# create an index for the sacCer3.fa reference.
bwa index sacCer3.fa

### Question 2.1 ###
grep -c "^>" sacCer3.fa
# there are 17 entries. 16 yeast chr and one mitochondira chr

# step 2.2 -------------------------------------------------------------------

for sample in A01_*.fastq
do
    sample=$(basename "${sample}" .fastq) # remove the extension
    echo "${sample}"
    bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}" sacCer3.fa ${sample}.fastq > ${sample}.sam # alignment file
    samtools sort -@ 4 -O bam -o ${sample}.bam ${sample}.sam #  sort the alignment file and convert it into a bam file
    samtools index ${sample}.bam
done

# step 2.3 -------------------------------------------------------------------

### Question 2.2 ###
less -S A01_09.sam # 20 description lines
echo "$(wc -l A01_09.sam | cut -b 1-9) - 20" | bc
# How many total read alignments are recorded in the SAM file? 669548

### Question 2.3 ###
grep -v '^@' A01_09.sam | grep -w "chrIII" | wc -l # count the lines with chrIII & filters out header lines
# How many of the alignments are to loci on chromosome III? 18195

# step 2.4 -------------------------------------------------------------------
# this step is done with step 2.2

# step 2.5 -------------------------------------------------------------------

### Question 2.4 ###
#  Does the depth of coverage appear to match that which you estimated in Step 1.3? Why or why not?
# Match. Generally, there are 3-8 reads cover each nucleotide, with the average 6. 

### Question 2.4 ###
# Set your window to chrI:113113-113343 (paste that string in the search bar and click enter). How many SNPs do you observe in this window? Are there any SNPs about which you are uncertain? Explain your answer.
# 3 SNP. The SNP on location 113,326 is not 100% to be a SNP. There are only two reads covering the loci, making the possibility of misread or heterozygous higher. 

### Question 2.4 ###
# Set your window to chrIV:825548-825931. What is the position of the SNP in this window? Does this SNP fall within a gene?
# chrIV:825,834. No, it's not on a gene, it is between SCC2 and SAS4.