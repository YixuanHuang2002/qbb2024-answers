#!/user/bin/env bash

wget https://www.dropbox.com/s/ogx7nbhlhrp3ul6/BYxRM.tar.gz
tar -xvzf BYxRM.tar.gz

### Question 1.1 ###
awk 'NR % 4 == 2 {print length ($0) }' A01_09.fastq > read_length.txt # count leagth of each read, all are 76
#sequencing reads length is 76 bp

### Question 1.2 ###
echo $(cat A01_09.fastq|wc -l)/4|bc # there are 4 lines for each read
# 669548 reads are present within the file

### Question 1.3 ###
expr \( 114 \* 669548 \) / 12070000
# average depth of coverage: 6

### Question 1.4 ###
du -sh *.fastq | sort -n
# largest: A01_62.fastq 149MB
# samllest: A01_27.fastq 110M

### Question 1.5 ###
fastqc A01_09.fastq
# What is the median base quality along the read? ~35
# How does this translate to the probability that a given base is an error? 10^(-3.5) = 0.032%
# the overall quality remains high, with little variation

