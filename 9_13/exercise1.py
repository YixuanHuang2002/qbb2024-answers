#!/usr/bin/env python3

# Step 1.2

import numpy as np

# Parameters
genome_size = 1000000  
read_length = 100  
coverage = 3  
num_reads = (coverage * genome_size) // read_length # have to be an integer

# create a zero coverage array
coverage_array = np.zeros(genome_size, dtype= int)

# generate random 
for i in range(num_reads):
    start = np.random.randint(0, genome_size - read_length )
    coverage_array[start:start + read_length] += 1

# Save coverage array
np.savetxt("coverage_3x.txt", coverage_array, fmt="%d")

# 1.5 repeat the analysis with 10x coverage(just change parameter)

genome_size = 1000000  
read_length = 100  
coverage = 10
num_reads = (coverage * genome_size) // read_length # answer of first questions - num of read


coverage_array = np.zeros(genome_size)


for _ in range(num_reads):
    start_position = np.random.randint(0, genome_size - read_length + 1)
    coverage_array[start_position:start_position + read_length] += 1

# Save the coverage data
np.savetxt("coverage_10x.txt", coverage_array)

# 1.6 repeat the analysis with 30x coverage

genome_size = 1000000  
read_length = 100  
coverage = 30
num_reads = (coverage * genome_size) // read_length # answer of first questions - num of read


coverage_array = np.zeros(genome_size)


for _ in range(num_reads):
    start_position = np.random.randint(0, genome_size - read_length + 1)
    coverage_array[start_position:start_position + read_length] += 1

# Save the coverage data
np.savetxt("coverage_30x.txt", coverage_array)
