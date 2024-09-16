#!/usr/bin/env python3

# 2.1

import numpy as np

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

# Define k value for k-mers
k = 3

# Create an empty list to store edges
edges = []

# Loop through each read
for read in reads:
    # Loop through the read to get all k-mers
    for i in range(len(read) - k + 1):
        # Extract the k-mer
        kmer = read[i:i+k]
        # Define source (first k-1 nucleotides) and destination (last k-1 nucleotides)
        source = kmer[:-1]  # First two characters
        destination = kmer[1:]  # Last two characters
        
        print(f"{source} -> {destination}")
        edges.append(f"{source} -> {destination}") # Append the edge in the format source -> destination

# Write edges to a file, one per line
with open("de_bruijn_edges.txt", "w") as f:
    for edge in edges:
        f.write(edge + "\n")


# 2.2
import scipy

# 2.3 

# Provided list of reads
reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 
         'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

# Define k value for k-mers
k = 3

# Create a set to store unique edges
edges = set()

# Loop through each read
for read in reads:
    # Loop through the read to get all k-mers
    for i in range(len(read) - k + 1):
        # Extract the k-mer
        kmer = read[i:i+k]
        # Define source (first k-1 nucleotides) and destination (last k-1 nucleotides)
        source = kmer[:-1]  # First two characters
        destination = kmer[1:]  # Last two characters
        # Add the edge to the set (set ensures uniqueness)
        edges.add((source, destination))

# Write DOT format to a file
with open("de_bruijn_graph.dot", "w") as f:
    f.write("digraph G {\n")  # Start of the DOT file
    for source, destination in sorted(edges):
        f.write(f'    "{source}" -> "{destination}";\n')
    f.write("}\n")  # End of the DOT file

print("De Bruijn graph in DOT format saved to de_bruijn_graph.dot")