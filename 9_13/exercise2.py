#!/usr/bin/env python3

# 2.1

import numpy as np

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

k = 3

# Create an empty set to store edges
edges = set()

for read in reads:
    for i in range(len(read) - k):
        kmer1 = read[i:i+k]
        kmer2 = read[i+1:i+1+k]
        edge = f"{kmer1} -> {kmer2}" # format source -> destination
        edges.add(edge)  # Append the edge in the format source -> destination


# Write edges to a file
with open("de_bruijn_edges.txt", "w") as f:
    for edge in edges:
        f.write(edge + "\n")


# 2.2
import scipy

# 2.3 

k = 3

# Create a set to store edges
edges = set()

for read in reads:
    for i in range(len(read) - k):
        kmer1 = read[i:i+k]
        kmer2 = read[i+1:i+1+k]
        edge = f"{kmer1} -> {kmer2}" # format source -> destination
        edges.add(edge)  # Append the edge in the format source -> destination


# Write DOT format to a file (the outcome is not correct)
#
#    f.write("digraph G {\n")  
#    for edge in edges:
#        kmer1, kmer2 = edge.split(' -> ')
#        f.write(f'    "{kmer1}" -> "{kmer2}";\n')
#    f.write("}\n")

def write_dot_file(edges, output_file):
    with open(output_file, 'w') as f:
        f.write("digraph G {\n")
        for edge in edges:
            kmer1, kmer2 = edge.split(' -> ')
            f.write(f'    "{kmer1}" -> "{kmer2}";\n')
        f.write("}\n")

write_dot_file(edges, 'de_bruijn_graph.dot')
