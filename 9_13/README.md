
# Step 1.1

Genome Size (G): 1,000,000 bp
Desired Coverage (C): 3x
Read Length (L): 100 bp
Number of Reads= 3×1,000,000 / 100 = 30000 Reads

# Step 1.4

# In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
numbers of zero/number of all the data (see my code in R)
The amount of genome with 0x coverage from your simulation should be close to the Poisson expectation of ~4.69%.
# How well does this match Poisson expectations? How well does the normal distribution fit the data?
For λ = 3, the Poisson probability of having 0x coverage is ~4.98%. My simulation aligns well with Poisson expectations.
The normal distribution fits well around the mean. But the Poisson distribution offers a better fit.

# Step 1.5
numbers of zero/number of all the data (see my code in R)
The amount of genome with 0x coverage : 0.24%

For λ = 10, the Poisson probability of having 0x coverage is ~0.00454%. My simulation aligns well with Poisson expectations.
Normal curve performes very well in 10x, and is better than the Poisson distribution.

# Step 1.6
numbers of zero/number of all the data (see my code in R)
The amount of genome with 0x coverage from your simulation should be close to the Poisson expectation of 9e-04

For λ = 30, the Poisson probability of having 0x coverage is 9.36e-14. For a Poisson distribution with 
λ=30, very few positions should have 0x coverage and most positions to have coverage close to the mean of 30. The simulation match this distribution.

At 30x coverage, the normal distribution and Poisson distribution are very similar. Both of them fit the data quite well. Poisson distribution may be a little bit better in 30x.

# step 2.4

(base) cmdb@QuantBio-14 ~ % conda activate graphviz
(graphviz) cmdb@QuantBio-14 ~ % cd ~/qbb2024-answers/9_13 
(graphviz) cmdb@QuantBio-14 9_13 % dot -Tpng de_bruijn_graph.dot -o ex2_digraph.png

# step 2.5
one possible answer: ATTGATTCTTATTCATTCATTT

# step 2.6
- Ensure complete coverage of k-mers.
- take longer k-mers.
- Collect more data (more reads), so that every 
