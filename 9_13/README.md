
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
The normal distribution fits well around the mean but diverges in the tails, especially for low-coverage regions, where the Poisson distribution offers a better fit.

# Step 1.5
numbers of zero/number of all the data (see my code in R)
The amount of genome with 0x coverage from your simulation should be close to the Poisson expectation of 0.24%

For λ = 10, the Poisson probability of having 0x coverage is ~0.00454%. Compare this value with the simulated percentage of 0x coverage calculated earlier. The simulated result should be close to this theoretical expectation.
Visually inspect the histogram and the overlaid lines. The normal curve should approximate the center but might diverge at the extremes. The Poisson distribution will provide a better fit in the low-coverage region.

# Step 1.6
numbers of zero/number of all the data (see my code in R)
The amount of genome with 0x coverage from your simulation should be close to the Poisson expectation of 9e-04%

For λ = 30, the Poisson probability of having 0x coverage is 9.36e-14. For a Poisson distribution with 
λ=30, very few positions should have 0x coverage and most positions to have coverage close to the mean of 30. In practice, the simulation should closely match this distribution.

At 30x coverage, the normal distribution becomes a good approximation because the Poisson distribution begins to resemble a normal distribution as the mean increases. The normal distribution should fit the data quite well, though it may not capture the discrete nature of the coverage depths as precisely as the Poisson distribution.

# step 2.4

(base) cmdb@QuantBio-14 ~ % conda activate graphviz
(graphviz) cmdb@QuantBio-14 ~ % cd ~/qbb2024-answers/9_13 
(graphviz) cmdb@QuantBio-14 9_13 % dot -Tpng de_bruijn_graph.dot -o ex2_digraph.png

# step 2.5
one possible answer: ATTGATTGATTCATTGATCATTTCATTGATT

# step 2.6
- Ensure complete coverage of k-mers.
- Use an Eulerian path or circuit to traverse the graph.
- Correct errors and ambiguities for accurate reconstruction.