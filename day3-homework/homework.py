#!/usr/bin/env python3

#############################################################################

# 1. The first step is to load the expression data into a nested list, the tissue names into a list, the gene IDs into a list, and the gene names into a list.

import sys
import numpy

# open File
fs = open(sys.argv[1],mode = 'r')
# skip 2 lines
fs.readline()
fs.readline()
# split colume by tabs and skips two entries
line = fs.readline()
fields = line.strip('\n').split('\t')
# create way to hold gene names
tissues = fields[2:] 
gene_names = []
gene_IDs = []
expression = []
# for each line
for line in fs:
    fields = line.strip('\n').split('\t')     # split line
    gene_IDs.append(fields[0])                 # save field 0 into gene id
    gene_names.append(fields[1])               # save field 1 into gene name 
    expression.append(fields[2:])              # save 2+ into expression values
fs.close() 

#############################################################################

#2. Because nested lists are clunky and difficult to work with, you will be converting the data into numpy arrays.

tissues = numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs)
gene_names = numpy.array(gene_names)
expression = numpy.array(expression,dtype = float)

#############################################################################

# 3. Let’s check the mean expression values for the first 10 genes using an approach you have already learned, nested for loops.

#############################################################################

# 4. Now let’s see how numpy arrays make working with data more streamlined by calculating the same mean expression values for the first 10 genes but using the build-in numpy function mean and printing them.

expression_10_mean = numpy.mean(expression[0:10,:],axis = 1)
print(expression_10_mean)

#############################################################################

# 5. To get a sense of the spread of data, calculate and compare the median and mean expression value of the entire dataset.

expression_median = numpy.median(expression)
expression_mean =  numpy.mean(expression)
print(expression_median)
print(expression_mean)
# observation: the median is 0.0271075, and the mean is 16.557814350910945. 
# The median and mean has a large gap.
# Most sample have very low expression level, but there are some samples having very high expression level that makes mean large.

#############################################################################

# 6. In order to work with a more normalized range of expression values, let’s apply a log-transformation to the data and check the mean and median values again.

pseudo_count = 1
expression_for_log = expression + pseudo_count

expression_log = numpy.log2(expression_for_log)

expression_log_mean = numpy.mean(expression_log)
expression_log_median = numpy.median(expression_log)

print(expression_log_mean)
print(expression_log_median)
# observation: the median is 0.03858718613570538, and the mean is 1.1150342022364093. 
# The median seems not change that much. But the mean is much smaller than the original value.
# The median and mean has a much smaller gap than the original. Make the data more easy to analysize and visualize.

#############################################################################

# 7. Now let’s find the expression gap for each gene between their highest and next highest expression level to identify highly tissue specific genes.

expression_copy = numpy.copy(expression_log)
expression_sort = numpy.sort(expression_copy,axis = 1)
expression_diff = expression_sort[:,-1] - expression_sort[:,-2]  

count = 0
for i in expression_diff:
    #print ("Q7: for gene " + gene_names[count] + ": " + str(expression_diff[count]) )
    count += 1


#############################################################################

# 8. Finally, using the expression difference array you just created, you can now identify genes that show high single-tissue specificity as defined as a difference of at least 10 (~1000-fold difference since the data are log2-transformed).

num_greater_10 = numpy.sum(expression_diff >= 10)

print("Q8: number of genes whose difference is greater than 10 is " + str(num_greater_10))

# Q8: number of genes whose difference is greater than 10 is 33

#############################################################################
