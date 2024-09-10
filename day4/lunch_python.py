#!/usr/bin/env python3

import sys
import numpy as np

#############################################################################

# 1. The first step is to load the gene-tissue pairs from the gene_tissue.tsv file.

fs = open(sys.argv[1],mode = 'r')  # open File
# create dicts to hold samples for gene_tissue pairs
relevant_samples = {}
# step through files
for lines in fs:
    # split lines into fields
    fields = lines.rstrip('\n').split('\t')
    # create key from gene
    key = fields[0]
    # create value from tissue
    value = fields[2]
    # initiate dict from keys with value
    relevant_samples[key] = value
fs.close()


# Mike's code 
filename = sys.argv[1]  
fs = open(filename, mode='r')  # Open the file in read mode

relevant_samples_1 = {}

for line in fs:
     # split lines into fields
    fields = line.rstrip("\n").split("\t")  
    # Create a key from geneID (fields[0]) and tissue type (fields[2])
    key = (fields[0], fields[2])
    # Initialize an empty list for each gene-tissue pair 
    relevant_samples_1[key] = []  
fs.close()  # Close the file after reading all lines


#############################################################################

# 2. Next you will need to figure out which tissue corresponds to which sampleIDs.

fs = open(sys.argv[2],mode = 'r')  # open File
#skip line
fs.readline()
# create dicts to hold samples for gene_tissue pairs
tissue_sample = {}
# step through files
for lines in fs:
    # split lines into fields
    fields = lines.rstrip('\n').split('\t')
    # create key from gene and tissue
    key = fields[6]
    value = fields[0]
    # initiate dict from keys with list to hold samples
    tissue_sample.setdefault(key,[])
    tissue_sample[key].append(value)
fs.close()

#############################################################################

# Question 3

# open File
fs3 = open(sys.argv[3],mode = 'r')  
# skip 2 line
fs3.readline()
fs3.readline()
# pull out the head line and divide it into a list
head = fs3.readline().rstrip('\n').split('\t')
# what we need is the name from the 3rd to the end
header = head[2:]
fs3.close()

#############################################################################

# Question 4 & 5

# create a new dict
tissue_columes = {}

# we pull every tissue sample we are interested in in the tissue_sample dict
for tissue ,samples in tissue_sample.items():
    # use tissue name as the key, and the value should be a list, then create a empty one
    tissue_columes.setdefault(tissue,[])
    # screening through all the samples in the tissue
    for sample in samples:
        # check if the sample is in the head line, if not, we don't want to conclude it because it's value should be a empty list
        if sample in header:
            # find the sample's position in list header
            position = header.index(sample)
            # add the corresponding positions
            tissue_columes[tissue].append(position)

#For each tissue type, see how many samples have expression data.

for tissue in tissue_columes:
    num = len(tissue_columes[tissue])
    print(tissue + ':' + str(num))

# Whole Blood:755
# Brain - Frontal Cortex (BA9):209
# Adipose - Subcutaneous:663
# Muscle - Skeletal:803
# Artery - Tibial:663
# Artery - Coronary:240
# Heart - Atrial Appendage:429
# Adipose - Visceral (Omentum):541
# Ovary:180
# Uterus:142
# Vagina:156
# Breast - Mammary Tissue:459
# Skin - Not Sun Exposed (Suprapubic):604
# Minor Salivary Gland:162
# Brain - Cortex:255
# Adrenal Gland:258
# Thyroid:653
# Lung:578
# Spleen:241
# Pancreas:328
# Esophagus - Muscularis:515
# Esophagus - Mucosa:555
# Esophagus - Gastroesophageal Junction:375
# Stomach:359
# Colon - Sigmoid:373
# Small Intestine - Terminal Ileum:187
# Colon - Transverse:406
# Prostate:245
# Testis:361
# Skin - Sun Exposed (Lower leg):701
# Nerve - Tibial:619
# Heart - Left Ventricle:432
# Pituitary:283
# Brain - Cerebellum:241
# Cells - Cultured fibroblasts:504
# Artery - Aorta:432
# Cells - EBV-transformed lymphocytes:174
# Brain - Cerebellar Hemisphere:215
# Brain - Caudate (basal ganglia):246
# Brain - Nucleus accumbens (basal ganglia):246
# Brain - Putamen (basal ganglia):205
# Brain - Hypothalamus:202
# Brain - Spinal cord (cervical c-1):159
# Liver:226
# Brain - Hippocampus:197
# Brain - Anterior cingulate cortex (BA24):176
# Brain - Substantia nigra:139
# Kidney - Cortex:85
# Brain - Amygdala:152
# Cervix - Ectocervix:9
# Fallopian Tube:9
# Cervix - Endocervix:10
# Bladder:21
# Kidney - Medulla:4
# Cells - Leukemia cell line (CML):0




# Which tissue types have the largest number of samples? The fewest?

# set the default max_lenth and tissue
max_length = 0
max_tissue = None


for tissue in tissue_columes.keys():
    # check everyline if the correspondin length is larger than the previous largest lenghth
    if len(tissue_columes[tissue]) > max_length: # if yes
        # update the max_length with the new number
        max_length = len(tissue_columes[tissue])
        # update the max_tissue with the new tissue
        max_tissue = tissue
print("Tissue types have the largest number of samples is " + max_tissue + ". The number of sample is " + str(max_length))
# Tissue types have the largest number of samples is Muscle - Skeletal. The number of sample is 803

# set the default min_lenth and tissue
min_length = max_length # chances are that min_length could be anything smaller than max_length, so we set it as max_length
min_tissue = None

for tissue in tissue_columes.keys():
    # check everyline if the correspondin length is smaller than the previous smallest lenghth
    if len(tissue_columes[tissue]) < min_length:
        # update the min_length with the new number
        min_length = len(tissue_columes[tissue])
        # update the min_tissue with the new tissue
        min_tissue = tissue
print("Tissue types have the fewest number of samples is " + min_tissue + ". The number of sample is " + str(min_length))
# Tissue types have the fewest number of samples is Cells - Leukemia cell line (CML). The number of sample is 0

#############################################################################
