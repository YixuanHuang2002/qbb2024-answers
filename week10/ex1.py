#!/usr/bin/env python
import imageio.v3 as iio
import numpy as np
import pandas as pd
import os
from skimage.filters import threshold_otsu

# Define directories and settings
image_directory = '/Users/cmdb/qbb2024-answers/week10/'

gene_names = ['APEX1', 'PIM2', 'SRSF1', 'POLR2B']
fields = ["field0", "field1"] 
channels = ["nascentRNA", "PCNA", "DAPI"] 

# Load images
combined_images = []
for gene in gene_names:
    for field in fields:
        combined_image = []
        for channel in channels:
            file_path = os.path.join(image_directory, f"{gene}_{field}_{channel}.tif")
            image = iio.imread(file_path).astype(np.uint16)
            image = (image - np.amin(image)) / (np.amax(image) - np.amin(image))  # Normalize
            combined_image.append(image)
        stacked_image = np.stack(combined_image, axis=-1)
        combined_images.append(stacked_image)

# Convert to numpy array
combined_images = np.array(combined_images)

# Binary masks
binary_masks = []
for image in combined_images:
    dapi_channel = image[:, :, 2]  # DAPI channel
    threshold = threshold_otsu(dapi_channel)
    mask = dapi_channel > threshold
    binary_masks.append(mask)

binary_masks = np.array(binary_masks)

# Define the find_labels function
def find_labels(mask):
    l = 0
    labels = np.zeros(mask.shape, np.int32)
    equivalence = [0]
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    for x in range(1, mask.shape[0]):
        if mask[x, 0]:
            if mask[x - 1, 0]:
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        if mask[x, -1]:
            if mask[x - 1, -1]:
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = np.array(equivalence)
    for i in range(1, len(equivalence))[::-1]:
        labels[np.where(labels == i)] = equivalence[i]
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels

# Generate labeled maps
label_maps = [find_labels(mask) for mask in binary_masks]
label_maps = np.array(label_maps)

# Filtering function
def filter_by_size(labels, minsize, maxsize):
    sizes = np.bincount(labels.ravel())
    for i in range(1, sizes.shape[0]):
        if sizes[i] < minsize or sizes[i] > maxsize:
            where = np.where(labels == i)
            labels[where] = 0
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels

# Filter label maps by size
minsize = 100
maxsize = 10000000
for i in range(len(label_maps)):
    label_maps[i] = filter_by_size(label_maps[i], minsize=minsize, maxsize=maxsize)

# Final filtering
final_label_maps = []
for label_map in label_maps:
    sizes = np.bincount(label_map.ravel())
    sizes[0] = 0
    mean_size = np.mean(sizes[sizes > 0])
    std_size = np.std(sizes[sizes > 0])
    minsize = max(mean_size - std_size, 100)
    maxsize = mean_size + std_size
    filtered_map = filter_by_size(label_map, minsize, maxsize)
    final_label_maps.append(filtered_map)

final_label_maps = np.array(final_label_maps)

# Calculate results
gene_list = []
field_list = []
nucleus_ids = []
nascent_rna_means = []
pcna_means = []
log2_ratios = []

for i, (label_map, image) in enumerate(zip(final_label_maps, combined_images)):
    gene_index = i // 2
    field_index = i % 2
    gene_name = gene_names[gene_index]
    field_name = fields[field_index]
    pcna_channel = image[:, :, 1]
    nascent_rna_channel = image[:, :, 0]
    num_labels = np.amax(label_map)
    for nucleus_id in range(1, num_labels + 1):
        pixels = np.where(label_map == nucleus_id)
        mean_pcna = pcna_channel[pixels].mean()
        mean_rna = nascent_rna_channel[pixels].mean()
        log2_ratio = np.log2(mean_rna / mean_pcna) if mean_pcna > 0 else np.nan
        gene_list.append(gene_name)
        field_list.append(field_name)
        nucleus_ids.append(nucleus_id)
        nascent_rna_means.append(mean_rna)
        pcna_means.append(mean_pcna)
        log2_ratios.append(log2_ratio)

results_df = pd.DataFrame({
    "Gene": gene_list,
    "Field": field_list,
    "Nucleus ID": nucleus_ids,
    "Mean Nascent RNA Signal": nascent_rna_means,
    "Mean PCNA Signal": pcna_means,
    "Log2 Ratio": log2_ratios
})

results_df.to_csv("nucleus_signals.csv", sep="\t", index=False)
print(results_df.head())