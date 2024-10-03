#!/usr/bin/env python3

# create empty lists
allele_frequencies = []
read_depths = []

for line in open("biallelic.vcf"):
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')

    # find the information
    info = fields[7]
    format = fields[8]
    samples = fields[9:]
    
    info_fields = dict(item.split('=') for item in info.split(';') if '=' in item)


    # Extract allele frequency (AF) from the INFO field, if present
    if 'AF' in info_fields:
        af_values = info_fields['AF'].split(',')
        allele_frequencies.append(float(af_values[0]))

    format_keys = format.split(':')

    for sample in samples:
        sample_data = sample.split(':')
        sample_dict = dict(zip(format_keys, sample_data))
        if 'DP' in sample_dict:
            read_depths.append(int(sample_dict['DP']))

# save the AF data
with open('AF.txt', 'w') as af_file:
    af_file.write('Allele_Frequency\n') 
    for af in allele_frequencies:
        af_file.write(f"{af}\n")

# save the DP data
with open('DP.txt', 'w') as dp_file:
    dp_file.write('Read_Depth\n')
    for dp in read_depths:
        dp_file.write(f"{dp}\n")