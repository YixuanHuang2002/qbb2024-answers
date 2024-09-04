# day2-lunch exercise


## Answer 1

- `head hg38-gene-metadata-feature.tsv | uniq`
- `cut -f7 hg38-gene-metadata-feature.tsv |sort | uniq -c`


- `cut -f7 hg38-gene-metadata-feature.tsv |sort | uniq -c| grep "protein_coding"`
- 19618 protein_coding genes

- biotype I want to learn: protein_coding genes, because I am interested in protein strcture, and I only care about the genes can be translated into a protein.


## Answer 2

- `cut -f1 hg38-gene-metadata-go.tsv | sort |uniq -c | sort -n`
- the most go_ids gene_id: ENSG00000168036, number of go_id: 273

- `grep -w "ENSG00000168036"  hg38-gene-metadata-go.tsv > ENSG00000168036_num.txt | sort -k3 > ENSG00000168036_num.txt `
- DNA replicate regulation and Cell duplicate, strctural related protein


## Answer 3

- `grep -w -e "IG...gene" -e "IG....gene" gene.gtf | cut -f 1 | sort | uniq -c`
-  Results: 91 chr14
  16 chr15
   6 chr16
  52 chr2
   1 chr21
  48 chr22

- `grep -w -e "IG.pseudogene" -e "IG...pseudogene" gene.gtf | cut -f 1 | sort | uniq -c`
- Results: 
1 chr1
   1 chr10
  84 chr14
   6 chr15
   8 chr16
   1 chr18
  45 chr2
  48 chr22
   1 chr8
   5 chr9

- compare: chr2 ,chr22, chr15, chr16 have many both IG genes and IG pseudogenes. But chr14, chr10, chr18,chr8, and chr9 that only have IG pseudogenes but no IG genes. chr21 has only IG genes. In summary, there are a lot of pseudogenes in chromosomes where IG genes do not exist.



## Answer 4

- Other colume (e.g. description of the genes) also contain the word "pseudogenes", so directly grep pseudogene may lead to many unrelated results.
- Betther pattern: specify which colume you are looking for `grep 'gene_type"pseudogene"' gene.gtf`

## Answer 5

- `cut -f 1,4,5,14 gene-tabs.gtf > gene-tabs.bed`