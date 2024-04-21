# Article Steps

## Input/Output

**Input:** Set of _genetic variants_, _LD blocks_.

- _Genetic variants_ - We need a list of all common SNPs from the ancestry we're looking at (-snp *.snp)
- 

**Output:** Statistical intersection of the loci with dataset in a compendium.

## Steps

1. Identify variants with linkage (r^2 > 0.8) to any input variant in major ancestry (EU, African, Asian--input data only includes EU) - done L1134 (can also include LD file).
2. Look at LD blocks and dataset and count number of observed intersections.
3. Look at LD blocks and dataset and calculate the expected intersection max_rep times. This produces a null distribution:
    - Find most strongly associated variant for each LD block
    - Find distance of each variant from this strong variant
    - Randomly choose a genomic variant w/ allele frequencies that match the reference variant from dbSNP

    - Genomic coordinates of artificial variants are created that are located at the same relative distances from this random variant using the distance vector. Members of this artificial LD block are intersected with each dataset, as was done for the observed intersections. This strategy accounts for the number of variants in the input LD block and their relative distances, while prohibiting ‘double counting’ due to multiple variants in the block intersecting the same dataset.
4. The expected intersection distributions are used to calculate z-scores and P-values for the observed intersection. The final reported P-values are Bonferroni-corrected (Pc) for the 1,544 TF datasets tested.
5. Relative risk = observed intersection / mean expected intersection
