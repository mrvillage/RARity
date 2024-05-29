# RARity

An method for estimating the contribution of rare coding variants to complex trait heritability.

[Paper](https://www.nature.com/articles/s41467-024-45407-8)

## Installation

```r
devtools::install_github("mrvillage/RARity")
```

## Usage

The main function is `rarity::rarity`.
- The first argument is the directory where the data is stored.
- The second argument is a vector of the names of the phenotype files.
- The third argument is a list of vectors of the names of the genes for each chromosome.
 - This list must have 22 elements, one for each chromosome.
 - It is expected that each gene name have a corresponding file in the data directory with the name `gene_name.rkyv.gz`.
- The function then returns a data frame with columns `trait_name`, `chr`, `gene`, `nb_individuals`, `nb_rvs`, `r2`, `adj_r2`, `adj_r2_per_var`, `block_var_r2`, `block_var_adj_r2`.

```r
rarity::rarity("data_dir", c("phenos1.rkyv.gz", "phenos2.rkyv.gz"), list(c("chr1_gene1", "chr1_gene2"), c("chr2_gene1", "chr2_gene2"), ...))
```
