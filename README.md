# RARity

An method for estimating the contribution of rare coding variants to complex trait heritability.

[Paper](https://www.nature.com/articles/s41467-024-45407-8)

## Installation

Requires the [Rust programming language](https://rust-lang.org).

```r
devtools::install_github("mrvillage/RARity")
# OR
install.packages("https://github.com/mrvillage/rarity/archive/refs/heads/master.tar.gz", repos=NULL) # use .zip for Windows
```

## Usage

The main function is `rarity::rarity`.
- The first argument is the directory where the data is stored.
- The second argument is a vector of the names of the phenotype files.
- The function then returns a data frame with columns `trait_name`, `chr`, `gene`, `nb_individuals`, `nb_rvs`, `r2`, `adj_r2`, `adj_r2_per_var`, `block_var_r2`, `block_var_adj_r2`.
- The function expects directories with the format `chr_01` to `chr_22` corresponding to each chromosome. If a chromosome is missing, the function will ignore it. Each file within the folder is read as a gene block.

```r
rarity::rarity("data_dir", c("phenos1.rkyv.gz", "phenos2.rkyv.gz"))
```

RARity also reads three environment variables:
- `RARITY_LOG` to set the log level (default: `info`).
- `RARITY_BLOCKS_PER_CHUNK` to set the number of blocks to process in parallel (default: `16`).
- `RAYON_NUM_THREADS` to set the number of threads to use (default: `num_cpus::get()`).
