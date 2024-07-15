# RARity

An method for estimating the contribution of rare coding variants to complex trait heritability.

[Paper](https://www.nature.com/articles/s41467-024-45407-8)

## Installation

Requires the [Rust programming language](https://rust-lang.org).

```sh 
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then install the package using the following command:

```r
devtools::install_github("mrvillage/RARity")
# OR
install.packages("https://github.com/mrvillage/rarity/archive/refs/heads/master.tar.gz", repos=NULL) # use .zip for Windows
```

## Usage

The main function is `rarity::rarity`.
- The first argument is the directory where the data is stored.
- The second argument is a vector of the names of the phenotype files.
- The function then returns a data frame with columns `pheno_file`, `trait_name`, `chr`, `gene`, `nb_individuals`, `nb_rvs`, `r2`, `adj_r2`, `adj_r2_per_var`, `block_var_r2`, `block_var_adj_r2`.
- The function expects directories with the format `chr_01` to `chr_22` corresponding to each chromosome. If a chromosome is missing, the function will ignore it. Each file within the folder is read as a gene block.
- Using `.RData` files is supported but STRONGLY discouraged as they cannot be loaded in parallel or in chunks so will all be loaded in sequence at the very beginning of execution. This will likely result in RAM issues for large datasets and will be significantly slower than other file types. The recommended file type is `.rkyv.gz` which is a compressed binary format that can be read in parallel and in chunks. For converting, see the functions available in [`lmutils.r`](https://github.com/mrvillage/lmutils.r).

```r
rarity::rarity("data_dir", c("phenos1.rkyv.gz", "phenos2.rkyv.gz"))
```

## Configuration

RARity exposes three global config options that can be set using environment variables or the `rarity` package functions:

- `RARITY_LOG`/`rarity::set_log_level` to set the log level (default: `info`).
- `RARITY_BLOCKS_PER_CHUNK`/`rarity::set_blocks_per_chunk` to set the number of blocks to process in parallel (default: `16`).
- `RAYON_NUM_THREADS`/`rarity::set_num_threads` to set the number of threads to use (default: `num_cpus::get()`).
- `RARITY_MIN_SUM`/`rarity::set_min_sum` to set the minimum sum of gene block columns to consider a gene block (default: `2.0`).
