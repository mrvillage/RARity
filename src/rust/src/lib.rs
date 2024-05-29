use extendr_api::prelude::*;
use lmutils::{get_r2s, IntoMatrix, Transform};
use rayon::prelude::*;

struct Results {
    trait_name: String,
    chr: usize,
    gene: String,
    nb_individuals: usize,
    nb_rvs: usize,
    r2: f64,
    adj_r2: f64,
    adj_r2_per_var: f64,
    block_var_r2: f64,
    block_var_adj_r2: f64,
}

impl Results {
    fn vec_to_df(results: Vec<Results>) -> Robj {
        data_frame!(
            trait_name = results
                .iter()
                .map(|x| x.trait_name.clone())
                .collect::<Vec<_>>(),
            chr = results.iter().map(|x| x.chr).collect::<Vec<_>>(),
            gene = results.iter().map(|x| x.gene.clone()).collect::<Vec<_>>(),
            nb_individuals = results.iter().map(|x| x.nb_individuals).collect::<Vec<_>>(),
            nb_rvs = results.iter().map(|x| x.nb_rvs).collect::<Vec<_>>(),
            r2 = results.iter().map(|x| x.r2).collect::<Vec<_>>(),
            adj_r2 = results.iter().map(|x| x.adj_r2).collect::<Vec<_>>(),
            adj_r2_per_var = results.iter().map(|x| x.adj_r2_per_var).collect::<Vec<_>>(),
            block_var_r2 = results.iter().map(|x| x.block_var_r2).collect::<Vec<_>>(),
            block_var_adj_r2 = results
                .iter()
                .map(|x| x.block_var_adj_r2)
                .collect::<Vec<_>>()
        )
    }
}

/// Run a RARity analysis.
/// `dir` is the directory where the data is stored.
/// `pheno_file` is a character vector of normalized phenotype file names, relative to `dir`.
/// `genes` is a list of character vectors, with the vector at index `i` containing the gene names
/// for chromosome `i`.
/// Returns a data frame with the results.
#[extendr]
pub fn rarity(dir: &str, phenos: &[Rstr], genes: List) -> Result<Robj> {
    let dir = std::path::Path::new(dir);
    let genes = genes
        .iter()
        .map(|(_, x)| x.as_string_vector())
        .collect::<Option<Vec<_>>>();
    if genes.is_none() {
        return Err(Error::from("genes must be a list character vectors"));
    }
    let genes = genes.unwrap();
    if genes.len() != 22 {
        return Err(Error::from("genes must have 22 elements"));
    }

    let phenos = phenos
        .iter()
        .map(|x| {
            lmutils::File::from_path(x.as_str())
                .unwrap()
                .into_matrix()
                .make_parallel_safe()
                .unwrap()
        })
        .collect::<Vec<_>>();
    let mut phenos = phenos
        .into_par_iter()
        .map(|x| x.to_owned().unwrap())
        .collect::<Vec<_>>();
    phenos
        .iter_mut()
        .for_each(|x| x.remove_column_by_name("eid"));
    let pheno_norm = phenos.iter().map(|x| x.as_mat_ref()).collect::<Vec<_>>();
    let results = genes
        .par_iter()
        .enumerate()
        .flat_map(|(chr, g)| {
            let chr = chr + 1;
            g.into_par_iter()
                .filter_map(|gene| {
                    let path = dir.join(format!("{}.rkyv.gz", gene));
                    if std::fs::metadata(&path).is_ok() {
                        let mut block: lmutils::OwnedMatrix<f64> = lmutils::File::from_path(path)
                            .unwrap()
                            .read_matrix(false)
                            .unwrap();
                        block.remove_column_by_name("eid");
                        let block = block.into_matrix();
                        let block = block
                            .nan_to_mean()
                            .min_sum(2.0)
                            .standardization()
                            .transform()
                            .unwrap();
                        let block = block.as_mat_ref().unwrap();
                        if block.ncols() == 0 {
                            return None;
                        }
                        Some(
                            pheno_norm
                                .par_iter()
                                .zip(&phenos)
                                .flat_map(|(pheno_norm, pheno)| {
                                    let r2s = get_r2s(block, *pheno_norm);
                                    let nb_individuals = block.nrows();
                                    let nb_rvs = block.ncols();
                                    let traits = pheno.colnames().unwrap();
                                    r2s.into_par_iter().enumerate().map(move |(i, x)| {
                                        // block_lcl_r2 and block_ucl_r2 are much more complicated to calculate so will not be returned and can be calculated easily in R with the ci.R2 function in MBESS
                                        let r2 = x.r2();
                                        let adj_r2 = x.adj_r2();
                                        let adj_r2_per_var = adj_r2 / nb_rvs as f64;
                                        let n = nb_individuals as f64;
                                        let m = nb_rvs as f64;
                                        let block_var_r2 = ((4.0 * r2)
                                            * (1.0 - r2).powi(2)
                                            * (n - m - 1.0).powi(2))
                                            / ((n.powi(2) - 1.0) * (n + 3.0));
                                        let block_var_adj_r2 =
                                            ((n - 1.0) / (n - m - 1.0)).powi(2) * block_var_r2;
                                        Results {
                                            trait_name: traits[i].to_string(),
                                            chr,
                                            gene: gene.to_string(),
                                            nb_individuals,
                                            nb_rvs,
                                            r2,
                                            adj_r2,
                                            adj_r2_per_var,
                                            block_var_r2,
                                            block_var_adj_r2,
                                        }
                                    })
                                })
                                .collect::<Vec<_>>(),
                        )
                    } else {
                        None
                    }
                })
                .flatten()
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Ok(Results::vec_to_df(results))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rarity;
    fn rarity;
}
