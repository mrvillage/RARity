use extendr_api::prelude::*;
use lmutils::{get_r2s, IntoMatrix, Transform};
use log::{debug, info};
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
/// for chromosome `i`.
/// Returns a data frame with the results.
/// @export
#[extendr]
pub fn rarity(dir: &str, phenos: &[Rstr]) -> Result<Robj> {
    let _ =
        env_logger::Builder::from_env(env_logger::Env::default().filter_or("RARITY_LOG", "info"))
            .try_init();

    let dir = std::path::Path::new(dir);

    info!("Reading data from {}", dir.display());

    let genes = (1..=22)
        .map(|chr| {
            let chr_dir = dir.join(format!("chr_{:02}", chr));
            if std::fs::metadata(&chr_dir).is_ok() {
                std::fs::read_dir(chr_dir)
                    .unwrap()
                    .filter_map(|x| x.ok())
                    .map(|x| x.path())
                    .filter(|x| x.is_file())
                    .map(|x| {
                        let file = lmutils::File::from_path(x.as_path())
                            .unwrap()
                            .into_matrix()
                            .make_parallel_safe()
                            .unwrap();
                        (x, file)
                    })
                    .collect::<Vec<_>>()
            } else {
                vec![]
            }
        })
        .collect::<Vec<_>>();

    info!("Reading phenotypes");

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

    info!("Calculating RARity");

    let blocks_per_chunk = std::env::var("RARITY_BLOCKS_PER_CHUNK")
        .unwrap_or("64".to_string())
        .parse::<usize>()
        .unwrap();

    let genes = genes
        .into_par_iter()
        .enumerate()
        .flat_map(|(chr, g)| {
            let chr = chr + 1;
            g.into_par_iter().map(move |g| (chr, g.0, g.1))
        })
        .collect::<Vec<_>>();

    let mut chunks = vec![];
    let mut chunk = vec![];
    for gene in genes {
        chunk.push(gene);
        if chunk.len() == blocks_per_chunk {
            chunks.push(chunk);
            chunk = vec![];
        }
    }

    let mut results = vec![];
    for chunk in chunks {
        results.par_extend(
            chunk
                .into_par_iter()
                .filter_map(|(chr, path, mat)| {
                    let gene = path.file_name().unwrap().to_str().unwrap();
                    info!("Processing gene {}", gene);
                    if std::fs::metadata(&path).is_ok() {
                        let mut block: lmutils::OwnedMatrix<f64> = mat.to_owned().unwrap();
                        debug!("Removing eid column");
                        block.remove_column_by_name("eid");
                        let block = block.into_matrix();
                        debug!("Normalizing block");
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
                        let results = pheno_norm
                            .par_iter()
                            .zip(&phenos)
                            .flat_map(|(pheno_norm, pheno)| {
                                info!("Calculating R2 for gene {}", gene);
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
                                    let block_var_r2 =
                                        ((4.0 * r2) * (1.0 - r2).powi(2) * (n - m - 1.0).powi(2))
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
                            .collect::<Vec<_>>();
                        info!("Processed gene {}", gene);
                        Some(results)
                    } else {
                        info!("Gene {} not found", gene);
                        None
                    }
                })
                .flatten(),
        );
    }

    info!("Returning results");

    Ok(Results::vec_to_df(results))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rarity;
    fn rarity;
}
