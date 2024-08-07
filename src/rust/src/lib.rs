use extendr_api::prelude::*;
use lmutils::{get_r2s, IntoMatrix, Transform};
use log::{debug, error, info};
use rayon::prelude::*;

struct Results {
    pheno_file: String,
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
            pheno_file = results
                .iter()
                .map(|x| x.pheno_file.clone())
                .collect::<Vec<_>>(),
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

    let genes = (0..=22)
        .map(|chr| {
            let chr_dir = if chr == 0 {
                dir.to_owned()
            } else {
                dir.join(format!("chr_{:02}", chr))
            };
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

    let pheno_files = phenos.iter().map(|x| x.to_string()).collect::<Vec<_>>();
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
    let phenos = phenos
        .into_par_iter()
        .map(|x| x.to_owned().unwrap())
        .collect::<Vec<_>>();
    let nrows = phenos[0].rows();
    for pheno in phenos.iter() {
        if pheno.rows() != nrows {
            return Err(Error::from("Phenotypes must have the same number of rows"));
        }
        if pheno.data().par_iter().any(|x| x.is_nan()) {
            return Err(Error::from("Phenotypes must not contain NaN values"));
        }
    }
    let phenos = phenos
        .into_par_iter()
        .map(|mut x| {
            x.remove_column_by_name_if_exists("eid");
            x.remove_column_by_name_if_exists("IID");
            x.into_matrix()
                .standardization()
                .transform()
                .unwrap()
                .to_owned()
                .unwrap()
        })
        .collect::<Vec<_>>();
    let pheno_norm = phenos.iter().map(|x| x.as_mat_ref()).collect::<Vec<_>>();
    let traits = phenos
        .iter()
        .map(|x| {
            x.colnames()
                .map(|x| x.to_vec())
                .unwrap_or_else(|| (1..=x.cols()).map(|x| x.to_string()).collect::<Vec<_>>())
        })
        .collect::<Vec<_>>();

    info!("Calculating RARity");

    let genes = std::sync::Mutex::new(
        genes
            .into_par_iter()
            .enumerate()
            .flat_map(|(chr, g)| {
                let chr = chr + 1;
                g.into_par_iter().map(move |g| (chr, g.0, g.1))
            })
            .collect::<Vec<_>>(),
    );

    let blocks_per_chunk = std::env::var("RARITY_BLOCKS_PER_CHUNK")
        .unwrap_or("16".to_string())
        .parse::<usize>()
        .unwrap()
        .clamp(1, genes.lock().unwrap().len());

    let results = std::sync::Mutex::new(vec![]);

    let min_sum = std::env::var("RARITY_MIN_SUM")
        .unwrap_or("2.0".to_string())
        .parse::<f64>()
        .unwrap();

    std::thread::scope(|s| {
        for _ in 0..blocks_per_chunk {
            s.spawn(|| loop {
                let mut genes = genes.lock().unwrap();
                let gene = genes.pop();
                drop(genes);
                if let Some((chr, path, mat)) = gene {
                    rayon::scope(|s| {
                        s.spawn(|_| {
                            let gene = path.file_name().unwrap().to_str().unwrap();
                            info!("Processing gene {}", gene);
                            if std::fs::metadata(&path).is_ok() {
                                let mut block: lmutils::OwnedMatrix<f64> = mat.to_owned().unwrap();
                                if block.rows() != nrows {
                                    error!("Block {gene} has different number of rows than phenotypes, expected {}, found {}", nrows, block.rows());
                                    return;
                                }
                                debug!("Removing eid column");
                                block.remove_column_by_name_if_exists("eid");
                                block.remove_column_by_name_if_exists("IID");
                                let block = block.into_matrix();
                                debug!("Normalizing block");
                                let block = block
                                    .nan_to_mean()
                                    .min_sum(min_sum)
                                    .standardization()
                                    .transform()
                                    .unwrap();
                                let block = block.as_mat_ref().unwrap();
                                if block.ncols() == 0 {
                                    return;
                                }
                                let res = pheno_norm
                                    .par_iter()
                                    .zip(&pheno_files)
                                    .zip(&traits)
                                    .flat_map(|((pheno_norm, pheno_file), traits)| {
                                        info!("Calculating R2 for gene {}", gene);
                                        let r2s = get_r2s(block, *pheno_norm);
                                        let nb_individuals = block.nrows();
                                        let nb_rvs = block.ncols();
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
                                                pheno_file: pheno_file.to_string(),
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
                                results.lock().unwrap().extend(res);
                                info!("Processed gene {}", gene);
                            } else {
                                info!("Gene {} not found", gene);
                            };
                        })
                    });
                } else {
                    break;
                }
            });
        }
    });

    info!("Returning results");

    Ok(Results::vec_to_df(results.into_inner().unwrap()))
}

/// Set the log level.
/// `level` is the log level.
/// @export
#[extendr]
pub fn set_log_level(level: &str) {
    std::env::set_var("RARITY_LOG", level);
}

/// Set the number of blocks per chunk.
/// `blocks_per_chunk` is the number of blocks per chunk.
/// @export
#[extendr]
pub fn set_blocks_per_chunk(blocks_per_chunk: u32) {
    std::env::set_var("RARITY_BLOCKS_PER_CHUNK", blocks_per_chunk.to_string());
}

/// Set the number of threads.
/// `num_threads` is the number of threads.
/// @export
#[extendr]
pub fn set_num_threads(num_threads: u32) {
    std::env::set_var("RAYON_NUM_THREADS", num_threads.to_string());
}

/// Set the min sum for the RARity analysis.
/// `min_sum` is the min sum.
/// @export
#[extendr]
pub fn set_min_sum(min_sum: f64) {
    std::env::set_var("RARITY_MIN_SUM", min_sum.to_string());
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rarity;
    fn rarity;
    fn set_log_level;
    fn set_blocks_per_chunk;
    fn set_num_threads;
    fn set_min_sum;
}
