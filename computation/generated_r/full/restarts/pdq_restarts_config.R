list(task = "pdq_restarts", seed = 20260917L, rng = "L'Ecuyer-CMRG/Inversion/Rejection", 
    seed_key_order = "master, experiment, zero-based population, n, zero-based replication", 
    seed_fold = "h <- (48271*h + key + 1) mod 2147483647", language = "R version 4.4.0 (2024-04-24)", 
    platform = "aarch64-apple-darwin20", source_md5 = c(`computation/R/nelder_mead.R` = "82d147eeae3c9eb2019cdf7d68762622", 
    `computation/R/mad_methods.R` = "242501b6cd9a960a7a711aee1a619453", 
    `computation/R/studies.R` = "fedbdf3b5e2c6dbcde65cbb78f76bb0c"
    ), args = list(implementation_md5 = "f447ac3b61b3726a959137c2c505e4c8", 
        ns = c(200L, 1000L, 5000L, 50000L), replications = 1000L, 
        results_dir = "computation/generated_r/full"))
