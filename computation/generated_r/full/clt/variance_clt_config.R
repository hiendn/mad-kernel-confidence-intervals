list(task = "variance_clt", seed = 20260917L, rng = "L'Ecuyer-CMRG/Inversion/Rejection", 
    seed_key_order = "master, experiment, zero-based population, n, zero-based replication", 
    seed_fold = "h <- (48271*h + key + 1) mod 2147483647", language = "R version 4.4.0 (2024-04-24)", 
    platform = "aarch64-apple-darwin20", source_md5 = c(`computation/R/nelder_mead.R` = "82d147eeae3c9eb2019cdf7d68762622", 
    `computation/R/mad_methods.R` = "242501b6cd9a960a7a711aee1a619453", 
    `computation/R/studies.R` = "c480dd6988e46d9dd06c6c4075b8616b"
    ), args = list(replications = 3000L, ns = c(200L, 1000L, 
    5000L), pop_indices = c(0L, 2L)))
