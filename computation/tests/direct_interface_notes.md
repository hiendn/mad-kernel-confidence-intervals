# Direct interval interface

`computation/R/mad_interval.R` implements the four direct Gaussian-kernel rules. It returns unnormalised midpoint MAD, raw and reported root-n variance, sampling variance, standard error, additive and log intervals, bandwidth, cutpoints and local estimates. Reporting rescaling leaves raw cutpoints and bandwidth unchanged.

Failures return `success = FALSE`, a reason and unavailable interval fields; JSON has no nonfinite numeric tokens. Closed ECDF endpoints are deliberate. Exact reflection equivariance is checked only away from empirical endpoint ties.

Run `Rscript computation/tests/test_direct_interface.R`. The current 70 R checks include CLI execution, finite endpoints, scaling, affine transformations, endpoint ties, invalid data, zero MAD and numerical overflow/underflow. They do not fit GLDs or rerun Monte Carlo studies.
