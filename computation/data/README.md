# Input data

`galaxies.csv` contains the 82 galaxy velocities in km/s used in the manuscript illustration. The study driver converts these to thousands of km/s. The data values are unchanged from the manuscript's supplied input.

`benchmark/` contains the eight fixed vectors used for descriptive timings: normal, lognormal, shoulder and bimodal populations at n=200 and n=1000. These vectors support timing reproduction only; the Monte Carlo drivers generate samples from their documented seeds.

The repository MANIFEST_SHA256.json records the hashes of these input files.
