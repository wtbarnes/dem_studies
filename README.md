# dem_studies
Comparison of various differential emission measure (DEM) codes. 

## Methods
These various methods include, but are not limited to,

* MCMC method of [Kashyap and Drake (1998)](https://ui.adsabs.harvard.edu/#abs/1998ApJ...503..450K/abstract) -- see POA package
* Regularized DEM method of [Hanna et al. (2012)](https://ui.adsabs.harvard.edu/#abs/2012A&A...539A.146H/abstract) -- [code](https://github.com/ianan/demreg)
* Bayesian inference method of [Warren et al. (2016)](https://ui.adsabs.harvard.edu/#abs/2016arXiv161005972W/abstract) -- no code publicly available
* XRT DEM (see SSW)
* AIA DEM (see SSW)
* CHIANTI DEM (see SSW/Chianti)

## Helpful DEM Papers
Some helpful papers on DEM inversion/interpretation

* [Aschwanden et al. (2015)](https://ui.adsabs.harvard.edu/#abs/2015SoPh..290.2733A/abstract) -- benchmarks of various methods
* [Landi et al. (2012)](https://ui.adsabs.harvard.edu/#abs/2012A&A...538A.111L/abstract) -- DEM reconstruction summary

## Roadmap

0. Choose heating parameters and AR
1. Model loops and emission for the needed EIS channels
2. Choose region/time to evaluate DEM at
3. Extract intensities for various lines
4. Apply inversion method
