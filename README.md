# PIGEON
`PIGEON` (**P**olygen**I**c **G**ene-**E**nvironment interacti**ON**) is a Python3-based command line tool for estimating polygenic gene-environment (GxE) interactions using GWIS (and GWAS) summary statistics.

The `PIGEON` software can be used to
1. estimate the GxE variance components (proportion of phenotypic variance explained by genome-wide GxE interactions).
2. perform hypothesis-free scans for PGSxE across many PGS (without the need to calculate PGS).

![](https://github.com/qlu-lab/PIGEON/blob/main/figure/PIGEON_Fig1.png)


## Manual

`PIGEON` can be downloaded via `git clone https://github.com/qlu-lab/PIGEON`

Please see the [wiki](https://github.com/qlu-lab/PIGEON/wiki) for the short tutorials describing the two basic functions (estimating GxE variance components and hypothesis-free scans for PGSxE), as well as the detailed manual of `PIGEON`.

## Version History
* Dec 12, 2022: Initial release.


## Citation

If you use `PIGEON`, please cite 

Miao, J., Song, G., Wu, Y., Hu, J., Wu, Y., Basu, S., Andrews, J. S., Schaumberg, K., Fletcher, J. M., Schmitz, L. L., & Lu, Q. (2022). [Reimagining Gene-Environment Interaction Analysis for Human Complex Traits](https://doi.org/10.1101/2022.12.11.519973 ). bioRxiv, 2022.2012.2011.519973. 

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jiacheng Miao at jmiao24@wisc.edu.

## "Birds" familial links
* [QUAIL](https://github.com/qlu-lab/QUAIL) (**QUA**ntile **I**ntegral **L**inear model) is a quantile regression-based framework to estimate genetic effects on the variance of quantitative traits.

## Acknowledgment

Part of the code is adapted from [LDSC](https://github.com/bulik/ldsc). We thank Dr. Bulik-Sullivan for sharing their code.
