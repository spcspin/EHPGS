# EHPGS
EHPGS is a  GS-based approach to identify potential parental lines and superior hybrid combinations from a breeding population composed of hybrids produced by a half diallel mating design.

## Installation

The EHPGS can be installed from GitHub:

``` r
library(devtools)
install_github("spcspin/EHPGS", dependencies = TRUE, force = TRUE)
```

## Main functions

- `kinship`: Function for calculate relationship matrix.
- `BGS`: Function for Bayesian Gibbs sampling algorithm.
- `EHPGS`: Function for evaluation of hybrid performance in plant breeding via genomic selection.

## Example dataset
The test dataset provided in this package is the pumpkin dataset which was published by Wu et al. (2019) <doi.org/10.3835/plantgenome2018.10.0082>.

## Authors

- Szu-Ping Chen
    - Author, maintainer
    - E-mail: R09621108@ntu.edu.tw
    - Department of Agronomy, National Taiwan University, Taipei, Taiwan
- Chih-Wei Tung
    - Author
    - E-mail: chihweitung@ntu.edu.tw
    - Department of Agronomy, National Taiwan University, Taipei, Taiwan
- Pei-Hsien Wang
    - Author
    - E-mail: R09621206@ntu.edu.tw
    - Department of Agronomy, National Taiwan University, Taipei, Taiwan    
- Chen-Tuo Liao
    - Author, thesis advisor
    - E-mail: ctliao@ntu.edu.tw
    - Department of Agronomy, National Taiwan University, Taipei, Taiwan

## Tutorial 



