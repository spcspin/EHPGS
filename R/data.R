#' Genotype information of parental lines
#'
#' Genotype information of parental lines. This data was published by Guo et al. (2019) <doi:10.1016/j.molp.2018.12.022>
#'
#' @format A numeric matrix with 24 rows (parental lines) and 46,134 columns (SNPs).
#'
#' @examples
#' data(parent.geno)
"parent.geno"

#' Genotype information of training set.
#'
#' Genotype information of F1s used to construct a prediction model. This data was published by Guo et al. (2019) <doi:10.1016/j.molp.2018.12.022>
#'
#' @format A numeric matrix with 276 rows (F1s of half diallel design) and 46,134 columns (SNPs).
#'
#' @examples
#' data(train.geno)
"train.geno"

#' Phenotype
#'
#' Grain yield (Mg/ha)  of samples. This data was published by Guo et al. (2019) <doi:10.1016/j.molp.2018.12.022>
#'
#' @format A numeric vector of grain yield values which were adjusted from the variability of the different locations and years.
#'
#' @examples
#' data(BLUP.pheno)
"BLUP.pheno"

#' Phenotype
#'
#' Grain yield (Mg/ha)  of samples at Clayton, NC. This data was published by Guo et al. (2019) <doi:10.1016/j.molp.2018.12.022>
#'
#' @format A numeric vector of grain yield values at Clayton, NC.
#'
#' @examples
#' data(NC.pheno)
"NC.pheno"

#' Phenotype
#'
#' Grain yield (Mg/ha)  of samples at Columbia, MO. This data was published by Guo et al. (2019) <doi:10.1016/j.molp.2018.12.022>
#'
#' @format A numeric vector of grain yield values at Columbia, MO.
#'
#' @examples
#' data(MO.pheno)
"MO.pheno"
