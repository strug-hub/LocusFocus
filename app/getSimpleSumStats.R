# Adapted from Fan Wang
# Script to obtain the simple sum P-values for a given set of GWAS p-values, and eQTL p-values for each tissue/gene pair
# Inputs: P_values_filename (GWAS p-values - for a set of SNPs - tab-separated, and all in one line)
#         ld_matrix_filename (the LD matrix filename for the set of SNPs input; the values per row must be tab-separated)
# Ouput: Returns a data.frame with the Simple Sum P-values, number of SNPs used and computation method (imhof or davies) used
# Example: getSimpleSumStats.R P_values_filename ld_matrix_filename

options(warn = -1)

# Check if required packages are installed:
list.of.packages <- c("argparser", "CompQuadForm", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(argparser, quietly = TRUE)
library(CompQuadForm, quietly = TRUE)
library(data.table, quietly = TRUE)
library(stringr, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(zeallot, quietly = TRUE)

######
# Parse arguments
######

p <- arg_parser("Calculate Simple Sum Statistic")
p <- add_argument(p, "P_values_filename", help = paste0(
  "Filename with GWAS and eQTL p-values - for a set of SNPs, each value tab-separated'\n'",
  "with 1st line being the GWAS p-values'\n'",
  "and each subsequent line is for eQTL p-values for each tissue/gene combination"
))
p <- add_argument(p, "ld_matrix_filename", help = paste0(
  "The LD matrix filename for the set of SNPs input;'\n'",
  "the values per row must be tab-separated; no header"
))
p <- add_argument(p, "--set_based_p", default = NULL, help = paste0(
  "For the Simple Sum method, a first-stage set-based Bonferroni p-value threshold'\n'",
  "is used for the set of secondary datasets with alpha 0.05'\n'",
  "(0.05 divided by the number of secondary datasets).'\n'",
  "Entering a value will override the default threshold."
))
p <- add_argument(p, "--outfilename", default = "SSPvalues.txt", help = "Output filename")
p <- add_argument(p, "--first_stage_only", flag = TRUE, help = "Whether to only perform and record first-stage set-based tests on both primary and secondary datasets.")
p <- add_argument(p, "--combine_lds", flag = TRUE, help = "Whether to combine LDs into one large sparse matrix.")
argv <- parse_args(p)
P_values_filename <- argv$P_values_filename
ld_matrix_filename <- argv$ld_matrix_filename
set_based_p <- argv$set_based_p
if (as.character(set_based_p) == "default") set_based_p <- NULL
outfilename <- argv$outfilename
first_stage_only <- argv$first_stage_only
combine_lds <- argv$combine_lds

session_data_dir <- dirname(ld_matrix_filename)
print(session_data_dir)
# test
# id <- "86bc6721-3b1e-45d6-bdf3-c12553d47e06"
# P_values_filename <- paste0('static/session_data/Pvalues-', id, '.txt')
# ld_matrix_filename <- paste0('static/session_data/ldmat-', id, '-001-002', '.txt')
# outfilename <- paste0('static/session_data/SSPvalues_setbasedtest-', id, '.txt')
# set_based_p <- "default"
# first_stage_only <- TRUE
# comnine_lds <- TRUE

ACONSTANT <- 6e-5

###############################################################################
############ FUNCTIONS
###############################################################################

set_based_test <- function(summary_stats, ld, num_genes, alpha = 0.05) {
  Z <- qnorm(summary_stats / 2)
  Zsq <- Z^2
  statistic <- sum(Zsq)
  m <- length(Zsq)
  eigenvalues <- eigen(ld)$values
  pv <- abs(imhof(statistic, eigenvalues)$Qq)
  if (is.null(set_based_p)) {
    if (pv < (alpha / num_genes)) {
      return(list(TRUE, pv))
    } else {
      return(list(FALSE, pv))
    }
  } else if (!is.na(as.numeric(set_based_p))) {
    if (pv < as.numeric(set_based_p)) {
      return(list(TRUE, pv))
    } else {
      return(list(FALSE, pv))
    }
  } else {
    stop(paste0("Provided set-based p-value (", set_based_p, ") is invalid."))
  }
}

get_p <- function(m, eigen_value, teststats, meth = "davies") {
  # l=length(eigen_value)
  # if(meth == 'davies'){
  #   pv<-abs(davies(teststats,eigen_value,h=rep(1,l),delta=rep(0,l))$Qq)
  # } else if(meth == 'imhof') {
  #   pv<-abs(imhof(teststats,eigen_value,h=rep(1,l),delta=rep(0,l))$Qq)
  # }
  if (meth == "davies") {
    pv <- davies(teststats, eigen_value)$Qq
  } else if (meth == "imhof") {
    pv <- imhof(teststats, eigen_value)$Qq
  }
  return(abs(pv))
}

get_a_diag <- function(eqtl_evid, m) {
  s <- sum(eqtl_evid)

  a_diag <- NULL
  if (s == 0 | s == m) {
    a_diag <- rep(1.0 / m, m)
  } else {
    t_bar <- mean(eqtl_evid)
    denom <- sum(eqtl_evid^2) - m * (t_bar^2)

    a_diag <- sapply(eqtl_evid, function(x) (x - t_bar) / denom)
  }

  return(a_diag)
}

get_eigenvalues <- function(eqtl_evid, ld.mat, m) {
  diag(ld.mat) <- diag(ld.mat) + ACONSTANT
  chol_Sigma <- chol(ld.mat)
  a_diag <- get_a_diag(eqtl_evid, m)

  matrix_A <- matrix(0, nrow = m, ncol = m)
  diag(matrix_A) <- a_diag

  matrix_mid <- chol_Sigma %*% matrix_A %*% t(chol_Sigma)
  eigenvalues <- eigen(matrix_mid)$values
  return(eigenvalues)
}

get_simple_sum_stats <- function(Zsq, eqtl_evid, m) {
  s <- sum(eqtl_evid)
  if (s == 0 | s == m) {
    return(mean(eqtl_evid))
  }

  reg <- lm(Zsq ~ eqtl_evid)
  SS <- summary(reg)$coefficients[2, 1]

  return(SS)
}

## if cut = 0, eQTL evidence would be -log10 transform of eQTL p-value;
## if cut < 0 (i.e. cut=0.05), eQTL evidence would be dichotomized eQTL p-value indicator by thresholds of eQTL p<cut.
get_eqtl_evid <- function(P, cut) {
  if (cut == 0) {
    covariate <- -log10(P)
  } else {
    covariate <- as.integer(P < cut)
  }
  return(covariate)
}


simple_sum_p <- function(P_gwas, P_eqtl, ld.mat, cut, m, meth = "davies") {
  ## need to match the GWAS SNP with the eQTL SNP and get m
  Z <- qnorm(P_gwas / 2)
  Zsq <- Z^2
  ## get eqtl evidence
  eqtl_evid <- get_eqtl_evid(P_eqtl, cut)
  # get Simple Sum statistic
  SS_stats <- get_simple_sum_stats(Zsq, eqtl_evid, m)
  ## get eigenvalues:
  eig_values <- get_eigenvalues(eqtl_evid, ld.mat, m)
  ## get Simple Sum p-values
  pv <- get_p(m, eig_values, SS_stats, meth = meth)
  return(pv)
}

drop_NA_from_LD <- function(P_mat, ld_mat) {
  if (!all(is.na(ld_mat))) {
    i <- 1
    while (any(is.na(ld_mat)) & i <= nrow(ld_mat)) {
      ldNA <- which(is.na(ld_mat[i, ]))
      if (!all(is.na(ldNA))) {
        ld_mat <- ld_mat[-ldNA, -ldNA]
        P_mat <- P_mat[, -ldNA, drop = FALSE]
      }
      i <- i + 1
    }
    return(list(P_mat=P_mat, ld_mat=ld_mat))
  } else {
    stop("LD matrix has all missing values")
  }
}

# Given a filename string, read all LD matrices and return a sparse, block diagonal matrix that combines all of them
# Filename string must be of the following format: `"ldmat-{UUID}-001-{end_index}.ld"`
# `{UUID}`: Unique identifier
# `{end_index}`: Total number of LDs; 3 digits with leading zeros
read_bdiag_LD <- function(ld_first_filename) {
  ld_filename_regex_pattern <- "(ldmat-.+-)([0-9]{3})-([0-9]{3})\\.txt$"

  matches <- str_match(ld_first_filename, ld_filename_regex_pattern)
  ld_prefix <- matches[2]
  start_index <- as.numeric(matches[3])
  end_index <- as.numeric(matches[4])
  ldmat_ <- fread(ld_first_filename, header = FALSE, stringsAsFactors = FALSE, na.strings = c("NaN", "nan", "NA", "-1"), sep = "\t")
  ldmat_ <- as.matrix(ldmat_)

  for (i in (start_index + 1):end_index) {
    if (i > end_index) {
      break
    }
    # load next LD and add it to our sparse matrix
    ld_filename <- sprintf("%s%03d-%03d.txt", ld_prefix, i, end_index)
    ldmat_next <- fread(file.path(session_data_dir, ld_filename), header = FALSE, stringsAsFactors = FALSE, na.strings = c("NaN", "nan", "NA", "-1"), sep = "\t")
    ldmat_next <- as.matrix(ldmat_next)
    ldmat_ <- bdiag(ldmat_, ldmat_next)
  }
  return(ldmat_)
}

############
# MAIN
############

# P-values returned can be negative and have the following meanings:
# -1: there was no eQTL data
# -2: fails the set_based_test(), so eQTL region is not significant after Bonferroni correction
# -3: could not compute the Simple Sum p-value; this is likely due to insufficient number of SNPs


### Load data

Pmat <- fread(P_values_filename, header = FALSE, stringsAsFactors = FALSE, na.strings = c("NaN", "nan", "NA", "-1"), sep = "\t")
# READ ldmat later
# filename = 'testdata/Pvalues.txt'
# Pmat <- fread(filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")
# filename = 'testdata/ldmat.txt'
# ldmat <- fread(filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")

# Columns for final returned result
Pss <- NULL
n <- NULL
comp_used <- NULL
first_stages <- NULL
first_stage_p <- NULL

Pmat <- as.matrix(Pmat)

if (nrow(Pmat) < 1) {
  stop("No secondary dataset P-values provided")
}

if (first_stage_only) {
  # only care about set based test result & P value
  # treat all lines in Pmat as if they came from .html
  num_lines <- nrow(Pmat)
  if (combine_lds) {
    ldmat <- read_bdiag_LD(ld_matrix_filename)
  } else {
    ldmat <- fread(ld_matrix_filename, header = FALSE, stringsAsFactors = FALSE, na.strings = c("NaN", "nan", "NA", "-1"), sep = "\t")
    ldmat <- as.matrix(ldmat)
  }
  for (i in 1:num_lines) {
    P_mat_i <- Pmat[i, ]
    ld_mat_i <- ldmat
    set_based_test_result <- "na"
    set_based_test_passed <- "na"
    set_based_test_p <- "na"
    # remove NA rows
    NArows <- which(is.na(P_mat_i))
    if (length(NArows) >= 1) {
      P_mat_i <- P_mat_i[-NArows]
      ld_mat_i <- ld_mat_i[-NArows, -NArows]
    }

    P_mat_i <- as.numeric(P_mat_i)

    if (length(P_mat_i) < 1) {
      first_stages <- c(first_stages, "na")
      first_stage_p <- c(first_stage_p, "na")
      next
    }

    # do pretest (set_based_test)
    t <- try({
      set_based_test_result <- set_based_test(P_mat_i, ld_mat_i, num_lines) # [1] is TRUE/FALSE, [2] is p value
      set_based_test_passed <- set_based_test_result[[1]]
      set_based_test_p <- set_based_test_result[[2]]
      first_stages <- c(first_stages, set_based_test_passed)
      first_stage_p <- c(first_stage_p, set_based_test_p)
    })
    if ("try-error" %in% class(t)) {
      print(t[1])
      first_stages <- c(first_stages, set_based_test_passed)
      first_stage_p <- c(first_stage_p, set_based_test_p)
    } # could not compute a SS p-value (SNPs not dense enough? can also get this if the LD matrix if not positive definite)
  }
  result <- data.frame(first_stages = first_stages, first_stage_p = first_stage_p)
} else {
  ldmat <- fread(ld_matrix_filename, header = FALSE, stringsAsFactors = FALSE, na.strings = c("NaN", "nan", "NA", "-1"), sep = "\t")
  ldmat <- as.matrix(ldmat)
  c(Pmat, ldmat) %<-% drop_NA_from_LD(Pmat, ldmat)

  # Normal simple sum here
  P_gwas <- Pmat[1, ]
  P_eqtl <- matrix(Pmat[2:nrow(Pmat), ], nrow = nrow(Pmat) - 1, ncol = ncol(Pmat))

  num_genes <- nrow(P_eqtl)
  for (i in 1:num_genes) {
    tempmat <- cbind(P_gwas, P_eqtl[i, ])
    ld_mat_i <- ldmat
    set_based_test_result <- "na"
    set_based_test_passed <- "na"
    set_based_test_p <- "na"

    # Remove NA rows
    NArows <- which(is.na(tempmat[, 1]) | is.na(tempmat[, 2]))
    if (length(NArows) >= 1) {
      tempmat <- tempmat[-NArows, ]
      ld_mat_i <- ld_mat_i[-NArows, -NArows]
    }

    P_gwas_i <- as.numeric(tempmat[, 1])
    P_eqtl_i <- as.numeric(tempmat[, 2])

    # Count SNPs
    snp_count <- nrow(tempmat)
    n <- c(n, snp_count)

    if (snp_count < 1) {
      Pss <- c(Pss, -1) # no eQTL data
      comp_used <- c(comp_used, "na")
      first_stages <- c(first_stages, "na")
      first_stage_p <- c(first_stage_p, set_based_test_p)
      next
    }

    # do pretest (set_based_test)
    t <- try({
      # unpack list
      c(set_based_test_passed, set_based_test_p) %<-% set_based_test(P_eqtl_i, ld_mat_i, num_genes) # [1] is TRUE/FALSE, [2] is p value
      if (isTRUE(set_based_test_passed)) {
        # if(TRUE) {
        P <- simple_sum_p(P_gwas = P_gwas_i, P_eqtl = P_eqtl_i, ld.mat = ld_mat_i, cut = 0, m = snp_count, meth = "davies")
        if (P == 0 | P < 0) {
          P <- simple_sum_p(P_gwas = P_gwas_i, P_eqtl = P_eqtl_i, ld.mat = ld_mat_i, cut = 0, m = snp_count, meth = "imhof")
          comp_used <- c(comp_used, "imhof")
          Pss <- c(Pss, P)
        } else {
          comp_used <- c(comp_used, "davies")
          Pss <- c(Pss, P)
        }
      } else {
        # always happens if first_stage_only, but we only look at first_stages afterwards anyways
        Pss <- c(Pss, -2) # not significant eQTL as per set-based test
        comp_used <- c(comp_used, "na")
      }
      first_stages <- c(first_stages, set_based_test_passed)
      first_stage_p <- c(first_stage_p, set_based_test_p)
    })
    if ("try-error" %in% class(t)) {
      print(t[1])
      Pss <- c(Pss, -3)
      comp_used <- c(comp_used, "na")
      first_stages <- c(first_stages, set_based_test_passed)
      first_stage_p <- c(first_stage_p, set_based_test_p)
    } # could not compute a SS p-value (SNPs not dense enough? can also get this if the LD matrix if not positive definite)
  }

  # print(Pss)
  # print(first_stages)
  # print(first_stage_p)

  result <- data.frame(Pss = Pss, n = n, comp_used = comp_used, first_stages = first_stages, first_stage_p = first_stage_p)
}


sessionid <- gsub(".txt", "", gsub("Pvalues-", "", P_values_filename))
write.table(result, outfilename, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
