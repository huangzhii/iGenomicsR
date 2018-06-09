quantile_rm_na <- function(x , probs){
  quantile(x , probs, na.rm=TRUE)
}
gene_filter <- function(exp_table, method, cutoff){
  res <- list()
  if(method == "var"){
    var_vector <- apply(exp_table, 1, var)
    res[["cutoff"]] <- paste("var cutoff:", quantile_rm_na(var_vector, probs=cutoff))
    filtered_table <- exp_table[var_vector > quantile_rm_na(var_vector, probs=cutoff[1]) & var_vector < quantile_rm_na(var_vector, probs=cutoff[2]), ]
  } else if (method == "cv"){
    cv_vector <- apply(exp_table, 1, function(x){ sd(x) / mean(x)})
    res[["cutoff"]] <- paste("cv cutoff:", quantile_rm_na(cv_vector, probs=cutoff))
    filtered_table <- exp_table[cv_vector > quantile_rm_na(cv_vector, probs=cutoff[1]) & cv_vector < quantile_rm_na(cv_vector, probs=cutoff[2]), ]
  } else if (method == "sumcov") {
    sumcov_table <- cov(t(exp_table))
    sumcov_vector <- apply(sumcov_table, 1, sum)
    res[["cutoff"]] <- paste("sumcov cutoff:", quantile_rm_na(sumcov_vector, probs=cutoff))
    filtered_table <- exp_table[sumcov_vector > quantile_rm_na(sumcov_vector, probs=cutoff[1]) & sumcov_vector < quantile_rm_na(sumcov_vector, probs=cutoff[2]), ]
  } else if (method == "sumcor") {
    sumcor_table <- cov(t(exp_table))
    sumcor_vector <- apply(sumcor_table, 1, sum)
    res[["cutoff"]] <- paste("sumcor cutoff:", quantile_rm_na(sumcor_vector, probs=cutoff))
    filtered_table <- exp_table[sumcor_vector > quantile_rm_na(sumcor_vector, probs=cutoff[1]) & sumcor_vector < quantile_rm_na(sumcor_vector, probs=cutoff[2]), ]
  } else if (method == "maxExp"){
    maxexp_vector <- apply(exp_table, 1, max)
    res[["cutoff"]] <- paste("maxexp cutoff:", quantile_rm_na(maxexp_vector, probs=cutoff))
    filtered_table <- exp_table[maxexp_vector > quantile_rm_na(maxexp_vector, probs=cutoff[1]) & maxexp_vector < quantile_rm_na(maxexp_vector, probs=cutoff[2]), ]
  } else if (method == "expCells"){
    nExpCells <- apply(exp_table, 1, function(x){sum(x >= cutoff)})
    res[["cutoff"]] <- paste("expCells cutoff: >=", nExpCell, "cells with FPKM >=", cutoff)
    filtered_table <- exp_table[nExpCells >= nExpCell, ]
  }
  res[["genes"]] <- rownames(filtered_table) 
  res
}