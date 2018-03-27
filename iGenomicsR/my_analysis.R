
run_analysis <- function(dataType, patients){
  patients_list <- list()
  for(i in colnames(patients)){
    patients_list[[i]] <- setdiff(patients[,i], "")
  }
  if(dataType == "mutation"){
    res <- get_related_genes_by_mutation(patients_list)
  } else if(dataType == "rna"){
    res <- get_related_genes_by_rna(patients_list)
  }else if(dataType == "protein"){
    res <- get_related_genes_by_protein(patients_list)
  } else if(dataType == "clinical"){
    res <- get_related_genes_by_clinical(patients_list)
  }
  res
}


get_related_genes_by_mutation <- function(PatList){
  d <- DB[["Mutation_gene"]][,unlist(PatList)]
  res <- matrix(NA, ncol=4, nrow=nrow(d))
  rownames(res) <- rownames(d)
  colnames(res) <- c(paste("nAltPat_", names(PatList), sep=""), "pvalue", "adj_pvalue")
  for(i in names(PatList)){
    res[,paste("nAltPat_",i, sep="")] <- apply(d[,PatList[[i]]], 1, sum)
  }
  for(i in rownames(d)){
    if(length(unique(d[i,,drop=TRUE]))>1){
      res[i, "pvalue"] <- fisher.test(unlist(d[i,,drop=TRUE]), unlist(PatList) %in% PatList[[1]])$p.value
    }
  }
  res <- res[!is.na(res[, "pvalue"]),]
  res[, "adj_pvalue"] <- p.adjust(res[, "pvalue"], method="BH")
  res <- res[order(res[,"adj_pvalue"]),]
  if(sum(res[,"adj_pvalue"] <= 0.05) <= 10){
    res <- head(res, n=10)
  } else(
    res <- subset(as.data.frame(res), adj_pvalue <= 0.05)
  )
  res
}

get_related_genes_by_rna <- function(PatList){
  d <- DB[["RNA"]][,unlist(PatList)]
  res <- matrix(NA, ncol=4, nrow=nrow(d))
  rownames(res) <- rownames(d)
  colnames(res) <- c(paste("mean_", names(PatList), sep=""), "pvalue", "adj_pvalue")
  
  # no need to test if there is a group with patients < 5
  patients_no_exp <- apply(d, 2, function(x){all(is.na(x))})
  patients_no_exp <- names(patients_no_exp)[patients_no_exp]
  if(any(lapply(PatList, function(x){length(setdiff(x, patients_no_exp))}) < 5 )){
    res[1,"pvalue"] <-  "No enough patients with RNA data"
    return(res[1,,drop=FALSE])
  }
  
  for(i in names(PatList)){
    res[,paste("mean_",i, sep="")] <- apply(d[,PatList[[i]]], 1, mean, na.rm=TRUE)
  }
  for(i in rownames(d)){
    if( sum(! is.na( d[i, PatList[[1]]])) >= 5 & sum(! is.na( d[i, PatList[[2]]] )) >= 5 )
    res[i, "pvalue"] <- wilcox.test(d[i,PatList[[1]]], d[i,PatList[[2]]])$p.value
  }
  res <- res[!is.na(res[, "pvalue"]),]
  res[, "adj_pvalue"] <- p.adjust(res[, "pvalue"], method="BH")
  res <- res[order(res[,"adj_pvalue"]),]
  if(sum(res[,"adj_pvalue"] <= 0.05) <= 10){
    res <- head(res, n=10)
  } else(
    res <- subset(as.data.frame(res), adj_pvalue <= 0.05)
  )
  res
  
}

get_related_genes_by_protein <- function(PatList){
  d <- DB[["Protein"]][,unlist(PatList)]
  res <- matrix(NA, ncol=4, nrow=nrow(d))
  rownames(res) <- rownames(d)
  colnames(res) <- c(paste("mean_", names(PatList), sep=""), "pvalue", "adj_pvalue")
  
  # no need to test if there is a group with patients < 5
  patients_no_exp <- apply(d, 2, function(x){all(is.na(x))})
  patients_no_exp <- names(patients_no_exp)[patients_no_exp]
  if(any(lapply(PatList, function(x){length(setdiff(x, patients_no_exp))}) < 5 )){
    print(lapply(PatList, function(x){length(setdiff(x, patients_no_exp))}))
    res[1,"pvalue"] <-  "No enough patients with protein data"
    return(res[1,,drop=FALSE])
  }
  
  for(i in names(PatList)){
    res[,paste("mean_",i, sep="")] <- apply(d[,PatList[[i]]], 1, mean, na.rm=TRUE)
  }
  for(i in rownames(d)){
    if( sum(! is.na( d[i, PatList[[1]]])) >= 5 & sum(! is.na( d[i, PatList[[2]]] )) >= 5 )
    res[i, "pvalue"] <- wilcox.test(d[i,PatList[[1]]], d[i,PatList[[2]]])$p.value
  }
  res <- res[!is.na(res[, "pvalue"]),]
  res[, "adj_pvalue"] <- p.adjust(res[, "pvalue"], method="BH")
  res <- res[order(res[,"adj_pvalue"]),]
  if(sum(res[,"adj_pvalue"] <= 0.05) <= 10){
    res <- head(res, n=10)
  } else(
    res <- subset(as.data.frame(res), adj_pvalue <= 0.05)
  )
  res
  
}

get_related_genes_by_clinical <- function(PatList){
  
  # test categorical data
  d <- DB[["Clinical"]][unlist(PatList),]
  res <- list()
  for(i in DB[["Clinical_cat_lab"]]){
    if(all(unique(d[,i]) %in% c("YES", "NO"))){
      x <- d[,i] == "YES"
      y <- rownames(d) %in% PatList[[1]]
      if(length(setdiff(unique(x), NA))>1 & length(setdiff(unique(y)))>1 ){
        temp <- fisher.test(x, y)
        res[[i]] <- c(sum(d[,i] == "YES" & rownames(d) %in% PatList[[1]], na.rm = TRUE),
                      sum(d[,i] == "YES" & rownames(d) %in% PatList[[2]], na.rm = TRUE),
                      temp$p.value)
      }    
      } else {
        for(l in setdiff(unique(d[,i]), NA)){
          x <- d[,i] == l
          y <- rownames(d) %in% PatList[[1]]
          if(length(setdiff(unique(x), NA)) >1 & length(setdiff(unique(x), NA)) > 1){
            temp <- fisher.test(x, y)
            res[[paste(i, l, sep="_")]] <- c(sum(d[,i] == l & rownames(d) %in% PatList[[1]], na.rm = TRUE),
                                             sum(d[,i] == l & rownames(d) %in% PatList[[2]], na.rm = TRUE),
                                             temp$p.value)
          } 
        }
      }
  }
  
  # test quantitative data
  for(i in DB[["Clinical_quan_lab"]]){
    x <- as.numeric(d[PatList[[1]],i])
    y <- as.numeric(d[PatList[[2]],i])
    if(length(setdiff(x, NA)) >= 5 & length(setdiff(y, NA)) >= 5){
      temp <- wilcox.test(x, y)
      res[[i]] <- c(mean(d[PatList[[1]],i], na.rm=TRUE),
                    mean(d[PatList[[2]],i], na.rm=TRUE),
                    temp$p.value)
    }
  }
  
  res <- t(as.data.frame(res))
  res <- cbind(res, p.adjust(res[,3], method="BH"))
  colnames(res) <- c(paste("nPatIn_", names(PatList), sep=""), "pvalue", "adj_pvalue")
  res <- res[order(res[,"adj_pvalue"]),]
  if(sum(res[,"adj_pvalue"] <= 0.05) <= 10){
    res <- head(res, n=10)
  } else(
    res <- subset(as.data.frame(res), adj_pvalue <= 0.05)
  )
  res
  
}