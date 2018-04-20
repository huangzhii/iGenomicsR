options(stringsAsFactors = FALSE)
library("plyr")
setwd("C:/projects/iGenomicsR_BRCA/data/")

DB <- list()

# read RNA data
d <- read.csv("CPTAC_RNAseq_exp.csv")[,1:86]
d <- subset(d, X != "'?'" & X != "''" )
DB[["RNA"]] <- ddply(d, .(X), function(x){
  x <- x[,-1]
  if(nrow(x)>1){
    apply(x, 2, mean)
    } else {as.vector(as.matrix(x))}})
rownames(DB[["RNA"]]) <- DB[["RNA"]][,1]
DB[["RNA"]] <- DB[["RNA"]][,-1]
sample_name <- colnames(d)[-1]
sample_name <- substring(sample_name, 3, 12)
sample_name <- gsub("\\.", "-", sample_name)
sample_name <- paste("TCGA-", sample_name, sep="")
colnames(DB[["RNA"]]) <- sample_name
rownames(DB[["RNA"]]) <- gsub("\\'", "", rownames(DB[["RNA"]]))

# read protein data
d <- read.csv("CPTAC_protein_exp_sort.csv")
d <- subset(d, X != "'?'" & X != "''" )
DB[["Protein"]] <- ddply(d, .(X), function(x){
  x <- x[,-1]
  if(nrow(x)>1){
    apply(x, 2, mean)
  } else {as.vector(as.matrix(x))}})
rownames(DB[["Protein"]]) <- DB[["Protein"]][,1]
DB[["Protein"]] <- DB[["Protein"]][,-1]
sample_name <- colnames(d)[-1]
sample_name <- substring(sample_name, 3, 12)
sample_name <- gsub("\\.", "-", sample_name)
sample_name <- paste("TCGA-", sample_name, sep="")
colnames(DB[["Protein"]]) <- sample_name
rownames(DB[["Protein"]]) <- gsub("\\'", "", rownames(DB[["Protein"]]))


# read clinical data
d <- read.csv("Clinical_TCGA_all_04_02_13.csv")
sample_name <- d[,1]
sample_name <- paste(sample_name, "-01", sep="")
rownames(d) <- sample_name
DB[["Clinical"]] <- d[,-1]
DB[["Clinical_quan_lab"]] <- c()
DB[["Clinical_cat_lab"]] <- c("ajcc_neoplasm_disease_lymph_node_stage", "ajcc_neoplasm_disease_stage", "breast_carcinoma_estrogen_receptor_status", "breast_carcinoma_progesterone_receptor_status", "person_neoplasm_cancer_status")

# read survival data
d <- read.csv("Clinical_TCGA_survival.csv")
rownames(d) <- paste(d[,1], "-01", sep="")
d <- apply(d, c(1,2), function(x){gsub("\\ ", "", x)})
get_DiseaseFreeMonths <- function(x){
  if(any(x[c("person_neoplasm_cancer_status", "days_to_last_followup")] %in% c("[NotAvailable]", "[Unknown]") )){return(NA)} else {
    if(x["person_neoplasm_cancer_status"] == "TUMORFREE"){
      as.numeric(x["days_to_last_followup"]) / 30
    } else {
      0
    }
  }
}
get_OverallSurvivalMonths <- function(x){
    if(x["vital_status"] == "LIVING"){
      if(any(x[c("days_to_last_followup")] %in% c("[NotAvailable]", "[Unknown]") )){return(NA)} else {
        as.numeric(x["days_to_last_followup"]) / 30
      } }    else {
        as.numeric(x["days_to_death"]) / 30
      } 
}
d <- data.frame(DiseaseFreeMonths = apply(d, 1, get_DiseaseFreeMonths),
                DiseaseFreeStatus = d[,"person_neoplasm_cancer_status"],
                OverallSurvivalMonths = apply(d, 1, get_OverallSurvivalMonths),
                OverallSurvivalStatus = d[,"vital_status"])
DB[["Clinical"]] <- data.frame(DB[["Clinical"]], d[rownames(DB[["Clinical"]]),])

# download mutation data
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[24,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[11,1]
d <- list()
for(i in rownames(DB[["RNA"]])[3751:nrow(DB[["RNA"]])]){
  x <- getProfileData(mycgds,i,mygeneticprofile,mycaselist)
  if(ncol(x) == 1 ){
    d[[i]] <- x
  }
}
z <- t(as.data.frame(d))
DB[["Mutation_gene"]] <- apply(z, c(1,2), function(x){!is.na(x)})
colnames(DB[["Mutation_gene"]]) <- gsub("\\.", "-", colnames(DB[["Mutation_gene"]]))

# read image data
d1 <- read.csv("CPTAC_Image_epithelial.csv", header = TRUE, row.names = 1)
d2 <- read.csv("CPTAC_Image_stromal_40.csv", header = TRUE, row.names = 1)
all(colnames(d1) == colnames(d2))
d <- rbind(d1, d2)
sample_name <- colnames(d)
sample_name <- substring(sample_name, 3, 12)
sample_name <- gsub("\\.", "-", sample_name)
sample_name <- paste("TCGA-", sample_name, sep="")
colnames(d) <- sample_name
DB[["Image"]] <- d

# take intersection of samples
samples <- intersect(intersect(colnames(DB[["RNA"]]), colnames(DB[["Protein"]])), intersect(rownames(DB[["Clinical"]]), colnames(DB[["Mutation_gene"]])))
DB[["RNA"]] <- apply(DB[["RNA"]][,samples], c(1,2), as.numeric)
DB[["Protein"]] <- apply(DB[["Protein"]][,samples], c(1,2), as.numeric)
DB[["Clinical"]] <- DB[["Clinical"]][samples,]
DB[["Mutation_gene"]] <- apply(DB[["Mutation_gene"]][,samples], c(1,2), as.numeric)
DB[["Image"]] <- DB[["Image"]][,samples]

# test mutation profile association
test_res <- list()
all_genes <- rownames(DB[["Mutation_gene"]])[apply(DB[["Mutation_gene"]], 1, sum)>5]
for(i in 1:(length(all_genes)-1)){
  for(j in (i+1):length(all_genes)){
    gi <- all_genes[i]
    gj <- all_genes[j]
    ni <- sum(as.numeric(DB[["Mutation_gene"]][gi,]), na.rm=TRUE)
    nj <- sum(as.numeric(DB[["Mutation_gene"]][gj,]), na.rm=TRUE)
    nij <- sum(as.numeric(DB[["Mutation_gene"]][gi,]) == 1 & as.numeric(DB[["Mutation_gene"]][gj,]) == 1, na.rm=TRUE)
    if( ni >= 10 & nj >= 10){
      temp <- fisher.test(DB[["Mutation_gene"]][gi,], DB[["Mutation_gene"]][gj,])
      test_res[[paste(i,j)]] <- c(gi, gj, ni, nj, nij, temp$estimate, temp$p.value)
    }
  }
}
test_res <- t(as.data.frame(test_res))
colnames(test_res) <- c("Gene1", "Gene2", "nAltPatGene1", "nAltPatGene2", "nSharedAltPat", "oddsRatio", "pvalue")
test_res <- as.data.frame(test_res)
test_res[["adj_pvalue"]] <- p.adjust(test_res[,"pvalue",drop=TRUE], method="BH")
test_res <- test_res[order(as.numeric(test_res[,"pvalue"])),]
DB[["mutation_gene_test"]] <- subset(test_res, adj_pvalue < 0.05)



# Clinical_survival
DB[["Clinical_survival"]] <- c("DiseaseFreeMonths", "DiseaseFreeStatus", "OverallSurvivalMonths", "OverallSurvivalStatus")



save(DB, file="Kun_20170211.RData")
