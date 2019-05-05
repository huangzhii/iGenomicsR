options(stringsAsFactors = FALSE)

# ensembl to gene symbol
ensembl2symbol <- read.table("ensembl2symbol.txt", sep="\t", header=TRUE)


save(ensembl2symbol, file="gene_id_db.RData")
