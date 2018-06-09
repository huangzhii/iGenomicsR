options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)
library(cowplot)

bin2dec <- function(x){
  new_x <- x
  for(i in 1:(length(x)-1)){
    new_x[1:i] <- new_x[1:i] * 2
  }
  sum(new_x)
}

fun_order_samples <- function(mutation_genes, rna_genes, protein_genes, clinical_lab, order_by, selected_samples,clust_para=list(method="hc"), order_clin_feature, image_features){
  res <- list()
  
  if(order_by == "mutation"){
    geneScore <- apply(DB[["Mutation_gene"]][mutation_genes,selected_samples,drop=FALSE], 1, sum, na.rm=TRUE)
    mutation_genes <- names(sort(geneScore, decreasing = TRUE))
    sampleScore <- apply(DB[["Mutation_gene"]][mutation_genes,selected_samples,drop=FALSE], 2, bin2dec)
    sample_order <- order(sampleScore, decreasing = TRUE)
    res[["ordered_samples"]] <- selected_samples[sample_order]
  } else if (order_by == "rna"){
    if(clust_para[["method"]]=="hc"){
      res[["hc"]] <- hclust(dist(t(DB[["RNA"]][rna_genes,selected_samples,drop=FALSE])), method="average")
      res[["ordered_samples"]] <- selected_samples[res[["hc"]]$order]
    } else {
      cl_temp <- kmeans(t(DB[["RNA"]][rna_genes,selected_samples,drop=FALSE]), clust_para[["k"]])
      res[["ordered_samples"]] <- selected_samples[order(cl_temp$cluster)]
      res[["cluster"]] <- cl_temp$cluster
    }
  } else if (order_by == "protein"){
    if(length(protein_genes)==1){
      res[["ordered_samples"]] <- selected_samples[order(DB[["Protein"]][protein_genes,selected_samples])]
    } else {
      temp_exp <- DB[["Protein"]][protein_genes,selected_samples,drop=FALSE]
      temp_exp <- apply(temp_exp, c(1,2), function(x){if(is.na(x)){0} else {x}})
      res[["hc"]] <- hclust(dist(t(temp_exp)), method="average")
      res[["ordered_samples"]] <- selected_samples[res[["hc"]]$order]
    }
  } else if (order_by == "image"){
    if(length(image_features)==1){
      res[["ordered_samples"]] <- selected_samples[order(DB[["Image"]][image_features,selected_samples])]
    } else {
      temp_exp <- DB[["Image"]][image_features,selected_samples,drop=FALSE]
      temp_exp <- apply(temp_exp, c(1,2), function(x){if(is.na(x)){0} else {x}})
      res[["hc"]] <- hclust(dist(t(temp_exp)), method="average")
      res[["ordered_samples"]] <- selected_samples[res[["hc"]]$order]
    } } else {
      if(order_clin_feature %in% DB[["Clinical_quan_lab"]]){
        to_order <- as.numeric(DB[["Clinical"]][selected_samples, order_clin_feature, drop=TRUE])
      } else {
        to_order <- DB[["Clinical"]][selected_samples, order_clin_feature, drop=TRUE]
      }
      res[["ordered_samples"]] <- selected_samples[order(to_order)]
    }
  res
}

fun_order_Genes <- function(genes, order_by, selected_samples){
  t <- try(
    if(order_by == "Mutation_gene"){
      geneScore <- apply(DB[["Mutation_gene"]][genes,selected_samples,drop=FALSE], 1, sum, na.rm=TRUE)
      ordered_genes <- names(sort(geneScore, decreasing = FALSE))
    } else if (order_by == "RNA"){
      gene_order <- hclust(as.dist(1-cor(t(DB[["RNA"]][genes,selected_samples,drop=FALSE]))), method="average")$order
      ordered_genes <- genes[gene_order]
    } else if (order_by == "Image"){
      gene_order <- hclust(as.dist(1-cor(t(DB[["Image"]][genes,selected_samples,drop=FALSE]))), method="average")$order
      ordered_genes <- genes[gene_order]
    } else if (order_by == "Protein"){
      gene_order <- hclust(as.dist(1-cor(t(DB[["Protein"]][genes,selected_samples,drop=FALSE]))), method="average")$order
      ordered_genes <- genes[gene_order]
    })
  if("try-error" %in% class(t)) {
    removeModal()
    not_in_list = genes[!(genes %in% rownames(DB[[order_by]]))]
    not_in_list = paste(not_in_list, collapse =', ')
    print(sprintf("Genes [ %s ] not found in DB.", not_in_list))
    showModal(modalDialog(
      title = "Input genes/features not found.", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p(sprintf("Error in fun_order_Genes: Genes [ %s ] are not belong into the database you previously defined. Please remove them from your input.", not_in_list)),
          style = "margin: auto; text-align: center"
      )
    ))
    return()
  }
  else {
    return(ordered_genes)
  }
  
}

########################################################################
# mutation based 
########################################################################
my_heatmap_mutation <- function(mutation_genes, rna_genes, protein_genes, clinical_lab,
                                order_by, rna_size=1, protein_size = 1, selected_samples=colnames(DB[["RNA"]]),
                                clust_para=list(method="hc"), order_clin_feature, image_features=vector(),
                                show.RNA.name = 1, show.protein.name = 1, sort.mutation=F, sort.rna=F,
                                sort.image=F, sort.protein=F){
  # remove duplicated inputs
  if(length(unique(mutation_genes)) != length(mutation_genes)){
    removeModal()
    showModal(modalDialog(
      title = "Warning", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p("Exists duplicated mutation genes. We help you removed it."),
          style = "margin: auto; text-align: center"
      )
    ))
    mutation_genes <- unique(mutation_genes)
  }
  if(length(unique(rna_genes)) != length(rna_genes)){
    removeModal()
    showModal(modalDialog(
      title = "Warning", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p("Exists duplicated RNA genes. We help you removed it."),
          style = "margin: auto; text-align: center"
      )
    ))
    rna_genes <- unique(rna_genes)
  }
  if(length(unique(protein_genes)) != length(protein_genes)){
    removeModal()
    showModal(modalDialog(
      title = "Warning", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p("Exists duplicated protein genes. We help you removed it."),
          style = "margin: auto; text-align: center"
      )
    ))
    protein_genes <- unique(protein_genes)
  }
  if(length(unique(image_features)) != length(image_features)){
    removeModal()
    showModal(modalDialog(
      title = "Warning", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p("Exists duplicated image features. We help you removed it."),
          style = "margin: auto; text-align: center"
      )
    ))
    image_features <- unique(image_features)
  }
  
  
  sample_order_res <- fun_order_samples(mutation_genes, rna_genes, protein_genes, clinical_lab, order_by, selected_samples, clust_para,order_clin_feature, image_features)
  ordered_samples <- sample_order_res[["ordered_samples"]] 
  
  if(length(mutation_genes) > 1 & sort.mutation){
    mutation_genes <- fun_order_Genes(mutation_genes, "Mutation_gene", ordered_samples)
  }
  if(length(rna_genes) > 2 & sort.rna){
    rna_genes <- fun_order_Genes(rna_genes, "RNA", ordered_samples)
  }
  if(length(image_features) > 2 & sort.image){
    image_features <- fun_order_Genes(image_features, "Image", ordered_samples)
  }
  if(length(protein_genes) > 2 & sort.protein){
    protein_genes <- fun_order_Genes(protein_genes, "Protein", ordered_samples)
  }
  
  PL <- list()
  plot_heights <- vector()
  my_theme <- theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
                    axis.text=element_text(size=12),
                    plot.margin=unit(c(0,0,-0.3,0), "cm"),
                    # legend.key.size=unit(0.5, "cm"),
                    legend.direction = "horizontal",
                    legend.position = "right", #("none", "left", "right", "bottom", "top", or two-element numeric vector)
                    legend.box = "horizontal",
                    legend.title = element_blank()
  )
  
  res <- list()
  # heatmap for image feature
  if(length(image_features)>0){
    image_features = rev(image_features) # This is for descending order of names in y axis
    t <- try(image <- DB[["Image"]][image_features,ordered_samples,drop=FALSE])
    if("try-error" %in% class(t)) {
      removeModal()
      not_in_list = image_features[!(image_features %in% rownames(DB[["Image"]]))]
      not_in_list = paste(not_in_list, collapse =', ')
      print(sprintf("genes [ %s ] not found in DB.", not_in_list))
      showModal(modalDialog(
        title = "Input genes/features not found.", footer = modalButton("OK"), easyClose = TRUE,
        div(class = "busy",
            p(sprintf("genes [ %s ] are not belong into the database you previously defined. Please remove them from your input.", not_in_list)),
            style = "margin: auto; text-align: center"
        )
      ))
      return()
    }
    res[["image_data"]] <- t(image)
    image <- t(apply(image, 1, scale))
    colnames(image) <- ordered_samples
    image <- melt(image)
    image[["Var2"]] <- factor(image[["Var2"]], levels=ordered_samples, ordered = TRUE)
    image[["Var1"]] <- factor(image[["Var1"]], levels=image_features, ordered = TRUE)
    # save(image, file = "~/Desktop/image.Rdata")
    image.value.tmp = image$value
    image.value.tmp[is.na(image.value.tmp)] <- 0
    PL[["image_plot"]] <- ggplot(image,  aes(Var2, Var1)) + 
      geom_tile(aes(fill = value)) + 
      scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(image.value.tmp),max(image.value.tmp),(max(image.value.tmp)-min(image.value.tmp))/5), digits=1)) +
      guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical") +
      theme_bw() + my_theme +
      labs(x="", y="")
    plot_heights <- c(plot_heights, length(image_features))
  }
  # heatmap for RNA
  if(length(rna_genes) >0 ){
    print("--- 1 ---")
    rna_genes = rev(rna_genes) # This is for descending order of names in y axis
    t <- try(res[["rna_data"]] <- t(DB[["RNA"]][rna_genes,ordered_samples,drop=FALSE]))
    if("try-error" %in% class(t)) {
      removeModal()
      not_in_list = rna_genes[!(rna_genes %in% rownames(DB[["RNA"]]))]
      not_in_list = paste(not_in_list, collapse =', ')
      print(sprintf("genes [ %s ] not found in DB.", not_in_list))
      showModal(modalDialog(
        title = "Input genes/features not found.", footer = modalButton("OK"), easyClose = TRUE,
        div(class = "busy",
            p(sprintf("genes [ %s ] are not belong into the database you previously defined. Please remove them from your input.", not_in_list)),
            style = "margin: auto; text-align: center"
        )
      ))
      return()
    }
    # colnames(res[["rna_data"]]) <- paste("rna_", colnames(res[["rna_data"]]), sep="")
    rna <- t(apply(DB[["RNA"]][rna_genes,ordered_samples,drop=FALSE], 1, scale))
    colnames(rna) <- ordered_samples
    rna <- melt(rna)
    rna[["Var2"]] <- factor(rna[["Var2"]], levels=ordered_samples, ordered = TRUE)
    rna[["Var1"]] <- factor(rna[["Var1"]], levels=rna_genes, ordered = TRUE)
    rna[["value"]][ rna[["value"]] > 2] <- 2
    rna[["value"]][ rna[["value"]] < -2] <- -2
    
    if(rna_size == 1){
      print("--- 2 ---")
      PL[["rna_plot"]] <- ggplot(rna,  aes(Var2, Var1)) + 
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(rna$value),max(rna$value),(max(rna$value)-min(rna$value))/5), digits=1)) +
        guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical") +
        theme_bw() + my_theme + labs(x="", y="")
      print("length of rna genes:")
      print(length(rna_genes))
      plot_heights <- c(plot_heights, length(rna_genes))
      
      if(clust_para[["method"]]=="km"){
        print("--- 3 ---")
        d <- data.frame(patient_id=ordered_samples, cluster=paste("Group ",sort(sample_order_res[["cluster"]]), sep=""))
        res[["rna_group"]] <- d[,"cluster",drop=FALSE]
        d[["patient_id"]] <- factor(d[["patient_id"]], levels=ordered_samples, ordered = TRUE)
        PL[["sample_group"]] <- ggplot(d, aes(patient_id, "cluster")) + 
          geom_tile(aes(fill = cluster), colour = "white") +
          theme_bw() + labs(x="", y="") + my_theme
        plot_heights <- c(plot_heights, 1)
      }
    } else {
      print("--- 4 ---")
      if(show.RNA.name == 1){
        PL[["rna_plot"]] <- ggplot(rna,  aes(Var2, Var1)) + 
          geom_tile(aes(fill = value)) + 
          scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(rna$value),max(rna$value),(max(rna$value)-min(rna$value))/5), digits=1)) +
          guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical")  +
          theme_bw() + labs(x="", y="") + my_theme + theme(plot.margin=unit(c(0,0,0,0), "cm"))
      }
      else{
        PL[["rna_plot"]] <- ggplot(rna,  aes(Var2, Var1)) + 
          geom_tile(aes(fill = value)) + 
          scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(rna$value),max(rna$value),(max(rna$value)-min(rna$value))/5), digits=1)) +
          guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical")  +
          theme_bw() + labs(x="", y="") + my_theme + theme(plot.margin=unit(c(0,0,0,0), "cm")) +
          theme(axis.text.y = element_blank())
      }
      plot_heights <- c(plot_heights, length(rna_genes) * rna_size)
      
      d <- data.frame(patient_id=ordered_samples, divide_patients_by_cuttree(sample_order_res[["hc"]])[ordered_samples,])
      res[["rna_group"]] <- d[,2:ncol(d),drop=FALSE]
      d <- melt(d)
      d[["patient_id"]] <- factor(d[["patient_id"]], levels=ordered_samples, ordered = TRUE)
      d[["value"]] <- as.factor(d[["value"]])
      PL[["sample_group"]] <- ggplot(d, aes(patient_id, variable)) + 
        geom_tile(aes(fill = value), colour = "white") +
        theme_bw() + labs(x="", y="") + my_theme
      plot_heights <- c(plot_heights, 3)
    }
  }
  
  # heatmap for mutations
  if(length(mutation_genes) > 0){
    t <- try(res[["mutation_data"]] <- t(DB[["Mutation_gene"]][mutation_genes,ordered_samples,drop=FALSE]))
    if("try-error" %in% class(t)) {
      removeModal()
      not_in_list = mutation_genes[!(mutation_genes %in% rownames(DB[["Mutation_gene"]]))]
      not_in_list = paste(not_in_list, collapse =', ')
      print(sprintf("genes [ %s ] not found in DB.", not_in_list))
      showModal(modalDialog(
        title = "Input genes/features not found.", footer = modalButton("OK"), easyClose = TRUE,
        div(class = "busy",
            p(sprintf("genes [ %s ] are not belong into the database you previously defined. Please remove them from your input.", not_in_list)),
            style = "margin: auto; text-align: center"
        )
      ))
      return()
    }
    colnames(res[["mutation_data"]]) <- paste("mutation_", colnames(res[["mutation_data"]]), sep="")
    mutation_genes = rev(mutation_genes) # This is for descending order of names in y axis
    mutations <- melt(t(DB[["Mutation_gene"]][mutation_genes,ordered_samples,drop=FALSE]))
    mutations[,"value"] <- as.factor(mutations[,"value"])
    mutations[["Var2"]] <- factor(mutations[["Var2"]], levels=mutation_genes, ordered = TRUE)
    mutations[["Var1"]] <- factor(mutations[["Var1"]], levels=ordered_samples, ordered = TRUE)
    PL[["mutation_plot"]] <- ggplot(mutations, aes(Var1, Var2)) + 
      geom_tile(aes(fill = value), colour = "grey50") +
      scale_fill_manual(values = c("white","steelblue")) +
      theme_bw() + my_theme +
      labs(x="", y="")
    plot_heights <- c(plot_heights, length(mutation_genes))
  }
  
  
  # heatmap for protein
  if(length(protein_genes)>0){
    print("--- 5 ---")
    # plot both epithelial expression and stroma expression
    protein_genes = rev(protein_genes) # This is for descending order of names in y axis
    t <- try(protein <- DB[["Protein"]][protein_genes,ordered_samples,drop=FALSE])
    if("try-error" %in% class(t)) {
      removeModal()
      not_in_list = protein_genes[!(protein_genes %in% rownames(DB[["Protein"]]))]
      not_in_list = paste(not_in_list, collapse =', ')
      print(sprintf("genes [ %s ] not found in DB.", not_in_list))
      showModal(modalDialog(
        title = "Input genes/features not found.", footer = modalButton("OK"), easyClose = TRUE,
        div(class = "busy",
            p(sprintf("genes [ %s ] are not belong into the database you previously defined. Please remove them from your input.", not_in_list)),
            style = "margin: auto; text-align: center"
        )
      ))
      return()
    }
    res[["protein_data"]] <- t(protein)
    protein <- t(apply(protein, 1, scale))
    colnames(protein) <- ordered_samples
    protein <- melt(protein)
    protein[["Var2"]] <- factor(protein[["Var2"]], levels=ordered_samples, ordered = TRUE)
    protein[["Var1"]] <- factor(protein[["Var1"]], levels=protein_genes, ordered = TRUE)
    print("levels: (order from bottom to top in plot.")
    print(levels(protein$Var1))
    if(show.protein.name == 1){
      print("--- 6 ---")
      PL[["protein_plot"]] <- ggplot(protein,  aes(Var2, Var1)) + 
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(protein$value),max(protein$value),(max(protein$value)-min(protein$value))/5), digits=1)) +
        guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical") +
        theme_bw() + my_theme +
        labs(x="", y="") + theme(plot.margin=unit(c(0,0,0,0), "cm"))
    }
    else{
      print("--- 7 ---")
      PL[["protein_plot"]] <- ggplot(protein,  aes(Var2, Var1)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "blue4" , mid="white", high = "red", breaks=round(seq(min(protein$value),max(protein$value),(max(protein$value)-min(protein$value))/5), digits=1)) +
        guides(fill = guide_colorbar(barwidth = 10, barheight = 1), direction = "vertical") +
        theme_bw() + my_theme +
        labs(x="", y="") + theme(plot.margin=unit(c(0,0,0,0), "cm")) +
        theme(axis.text.y = element_blank())
    }
    print(length(protein_genes))
    if(length(protein_genes)>50){
      plot_heights <- c(plot_heights, length(protein_genes)*protein_size)
    }
    else{
      plot_heights <- c(plot_heights, length(protein_genes))
    }
  }
  
  t <- try(res[["clin_data"]] <- DB[["Clinical"]][ordered_samples,clinical_lab, drop=FALSE])
  if("try-error" %in% class(t)) {
    removeModal()
    not_in_list = ordered_samples[!(ordered_samples %in% rownames(DB[["Clinical"]]))]
    not_in_list = paste(not_in_list, collapse =', ')
    print(sprintf("Samples [ %s ] not found in DB.", not_in_list))
    showModal(modalDialog(
      title = "Input samples not found.", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p(sprintf("Samples [ %s ] are not belong into the database (Clinical) you previously defined. Please remove them from your input.", not_in_list)),
          style = "margin: auto; text-align: center"
      )
    ))
    return()
  }
  # heatmap for categorical clinical data
  clinical_cat_lab <- intersect(DB[["Clinical_cat_lab"]], clinical_lab)
  clinical_cat <- list()
  if(length(clinical_cat_lab)>0){
    for( i in clinical_cat_lab){
      clinical_cat[[i]] <- melt(t(DB[["Clinical"]][ordered_samples,i, drop=FALSE]))
      clinical_cat[[i]][,"value"] <- as.factor(clinical_cat[[i]][,"value"])
      clinical_cat[[i]][["Var2"]] <- factor(clinical_cat[[i]][["Var2"]], levels=ordered_samples, ordered = TRUE)
      PL[[i]] <- ggplot(clinical_cat[[i]], aes(Var2, Var1)) +
        geom_tile(aes(fill = value), colour = "white") +
        theme_bw() + my_theme +
        labs(x="", y="")
      plot_heights <- c(plot_heights, 1 ) 
    }
  }
  
  # heatmap for quantitative clinical data
  clinical_quan_lab <- intersect(DB[["Clinical_quan_lab"]], clinical_lab)
  if(length(clinical_quan_lab)>0){
    clinical_quan <- melt(t(DB[["Clinical"]][ordered_samples,clinical_quan_lab, drop=FALSE]))
    clinical_quan[["Var2"]] <- factor(clinical_quan[["Var2"]], levels=ordered_samples, ordered = TRUE)
    clinical_quan[,"value"] <- as.numeric(clinical_quan[,"value"])
    PL[["clinical_quan_plot"]] <- ggplot(clinical_quan,  aes(Var2, Var1)) + 
      geom_tile(aes(fill = value)) + 
      scale_fill_gradient(low = "white" , high = "red") + 
      theme_bw() + my_theme +
      labs(x="", y="")
    plot_heights <- c(plot_heights, length(clinical_quan_lab) )
  }
  
  d <- as.data.frame(res)
  res <- list()
  res[["table"]] <- d
  # plot all plots
  res[["plot"]] <- plot_grid(plotlist=PL, align="v", ncol=1, rel_heights=plot_heights)
  # return sample order also
  res[["sample_order_res"]] <- sample_order_res
  
  res
}


########################################################################
# rna based
########################################################################
parse_criteria <- function(rules){
  to_return <- c(0,1)
  for(i in rules){
    split1 <- gsub(" ", "", unlist(strsplit(i, ">", fixed = TRUE)))
    split2 <- gsub(" ", "", unlist(strsplit(i, "<", fixed = TRUE)))
    if(length(split1) > 1){ to_return[1] <- max(to_return[1], as.numeric(split1[2]))}
    if(length(split2) > 1){ to_return[2] <- min(to_return[2], as.numeric(split2[2]))}
  }
  to_return
}
my_heatmap_rna <- function(mode, clust_para=list(method="hc"), mutation_genes, image_features = vector(),
                           protein_genes, clinical_lab, rna_criteria, rna_genes, show.RNA.name = 1,
                           show.protein.name = 1, sort.mutation=F, sort.rna=F,
                           sort.image=F, sort.protein=F){
  rna_samples <- colnames(DB[["RNA"]])[!apply(DB[["RNA"]], 2, function(x){all(is.na(x))})]
  
  # order by user specified genes
  if(mode==0){
    my_heatmap_mutation(mutation_genes=mutation_genes, image_features=image_features, rna_genes=rna_genes, 
                        protein_genes=protein_genes, clinical_lab=clinical_lab, 
                        order_by="rna", selected_samples = rna_samples,
                        clust_para=clust_para, show.RNA.name = show.RNA.name, show.protein.name = show.protein.name,
                        sort.mutation, sort.rna, sort.image, sort.protein)
  } else {
    rna_genes <- rownames(DB[["RNA"]])
    maxExp_filters <- rna_criteria[grep("maxExp", rna_criteria)]
    if(length(maxExp_filters)>0){
      temp_res <- gene_filter(DB[["RNA"]][rna_genes,rna_samples], "maxExp", parse_criteria(maxExp_filters))
      rna_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    var_filters <- rna_criteria[grep("var", rna_criteria)]
    if(length(var_filters)>0){
      temp_res <- gene_filter(DB[["RNA"]][rna_genes,rna_samples], "var", parse_criteria(var_filters))
      rna_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    cv_filters <- rna_criteria[grep("cv", rna_criteria)]
    if(length(cv_filters)>0){
      temp_res <- gene_filter(DB[["RNA"]][rna_genes,rna_samples], "cv", parse_criteria(cv_filters))
      rna_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    my_heatmap_mutation(mutation_genes=mutation_genes, image_features=image_features, rna_genes=rna_genes, 
                        protein_genes=protein_genes, clinical_lab=clinical_lab, 
                        order_by="rna", rna_size=0.01, selected_samples = rna_samples,
                        clust_para=clust_para, show.RNA.name = show.RNA.name, show.protein.name = show.protein.name,
                        sort.mutation, sort.rna, sort.image, sort.protein)
  }
}

divide_patients_by_cuttree <- function(hclust_tree){
  d <- list()
  for(i in 2:10){
    d[[paste("G", i, sep="")]] <- cutree(hclust_tree, k=i)
  }
  d <- as.data.frame(d)
  d
}

########################################################################
# protein based
########################################################################
my_heatmap_protein <- function(mode, mutation_genes, image_features = vector(), rna_genes,
                               clinical_lab, protein_criteria, protein_genes, show.RNA.name=1,
                               show.protein.name = 1, sort.mutation=F, sort.rna=F,
                               sort.image=F, sort.protein=F){
  protein_samples <- colnames(DB[["Protein"]])[!apply(DB[["Protein"]], 2, function(x){any(is.na(x))})]
  print(length(protein_samples))
  # order by user specified genes
  if(mode==0){
    my_heatmap_mutation(mutation_genes=mutation_genes, image_features=image_features, protein_genes=protein_genes, 
                        rna_genes=rna_genes, clinical_lab=clinical_lab, protein_size=0.01,
                        order_by="protein", selected_samples = protein_samples, show.RNA.name = show.RNA.name, show.protein.name = show.protein.name,
                        sort.mutation, sort.rna, sort.image, sort.protein)
  } else {
    protein_genes <- rownames(DB[["Protein"]])
    maxExp_filters <- protein_criteria[grep("maxExp", protein_criteria)]
    if(length(maxExp_filters)>0){
      temp_res <- gene_filter(DB[["Protein"]][protein_genes,protein_samples], "maxExp", parse_criteria(maxExp_filters))
      protein_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    var_filters <- protein_criteria[grep("var", protein_criteria)]
    if(length(var_filters)>0){
      temp_res <- gene_filter(DB[["Protein"]][protein_genes,protein_samples], "var", parse_criteria(var_filters))
      protein_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    cv_filters <- protein_criteria[grep("cv", protein_criteria)]
    if(length(cv_filters)>0){
      temp_res <- gene_filter(DB[["Protein"]][protein_genes,protein_samples], "cv", parse_criteria(cv_filters))
      protein_genes <- temp_res[["genes"]]
      print(temp_res[["cutoff"]])
    }
    my_heatmap_mutation(mutation_genes=mutation_genes, image_features=image_features, rna_genes=rna_genes, 
                        protein_genes=protein_genes, clinical_lab=clinical_lab, protein_size=0.05,
                        order_by="protein", selected_samples = protein_samples, show.RNA.name = show.RNA.name, show.protein.name = show.protein.name,
                        sort.mutation, sort.rna, sort.image, sort.protein)
  }
}
