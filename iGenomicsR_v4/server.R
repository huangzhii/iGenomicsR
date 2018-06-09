library(data.table)
library(DT)
library(RColorBrewer)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(gplots)
library(plyr)
library(GGally)
library(survival)
source("utils.R")
source("my_heatmap.R")
source("my_analysis.R")
source("gene_filter.R")

options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB
options(shiny.sanitize.errors = FALSE)
options(stringsAsFactors = FALSE)


function(input, output, session) {
  output$check1 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check2 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check3 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check4 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check5 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  
  ########################################################################
  # data upload panel 
  ########################################################################
  observeEvent(input$action1,{
    hideTab(inputId = "data.navigator", target = "Mutation")
    hideTab(inputId = "data.navigator", target = "Image features")
    hideTab(inputId = "data.navigator", target = "RNA expression")
    hideTab(inputId = "data.navigator", target = "Protein expression")
    hideTab(inputId = "data.integration", target = "Mutation")
    hideTab(inputId = "data.integration", target = "Image features")
    hideTab(inputId = "data.integration", target = "RNA expression")
    hideTab(inputId = "data.integration", target = "Protein expression")
    
    analysisDataTypeList = list()
    if("Mutation_gene" %in% names(DB)){
      output$hasMutationData <-  reactive(TRUE)
      outputOptions(output, "hasMutationData", suspendWhenHidden = FALSE)
      analysisDataTypeList["Mutation"] = 0
      showTab(inputId = "data.navigator", target = "Mutation")
      showTab(inputId = "data.integration", target = "Mutation")
    }
    if("Image" %in% names(DB)){
      output$hasImageData <-  reactive(TRUE)
      outputOptions(output, "hasImageData", suspendWhenHidden = FALSE)
      analysisDataTypeList["Image feature"] = 4
      showTab(inputId = "data.navigator", target = "Image features")
      showTab(inputId = "data.integration", target = "Image features")
    }
    if("RNA" %in% names(DB)){
      output$hasRNAData <-  reactive(TRUE)
      outputOptions(output, "hasRNAData", suspendWhenHidden = FALSE)
      analysisDataTypeList["RNA expression"] = 1
      showTab(inputId = "data.navigator", target = "RNA expression")
      showTab(inputId = "data.integration", target = "RNA expression")
    }
    if("Protein" %in% names(DB)){
      output$hasProteinData <-  reactive(TRUE)
      outputOptions(output, "hasProteinData", suspendWhenHidden = FALSE)
      analysisDataTypeList["Protein expression"] = 2
      showTab(inputId = "data.navigator", target = "Protein expression")
      showTab(inputId = "data.integration", target = "Protein expression")
    }
    
    analysisDataTypeList["Clinical data"] = 3
    analysisDataTypeList["Survival"] = 5
    
    output$AnalysisDataTypeUI<-renderUI({
      awesomeRadio("AnalysisDataType", "", analysisDataTypeList)
    })
    
    if(length(DB) >= 3 & !is.null(DB[["Clinical"]])){
      output$data_ready_flag <-  reactive(TRUE)
      outputOptions(output, "data_ready_flag", suspendWhenHidden = FALSE)
      sendSweetAlert(session, title = "File Upload Success", text = NULL, type = "success",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      output$dataUploadSummary <- renderTable({
        return(summarize_dataUpload())
      }, include.colnames=FALSE)
      
      
      output$sampleVenn <- renderPlot({
        # Venn Plot
        vennData <- list()
        for(i in setdiff(names(DB), "Clinical")){
          vennData[[i]] <- colnames(DB[[i]])
        }
        if("Clinical" %in% names(DB)){
          vennData[["Clinical"]] <- rownames(DB[["Clinical"]]) 
        }
        venn(vennData)
      }, height = 400, width = 600)
      session$sendCustomMessage("buttonCallbackHandler", "tab1")
    }
    else{
      sendSweetAlert(session, title = "Insufficient Input Data", text = "Please upload required clinical file and at least two omics data.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    }
  })
  
  
  loadData.mutation <- function(){
    if (is.null(input$csvfile_mutation)){
      return(NULL)}else{
        DB[["Mutation_gene"]] <<- read.table(input$csvfile_mutation$datapath, sep=input$sep, header=TRUE, row.names = 1)
        colnames(DB[["Mutation_gene"]]) <<- gsub(".", "-", colnames(DB[["Mutation_gene"]]), fixed = TRUE)
        DB[["Mutation_gene"]][is.na(DB[["Mutation_gene"]])] <<- 0
        session$sendCustomMessage("buttonCallbackHandler", "tab1")
      }
    return(dim( DB[["Mutation_gene"]] )) 
  }
  loadData.mRNA <- function(){
    if (is.null(input$csvfile_mRNA)){
      return(NULL)}else{
        DB[["RNA"]] <<- read.table(input$csvfile_mRNA$datapath, sep=input$sep, header=TRUE, row.names = 1)
        colnames(DB[["RNA"]]) <<- gsub(".", "-", colnames(DB[["RNA"]]), fixed = TRUE)
        DB[["RNA"]] <<- apply(DB[["RNA"]], c(1,2), as.numeric)
        DB[["RNA"]][is.na(DB[["RNA"]])] <<- 0
        session$sendCustomMessage("buttonCallbackHandler", "tab1")
      }
    return(dim( DB[["RNA"]] ))
  }
  loadData.protein <- function(){
    if (is.null(input$csvfile_protein)){
      return(NULL)}else{
        DB[["Protein"]] <<- read.table(input$csvfile_protein$datapath, sep=input$sep, header=TRUE, row.names = 1)
        colnames(DB[["Protein"]]) <<- gsub(".", "-", colnames(DB[["Protein"]]), fixed = TRUE)
        DB[["Protein"]] <<- apply(DB[["Protein"]], c(1,2), as.numeric)
        DB[["Protein"]][is.na(DB[["Protein"]])] <<- 0
        session$sendCustomMessage("buttonCallbackHandler", "tab1")
      }
    return(dim(DB[["Protein"]]))
  }
  loadData.clinical<- function(){
    if (is.null(input$csvfile_clinical)){
      return(NULL)}else{
        DB[["Clinical"]] <<- read.table(input$csvfile_clinical$datapath, sep=input$sep, header=TRUE, row.names = 1)
        session$sendCustomMessage("buttonCallbackHandler", "tab1")
      }
    return(dim(DB[["Clinical"]]))
  }
  loadData.image<- function(){
    if (is.null(input$csvfile_image)){
      return(NULL)}else{
        DB[["Image"]] <<- read.table(input$csvfile_image$datapath, sep=input$sep, header=TRUE, row.names = 1)
        DB[["Image"]] <<- t(DB[["Image"]])
        # Data Integration
        DB[["Image"]][is.na(DB[["Image"]])] <<- 0
        output$ImageInputFeaturesUI <- renderUI({
          selectInput(inputId="ImageInputFeatures", label="Select features here",
                      choices = rownames(DB[["Image"]]),
                      multiple = T)
        })
        
        output$ImageInputFeaturesSubUI_Mutation <- renderUI({
          selectInput(inputId="ImageInputFeaturesSelectionForMutation", label="Select features here",
                      choices = rownames(DB[["Image"]]),
                      multiple = T)
        })
        output$ImageInputFeaturesSubUI_RNA <- renderUI({
          selectInput(inputId="ImageInputFeaturesSelectionForRNA", label="Select features here",
                      choices = rownames(DB[["Image"]]),
                      multiple = T)
        })
        output$ImageInputFeaturesSubUI_Protein <- renderUI({
          selectInput(inputId="ImageInputFeaturesSelectionForProtein", label="Select features here",
                      choices = rownames(DB[["Image"]]),
                      multiple = T)
        })
        output$ImageInputFeaturesSubUI_Clinical <- renderUI({
          selectInput(inputId="ImageInputFeaturesSelectionForClinical", label="Select features here",
                      choices = rownames(DB[["Image"]]),
                      multiple = T)
        })
        session$sendCustomMessage("buttonCallbackHandler", "tab1")
      }
    return(dim(DB[["Image"]]))
  }
  
  observeEvent(loadData.mutation(),{output$check1 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.mRNA(),{output$check2 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.protein(),{output$check3 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.clinical(),{output$check4 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.image(),{output$check5 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  
  summarize_dataUpload <- eventReactive(input$uploadSummaryButton, {
    smartModal(error=F, title = "Summarizing Uploaded Data", content = "Summarizing uploaded data, please wait for a little while...")
    d <- list()
    d[["dataTypes"]] <- c("Uploaded data types", paste(names(DB), collapse = ";"))
    for(i in setdiff(names(DB), "Clinical")){
      d[[i]] <- c(paste("# of samples for", i), ncol(DB[[i]]))
    }
    if( "Clinical" %in% names(DB)){
      d[["clinFeatures"]] <- c("Clinical features", paste(names(DB[["Clinical"]]), collapse = ";"))
      for(i in names(DB[["Clinical"]])){
        d[[i]] <- c(paste("# of NA in", i), sum(is.na(toupper(DB[["Clinical"]][,i])) | DB[["Clinical"]][,i]==""))
      }
    }
    removeModal()
    return(t(as.data.frame(d)))
  })
  
  sampleSelection <- observeEvent(input$sampleSelectionButton, {
    allDataTypes <- colnames(DB) # [1] "Mutation_gene" "RNA" "Protein" "Clinical"
    if ("Clinical" %in% allDataTypes){
      samples <- rownames(DB[["Clinical"]])
    } 
    else if (length(setdiff(allDataTypes, "Clinical")) >= 1 ){
      samples <- colnames(DB[[setdiff(allDataTypes, "Clinical")[1]]])
      for(i in setdiff(allDataTypes, "Clinical")){
        samples <- intersect(samples, colnames(DB[[i]]))
      }
    }
    for(i in setdiff(allDataTypes, "Clinical")){
      DB[[i]] <<- DB[[i]][,samples]
    }
    if ("Clinical" %in% allDataTypes){
      DB[["Clinical"]] <<- DB[[i]][samples,]
    } 
    DB[["Clinical_cat_lab"]] <<- CatClin()
    DB[["Clinical_quan_lab"]] <<- QuanClin()
    assign("DB", DB, envir = .GlobalEnv)
    # print(names(DB))
    # 
    # print(head(DB[["Clinical"]]))
    # print(head(DB[["Clinical_cat_lab"]]))
    # print(head(DB[["Clinical_quan_lab"]]))
    
    sendSweetAlert(session, title = "Success", text = "Now only keep samples shared by all subjects", type = "success",
                   btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    if ("Mutation_gene" %in% names(DB)){session$sendCustomMessage("download_seleted_data_ready_mutation", "lalala")}
    if ("RNA" %in% names(DB)){session$sendCustomMessage("download_seleted_data_ready_rna", "lalala")}
    if ("Image" %in% names(DB)){session$sendCustomMessage("download_seleted_data_ready_image", "lalala")}
    if ("Protein" %in% names(DB)){session$sendCustomMessage("download_seleted_data_ready_protein", "lalala")}
    if ("Clinical" %in% names(DB)){session$sendCustomMessage("download_seleted_data_ready_clinical", "lalala")}
  })
  
  get_all_clin <- reactive({
    ClinList <- list()
    for(i in colnames(DB[["Clinical"]])){
      ClinList[[i]] <- i
    }
    ClinList
  })
  output$checkCatClinUI <- renderUI({
    if (is.null(input$csvfile_clinical))
      return()
    checkboxGroupInput('checkCatClin', 'Select categorical clinical feature', get_all_clin(), 
                       selected = "")
  })
  output$checkQuanClinUI <- renderUI({
    if (is.null(input$csvfile_clinical))
      return()
    checkboxGroupInput('checkQuanClin', 'Select quantitative clinical feature', get_all_clin(), 
                       selected = "")
  })
  CatClin <- reactive({input$checkCatClin})
  QuanClin <- reactive({input$checkQuanClin})
  
  output$downloadSelectedMutation <- downloadHandler(
    filename = function() { "SelectedMutationTable.csv" },
    content = function(file) {
      sampleSelection
      write.csv(DB[["Mutation_gene"]], file, row.names=TRUE)
    }) ###
  output$downloadSelectedImageFeature <- downloadHandler(
    filename = function() { "SelectedImageFeatureTable.csv" },
    content = function(file) {
      sampleSelection
      write.csv(DB[["Image"]], file, row.names=TRUE)
    }) ###
  output$downloadSelectedRNA <- downloadHandler(
    filename = function() { "SelectedRNATable.csv" },
    content = function(file) {
      sampleSelection
      write.csv(DB[["RNA"]], file, row.names=TRUE)
    }) ###
  output$downloadSelectedProtein <- downloadHandler(
    filename = function() { "SelectedProteinTable.csv" },
    content = function(file) {
      sampleSelection
      write.csv(DB[["Protein"]], file, row.names=TRUE)
    }) ###
  output$downloadSelectedClin <- downloadHandler(
    filename = function() { "SelectedClinicalTable.csv" },
    content = function(file) {
      sampleSelection
      write.csv(DB[["Clinical"]], file, row.names=TRUE)
    }) ###
  
  observeEvent(input$action2,{
    if(length(DB) >= 3 & !is.null(DB[["Clinical"]])){
      session$sendCustomMessage("buttonCallbackHandler", "tab2")
    }
    else{
      sendSweetAlert(session, title = "Insufficient Input Data", text = "Please upload required clinical file and at least two omics data.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    }
  })
  
  
  ########################################################################
  # Data navigator panel 
  ########################################################################
  
  observeEvent(input$action.navigator.mutation,{
    
    # output gene mutation test result
    output$geneMutationTestResTable <- DT::renderDataTable({
      # print(head(DB[["Mutation_gene"]]))
      DB[["mutation_gene_test"]] <<- run_gene_mutation_association(DB[["Mutation_gene"]])
      # print(head(DB[["mutation_gene_test"]]))
      d <- DB[["mutation_gene_test"]]
      d[,"oddsRatio"] <- format(as.numeric(d[,"oddsRatio"]),nsmall=2, digits=2)
      d[,"pvalue"] <- format(as.numeric(d[,"pvalue"]),scientific=TRUE, nsmall=2,digits=2)
      d[,"adj_pvalue"] <- format(as.numeric(d[,"adj_pvalue"]),scientific=TRUE, nsmall=2,digits=2)
      return(d)
    },selection="none",options=list(searching=F, ordering=F))#,extensions = 'Responsive'
    
    output$dowloadGeneMutationTestRes <- downloadHandler(
      filename = function() { "mutation_gene_association_test.csv" },
      content = function(file) {
        write.csv(DB[["mutation_gene_test"]], file, row.names=FALSE)
      })
    
    # output selected mutations
    output$selectedGeneMutationsTable <- DT::renderDataTable({
      DB[["Mutation_gene"]][unlist(strsplit(gsub(" ", "", input$MutationInputGenes), ",", fixed = TRUE)),]
    },selection="none",options=list(searching=F, ordering=F))#,extensions = 'Responsive'
  })
  
  observeEvent(input$action.navigator.RNA,{
    output$RNADotPlot <- renderPlot({
      genes <- c(input$navigator.RNA.expression.gene.1, input$navigator.RNA.expression.gene.2)
      d <- DB[["RNA"]][genes,,drop=FALSE]
      d <- d[,!apply(d, 2, function(z){any(is.na(z))})]
      x <- d[genes[1],,drop=TRUE]
      y <- d[genes[2],,drop=TRUE]
      corr <- cor(x, y)
      plot(x, y, xlab=genes[1], ylab=genes[2], main="RNA expression", sub=paste("correlation:", corr))
    }, height = 500, width = 500)
  })
  
  observeEvent(input$action.navigator.protein,{
    # save(DB, file = "~/Desktop/DB.Rdata")
    output$ProteinDotPlot1 <- renderPlot({
      genes <- c(input$navigator.protein.expression.gene.1, input$navigator.protein.expression.gene.2)
      d <- DB[["Protein"]][genes,,drop=FALSE]
      d <- d[,!apply(d, 2, function(z){any(is.na(z))})]
      d <- apply(d, c(1,2), as.numeric)
      x <- d[genes[1],,drop=TRUE]
      y <- d[genes[2],,drop=TRUE]
      corr <- cor(x, y)
      plot(x, y, xlab=genes[1], ylab=genes[2], main="Protein expression", sub=paste("correlation:", corr))
    }, height = 500, width = 500)
    
  })
  output$ImageFeaturesNavigatorPlot <- renderPlot({
    # dev.off()
    heatmap(DB[["Image"]], margins = c(8, 15) )
  })
  output$ImageFeaturesTable <- DT::renderDataTable({
    DT::datatable(DB[["Image"]], escape=F, selection = 'none', rownames = T,
                  options=list(searching=F, ordering=F)) #, extensions = 'Responsive'
    
  })
  output$ClinicalInfoTable <- DT::renderDataTable({
    DT::datatable(DB[["Clinical"]], escape=F, selection = 'none', rownames = T,
                  options=list(searching=F, ordering=F)) #, extensions = 'Responsive'
    
  })
  
  ########################################################################
  # mutation panel 
  ########################################################################
  
  
  # Oncoplot
  selectedClinFeature <- reactive({
    SCF <- list()
    for(i in c(CatClin(), QuanClin())){
      SCF[[i]] <- i
    }
    return(SCF)
  })
  
  output$OncoPlotClinUI <- renderUI({
    if (!input$OncoPlotHasClin)
      return()
    checkboxGroupInput('OncoPlotClin', '', selectedClinFeature(), selected = "") 
  })
  output$ImageheatClinUI <- renderUI({
    if (!input$ImageheatHasClin)
      return()
    checkboxGroupInput('ImageheatClin', '', selectedClinFeature(), selected = "")
  })
  output$RNAheatClinUI <- renderUI({
    if (!input$RNAheatHasClin)
      return()
    checkboxGroupInput('RNAheatClin', '', selectedClinFeature(), selected = "")
  })
  output$ProteinheatClinUI <- renderUI({
    if (!input$ProteinheatHasClin)
      return()
    checkboxGroupInput('ProteinheatClin', '', selectedClinFeature(), selected = "")
  })
  output$ClinheatClinUI <- renderUI({
    checkboxGroupInput('ClinheatClin', '', c(CatClin(), QuanClin()), 
                       selected = c(CatClin(), QuanClin()))
  })
  output$ClinheatSelectOrderFeatureUI <- renderUI({
    selectInput('ClinheatSelectOrderFeature', "Order samples by", c(CatClin(), QuanClin()), 
                selected = c(CatClin(), QuanClin())[1])
  })
  
  observeEvent(input$action.integration.mutation.inputgenes,{
    OncoPlot_res <<- my_heatmap_mutation(mutation_genes = gsub("\\s","", strsplit(input$genesToPullMutation,",")[[1]]),
                                         image_features = if(input$OncoPlotHasImage){
                                           gsub("\\s","", input$ImageInputFeaturesSelectionForMutation)
                                         },
                                         rna_genes = if(input$OncoPlotHasRna){
                                           gsub("\\s","", strsplit(input$MutationInputRna,",")[[1]])
                                         },
                                         protein_genes = if(input$OncoPlotHasProtein){
                                           gsub("\\s","", strsplit(input$MutationInputProteins,",")[[1]])
                                         },
                                         clinical_lab = input$OncoPlotClin,
                                         order_by="mutation",
                                         sort.mutation = input$do_hclust_mutation,
                                         sort.rna = input$do_hclust_mutation_rna,
                                         sort.image = input$do_hclust_mutation_image,
                                         sort.protein = input$do_hclust_mutation_protein)
    # save(OncoPlot_res, file= "~/Desktop/oncoplot.Rdata")
    height_of_plot <- length(gsub("\\s","", strsplit(input$genesToPullMutation,",")[[1]])) +
      length(input$OncoPlotClin)
    if(input$OncoPlotHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForMutation))
    }
    if(input$OncoPlotHasRna){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$MutationInputRna,",")[[1]]))
    }
    if(input$OncoPlotHasProtein){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$MutationInputProteins,",")[[1]]))
    }
    output$OncoPlot <- renderPlot({
      return(OncoPlot_res[["plot"]])
    }, width=input$myWidth1, height=input$myHeight1/20*height_of_plot)
    
  })
  
  # download ordered data for heatmap
  output$downloadOncoPlotData <- downloadHandler(
    filename = function() { "data_ordered_by_mutation.csv" },
    content = function(file) {
      write.csv( t(OncoPlot_res[["table"]]), file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # image panel 
  ########################################################################
  # Imageheat
  
  observeEvent(input$action.integration.image.inputgenes,{
    Imageheat_res <<- my_heatmap_mutation(image_features = gsub("\\s","", input$ImageInputFeatures),
                                          mutation_genes = if(input$ImageheatHasMutation){
                                            gsub("\\s","", strsplit(input$ImageInputMutations,",")[[1]])
                                          },
                                          rna_genes = if(input$ImageheatHasRna){
                                            gsub("\\s","", strsplit(input$ImageInputRna,",")[[1]])
                                          },
                                          protein_genes = if(input$ImageheatHasProtein){
                                            gsub("\\s","", strsplit(input$ImageInputProteins,",")[[1]])
                                          },
                                          clinical_lab = input$ImageheatClin,
                                          order_by="image",
                                          sort.mutation = input$do_hclust_image_mutation,
                                          sort.rna = input$do_hclust_image_rna,
                                          sort.image = input$do_hclust_image,
                                          sort.protein = input$do_hclust_image_protein)
    # save(Imageheat_res, file= "~/Desktop/oncoplot.Rdata")
    height_of_plot <- length(gsub("\\s","", input$ImageInputFeatures)) +
      length(input$ImageheatClin)
    if(input$ImageheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ImageInputMutations,",")[[1]]))
    }
    if(input$ImageheatHasRna){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ImageInputRna,",")[[1]]))
    }
    if(input$ImageheatHasProtein){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ImageInputProteins,",")[[1]]))
    }
    output$Imageheat <- renderPlot({
      return(Imageheat_res[["plot"]])
    }, width=input$myWidth2, height=input$myHeight2/20*height_of_plot)
    
  })
  
  # download ordered data for heatmap
  output$downloadImageheatData <- downloadHandler(
    filename = function() { "data_ordered_by_image_feature.csv" },
    content = function(file) {
      write.csv( t(Imageheat_res[["table"]]), file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # rna expression panel 
  ########################################################################
  # gene expression clustering and heatmap
  observeEvent(input$action.integration.RNA.denovo,{
    clust_para <- list()
    clust_para[["method"]] <- "hc"
    rna_RNAheatClustPara <- clust_para
    RNAheat_res <<- my_heatmap_rna(mode = 1, #denovo
                                   clust_para = rna_RNAheatClustPara,
                                   mutation_genes = if(input$RNAheatHasMutation){
                                     gsub("\\s","", strsplit(input$RNAheatInputMutation,",")[[1]])
                                   },
                                   image_features = if(input$RNAheatHasImage){
                                     gsub("\\s","", input$ImageInputFeaturesSelectionForRNA)
                                   },
                                   protein_genes = if(input$RNAheatHasProtein){
                                     gsub("\\s","", strsplit(input$RNAheatInputProteins,",")[[1]])
                                   },
                                   clinical_lab=input$RNAheatClin,
                                   rna_criteria = strsplit(input$RNAheatGeneCutoff,"and")[[1]],
                                   rna_genes = NULL,
                                   show.RNA.name = input$show.RNA.name,
                                   sort.mutation = input$do_hclust_rna_mutation,
                                   sort.rna = F,
                                   sort.image = input$do_hclust_rna_image,
                                   sort.protein = input$do_hclust_rna_protein)
    
    height_of_plot <- 40 + length(input$RNAheatClin)
    if(input$RNAheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$RNAheatInputMutation,",")[[1]]))
    }
    if(input$RNAheatHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForRNA))
    }
    if(input$RNAheatHasProtein){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$RNAheatInputProteins,",")[[1]]))
    }
    output$RNAheat <- renderPlot({
      return(RNAheat_res[["plot"]])
    }, height = input$myHeight3/40*height_of_plot, width = input$myWidth3)
    
    # gene clustering dendrogram
    output$RNAdendro <- renderPlot({
      if((clust_para[["method"]] == "hc")){
        plot(RNAheat_res[["sample_order_res"]][["hc"]], cex=0.5)
      }
    }, height = input$myHeight3/2, width = input$myWidth3)
  })
  
  observeEvent(input$action.integration.RNA.inputgenes,{
    clust_para <- list()
    clust_para[["method"]] <- c("hc", "km")[as.numeric(input$RNAheatClustMethod) + 1]
    if(clust_para[["method"]] == "km"){
      clust_para[["k"]] <- as.numeric(input$RNAheatKmeansK)
    }
    rna_RNAheatClustPara <- clust_para
    RNAheat_res <<- my_heatmap_rna(mode = 0, # gene lists
                                   clust_para = rna_RNAheatClustPara,
                                   mutation_genes = if(input$RNAheatHasMutation){
                                     gsub("\\s","", strsplit(input$RNAheatInputMutation,",")[[1]])
                                   },
                                   image_features = if(input$RNAheatHasImage){
                                     gsub("\\s","", input$ImageInputFeaturesSelectionForRNA)
                                   },
                                   protein_genes = if(input$RNAheatHasProtein){
                                     gsub("\\s","", strsplit(input$RNAheatInputProteins,",")[[1]])
                                   }, 
                                   clinical_lab=input$RNAheatClin,
                                   rna_criteria = NULL,
                                   rna_genes = gsub("\\s","", strsplit(input$RNAheatInputGenes,",")[[1]],
                                    sort.mutation = input$do_hclust_rna_mutation,
                                    sort.rna = F,
                                    sort.image = input$do_hclust_rna_image,
                                    sort.protein = input$do_hclust_rna_protein)
    )
    
    height_of_plot <- length(gsub("\\s","", strsplit(input$RNAheatInputGenes,",")[[1]])) + length(input$RNAheatClin)
    if(input$RNAheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$RNAheatInputMutation,",")[[1]]))
    }
    if(input$RNAheatHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForRNA))
    }
    if(input$RNAheatHasProtein){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$RNAheatInputProteins,",")[[1]]))
    }
    output$RNAheat <- renderPlot({
      return(RNAheat_res[["plot"]])
    }, height = input$myHeight3/20*height_of_plot, width = input$myWidth3)
    
    # gene clustering dendrogram
    if((clust_para[["method"]] == "hc")){
      output$RNAdendro <- renderPlot({
        plot(RNAheat_res[["sample_order_res"]][["hc"]], cex=0.5)
      }, height = input$myHeight3/2, width = input$myWidth3)
    }
    if((clust_para[["method"]] == "km")){
      output$RNAdendro <- renderPlot({
      }, height = 0, width = 0)
    }
  })
  
  # download ordered data for heatmap
  output$downloadRNAheatData <- downloadHandler(
    filename = function() { "data_ordered_by_rna.csv" },
    content = function(file) {
      write.csv(t(RNAheat_res[["table"]]), 
                file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # protein panel 
  ########################################################################
  # plot heatmap order by protein
  observeEvent(input$action.integration.protein.denovo,{
    Proteinheat_res <<- my_heatmap_protein(mode=1,
                                           mutation_genes=if(input$ProteinheatHasMutation){
                                             gsub("\\s","", strsplit(input$ProteinheatInputMutation,",")[[1]])
                                           },
                                           image_features = if(input$ProteinheatHasImage){
                                             gsub("\\s","", input$ImageInputFeaturesSelectionForProtein)
                                           },
                                           rna_genes=if(input$ProteinheatHasRNA){
                                             gsub("\\s","", strsplit(input$ProteinheatInputRNA,",")[[1]])
                                           }, 
                                           clinical_lab=input$ProteinheatClin,
                                           protein_criteria = strsplit(input$ProteinheatGeneCutoff,"and")[[1]], 
                                           protein_genes = NULL,
                                           show.protein.name = input$show.protein.name,
                                           sort.mutation = input$do_hclust_protein_mutation,
                                           sort.rna = input$do_hclust_protein_rna,
                                           sort.image = input$do_hclust_protein_image,
                                           sort.protein = F)
    
    height_of_plot <- 40 + length(input$ProteinheatClin)
    if(input$ProteinheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ProteinheatInputMutation,",")[[1]]))
    }
    if(input$ProteinheatHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForProtein))
    }
    if(input$ProteinheatHasRNA){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ProteinheatInputRNA,",")[[1]]))
    }
    output$Proteinheat <- renderPlot({
      Proteinheat_res[["plot"]]
    }, height = input$myHeight4/40*height_of_plot, width = input$myWidth4)
    # gene clustering dendrogram
    output$Proteindendro <- renderPlot({
      plot(Proteinheat_res[["sample_order_res"]][["hc"]], cex=0.5)
    }, height = input$myHeight4/2, width = input$myWidth4)
  })
  
  observeEvent(input$action.integration.protein.inputgenes,{
    Proteinheat_res <<- my_heatmap_protein(mode=0,
                                           mutation_genes=if(input$ProteinheatHasMutation){
                                             gsub("\\s","", strsplit(input$ProteinheatInputMutation,",")[[1]])
                                           },
                                           image_features = if(input$ProteinheatHasImage){
                                             gsub("\\s","", input$ImageInputFeaturesSelectionForProtein)
                                           },
                                           rna_genes=if(input$ProteinheatHasRNA){
                                             gsub("\\s","", strsplit(input$ProteinheatInputRNA,",")[[1]])
                                           }, 
                                           clinical_lab=input$ProteinheatClin,
                                           protein_criteria = NULL,
                                           protein_genes = gsub("\\s","", strsplit(input$ProteinheatInputGenes,",")[[1]]),
                                           sort.mutation = input$do_hclust_protein_mutation,
                                           sort.rna = input$do_hclust_protein_rna,
                                           sort.image = input$do_hclust_protein_image,
                                           sort.protein = F)
    
    height_of_plot <- length(gsub("\\s","", strsplit(input$ProteinheatInputGenes,",")[[1]])) + length(input$ProteinheatClin)
    if(input$ProteinheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ProteinheatInputMutation,",")[[1]]))
    }
    if(input$ProteinheatHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForProtein))
    }
    if(input$ProteinheatHasRNA){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ProteinheatInputRNA,",")[[1]]))
    }
    output$Proteinheat <- renderPlot({
      Proteinheat_res[["plot"]]
    }, height = input$myHeight4/20*height_of_plot, width = input$myWidth4)
    # gene clustering dendrogram
    output$Proteindendro <- renderPlot({
      plot(Proteinheat_res[["sample_order_res"]][["hc"]], cex=0.5)
    }, height = input$myHeight4/2, width = input$myWidth4)
    
  })
  
  # download ordered data for heatmap
  output$downloadProteinheatData <- downloadHandler(
    filename = function() { "data_ordered_by_protein.csv" },
    content = function(file) {
      write.csv(t(Proteinheat_res[["table"]]), 
                file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # Clinical panel 
  ########################################################################
  
  # plot heatmap order by clinical feature
  observeEvent(input$action.integration.clinical,{
    Clinheat_res <<- my_heatmap_mutation(mutation_genes = if(input$ClinheatHasMutation){
      gsub("\\s","", strsplit(input$ClinheatInputMutation,",")[[1]])
    },
    image_features = if(input$ClinheatHasImage){
      gsub("\\s","", input$ImageInputFeaturesSelectionForClinical)
    },
    rna_genes = if(input$ClinheatHasRNA){
      gsub("\\s","", strsplit(input$ClinheatInputRNA,",")[[1]])
    },
    protein_genes = if(input$ClinheatHasProtein){
      gsub("\\s","", strsplit(input$ClinheatInputProtein,",")[[1]])
    },
    clinical_lab = input$ClinheatClin,
    order_by = "clinical",
    order_clin_feature = input$ClinheatSelectOrderFeature,
    sort.mutation = input$do_hclust_clin_mutation,
    sort.rna = input$do_hclust_clin_rna,
    sort.image = input$do_hclust_clin_image,
    sort.protein = input$do_hclust_clin_protein)
    
    height_of_plot <- length(input$ClinheatClin)
    if(input$ClinheatHasMutation){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ClinheatInputMutation,",")[[1]]))
    }
    if(input$ClinheatHasImage){
      height_of_plot = height_of_plot + length(gsub("\\s","", input$ImageInputFeaturesSelectionForClinical))
    }
    if(input$ClinheatHasRNA){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ClinheatInputRNA,",")[[1]]))
    }
    if(input$ClinheatHasProtein){
      height_of_plot = height_of_plot + length(gsub("\\s","", strsplit(input$ClinheatInputProtein,",")[[1]]))
    }
    
    output$Clinheat <- renderPlot({
      Clinheat_res[["plot"]]
    }, height = input$myHeight5/20*height_of_plot, width = input$myWidth5)
  })
  # download ordered data for heatmap
  output$downloadClinheatData <- downloadHandler(
    filename = function() { "data_ordered_by_clinical.csv" },
    content = function(file) {
      write.csv( t(Clinheat_res[["table"]]), file, row.names=TRUE)
    }) ###
  
  
  
  ########################################################################
  # analysis panel
  ########################################################################
  get_patient_groups <- reactive({
    if(is.null(input$patientGroups1) | is.null(input$patientGroups2)){
      sendSweetAlert(session, title = "Error", text = "Insufficient Input Data.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    else if(input$patientGroups1 == "" | input$patientGroups2 == "") {
      get_patient_groups = NULL
    }
    tmp1<-matrix(strsplit(input$patientGroups1, "\n")[[1]])
    tmp2<-matrix(strsplit(input$patientGroups2, "\n")[[1]])
    maxlength <- max(length(tmp1), length(tmp2))
    length(tmp1) <- maxlength
    length(tmp2) <- maxlength
    
    # if(length(tmp1)!=length(tmp2)){
    #   sendSweetAlert(session, title = "Error", text = "Unbalanced patient groups. (Hint: each should include title in the first row)", type = "error",
    #                  btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    #   return()
    # }
    data = data.frame(cbind(tmp1, tmp2))
    colnames(data) <- c("Group 1", "Group 2")
    rownames(data) <- rep(1:dim(data)[1])
    
    all_patients <- as.vector(as.matrix(data))
    all_patients <- all_patients[all_patients != ""]
    all_patients <- all_patients[!is.na(all_patients)]
    timesOfPat <- table(all_patients)
    if(any(timesOfPat > 1)){
      get_patient_groups = paste(names(timesOfPat)[timesOfPat>1], "showed up multiple times")
      sendSweetAlert(session, title = "Error", text = "Same patient(s) discovered in two groups. Please make sure two groups of patients are mutually exclusive.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    else{
      get_patient_groups = data
    }
    return(get_patient_groups)
  })
  
  
  observeEvent(input$goAnalysisButton, {
    ## *** output user inputted patient groups ***
    output$inputtedPatientGroups <- DT::renderDataTable({ 
      if(input$patientGroups1 != "" & input$patientGroups1 != ""){
        out <- t(get_patient_groups())
        out <- cbind(apply(out, 1, function(x){paste(sum(!is.na(x!="")), "patients")}), out)
        out
      }
    }, selection="none",options=list(searching=F, ordering=F)) #,extensions = 'Responsive'
    
    if(length(DB) == 0){
      sendSweetAlert(session, title = "Error", text = "Insufficent input data.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    
    
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    progress$set(message = "Calculating... This may take a while", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    dataTypes <- list("0"= "mutation", "1"="rna", "2"="protein", "3"="clinical", "4"="image")
    dataTypes_upper <- list("0"= "Mutation", "1"="RNA", "2"="Protein", "3"="Clinical", "4"="Image Features")
    analysis.title = sprintf("Current Analysis Results: %s", dataTypes_upper[toString(input$AnalysisDataType)])
    
    output$analysisTitle <- renderUI({
      h4(analysis.title, style="color: STEELBLUE")
    })
    get_analysis_res <<- run_analysis(dataTypes[toString(input$AnalysisDataType)], get_patient_groups())
    get_analysis_res <<- data.frame(get_analysis_res)
    get_analysis_res <<- get_analysis_res[sort.list(get_analysis_res$pvalue), ]
    # save(get_analysis_res, file = "~/Desktop/get_analysis_res.Rdata")
    
    
    # output analysis result
    output$analysisResTable <- DT::renderDataTable({
      res <- get_analysis_res
    }, selection="none",options=list(searching=F, ordering=F)) #,extensions = 'Responsive'
    
    
  })
  # disease free survival
  output[["DFSurvivalPlot"]] <- renderPlot({
    PatList <- as.list(get_patient_groups())
    PatList <- lapply(PatList, function(x){setdiff(x, "")})
    clin_d <- data.frame(time=as.numeric(DB[["Clinical"]][unlist(PatList),"DiseaseFreeMonths"]),
                         event=DB[["Clinical"]][unlist(PatList),"DiseaseFreeStatus"]=="WITHTUMOR",
                         group = c(rep(names(PatList)[1], length(unlist(PatList[1]))), rep(names(PatList)[2], length(unlist(PatList[2]))))
    )
    clin_d <- clin_d[apply(clin_d, 1, function(x){!any(is.na(x))}),]
    survd <- survdiff(Surv(time, event, type="right") ~ group, data = clin_d)
    survf <- survfit(Surv(time,event) ~ group, data = clin_d)
    print(ggsurv(survf) + labs(title=paste("pvalue:", 1-pchisq(survd$chisq, 1)),
                               x='Time (Month)', y='Disease Free Survival'))
  }, height = 500, width = 700)
  
  
  #survival
  output[["SurvivalPlot"]] <- renderPlot({
    PatList <- as.list(get_patient_groups())
    PatList <- lapply(PatList, function(x){setdiff(x, "")})
    clin_d <- data.frame(time=as.numeric(DB[["Clinical"]][unlist(PatList),"OverallSurvivalMonths"]),
                         event=DB[["Clinical"]][unlist(PatList),"OverallSurvivalStatus"]=="DECEASED",
                         group = c(rep(names(PatList)[1], length(unlist(PatList[1]))), rep(names(PatList)[2], length(unlist(PatList[2]))))
    )
    clin_d <- clin_d[apply(clin_d, 1, function(x){!any(is.na(x))}),]
    survd <- survdiff(Surv(time, event, type="right") ~ group, data = clin_d)
    survf <- survfit(Surv(time,event) ~ group, data = clin_d)
    print(ggsurv(survf) + labs(title=paste("pvalue:", 1-pchisq(survd$chisq, 1)),
                               x='Time (Month)', y='Overall Survival'))
  }, height = 500, width = 700)
  
  output$dowloadAnalysisRes <- downloadHandler(
    filename = function() { "full_significant_genes.csv" },
    content = function(file) {
      write.csv(get_analysis_res,file, row.names=TRUE)
    })
  
  
  
}