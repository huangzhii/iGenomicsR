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

DB <- list()
all.file.uploaded <- F
function(input, output, session) {
  output$check1 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check2 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check3 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  output$check4 <- renderText({'<img src="./images/check_no.png", style="width:30px">'})
  
  ########################################################################
  # data upload panel 
  ########################################################################
  observeEvent(input$action1,{
    if(!is.null(DB[["Mutation_gene"]]) & !is.null(DB[["RNA"]]) & !is.null(DB[["Protein"]]) & !is.null(DB[["Clinical"]])){
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
      sendSweetAlert(session, title = "Insufficient Input Data", text = "Please upload required files.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    }
  })
  
  
  loadData.mutation <- function(){
    if (is.null(input$csvfile_mutation)){
      return(NULL)}else{
        DB[["Mutation_gene"]] <<- read.table(input$csvfile_mutation$datapath, sep=input$sep, header=TRUE, row.names = 1)
        colnames(DB[["Mutation_gene"]]) <<- gsub(".", "-", colnames(DB[["Mutation_gene"]]), fixed = TRUE)
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
  observeEvent(loadData.mutation(),{output$check1 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.mRNA(),{output$check2 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.protein(),{output$check3 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  observeEvent(loadData.clinical(),{output$check4 <- renderText({'<img src="./images/check_yes.png", style="width:30px">'})})
  
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
        d[[i]] <- c(paste("# of NA in", i), sum(toupper(DB[["Clinical"]][,i]) == "NA" | DB[["Clinical"]][,i]==""))
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
      DB[[i]] <- DB[[i]][,samples]
    }
    if ("Clinical" %in% allDataTypes){
      DB[["Clinical"]] <- DB[[i]][samples,]
    } 
    DB[["Clinical_cat_lab"]] <- CatClin()
    DB[["Clinical_quan_lab"]] <- QuanClin()
    assign("DB", DB, envir = .GlobalEnv)
    print(names(DB))
    
    print(head(DB[["Clinical"]]))
    print(head(DB[["Clinical_cat_lab"]]))
    print(head(DB[["Clinical_quan_lab"]]))
    
    sendSweetAlert(session, title = "Success", text = "Now only keep samples shared by all subjects", type = "success",
                   btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    session$sendCustomMessage("download_seleted_data_ready", "lalala")
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
  
  
  ########################################################################
  # mutation panel 
  ########################################################################
  
  observeEvent(input$action2,{
    if(!is.null(DB[["Mutation_gene"]]) & !is.null(DB[["RNA"]]) & !is.null(DB[["Protein"]]) & !is.null(DB[["Clinical"]])){
      session$sendCustomMessage("buttonCallbackHandler", "tab2")
    }
    else{
      sendSweetAlert(session, title = "Insufficient Input Data", text = "Please upload required files.", type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    }
  })
  
  
  
  
  
  
  
}