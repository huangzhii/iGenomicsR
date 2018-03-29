
shinyServer(function(input, output, session) {
  library(data.table)
  library(DT)
  library(RColorBrewer)
  library(GenomicRanges)
  library(ggplot2)
  library(reshape2)
  library(gplots)
  library(plyr)
  #library("survplot")
  library("GGally")
  library("survival")
  source("utils.R")
  source("my_heatmap.R")
  source("my_analysis.R")
  source("gene_filter.R")
  options(shiny.maxRequestSize=30*1024^2)
  
  options(stringsAsFactors = FALSE)
  
  # *** The plot dimensions ***
  heightSize <- reactive ({ input$myHeight })
  widthSize <- reactive ({ input$myWidth })

  ########################################################################
  # data upload panel 
  ########################################################################
  dataL <- reactive({
    DB <- list()
    mySep<-switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="")
    if (!is.null(input$uploadMutation))  {
      inFile <- input$uploadMutation
      DB[["Mutation_gene"]] <- read.table(inFile$datapath, sep=mySep, header=TRUE, row.names = 1)
      colnames(DB[["Mutation_gene"]]) <- gsub(".", "-", colnames(DB[["Mutation_gene"]]), fixed = TRUE)
       
    }
    if (!is.null(input$uploadRNA))  {
      inFile <- input$uploadRNA
      DB[["RNA"]] <- read.table(inFile$datapath, sep=mySep, header=TRUE, row.names = 1)
      colnames(DB[["RNA"]]) <- gsub(".", "-", colnames(DB[["RNA"]]), fixed = TRUE)
      DB[["RNA"]] <- apply(DB[["RNA"]], c(1,2), as.numeric)
    }
    if (!is.null(input$uploadProtein))  {
      inFile <- input$uploadProtein
      DB[["Protein"]] <- read.table(inFile$datapath, sep=mySep, header=TRUE, row.names = 1)
      colnames(DB[["Protein"]]) <- gsub(".", "-", colnames(DB[["Protein"]]), fixed = TRUE)
      DB[["Protein"]] <- apply(DB[["Protein"]], c(1,2), as.numeric)
    }
    if (!is.null(input$uploadClin))  {
      inFile <- input$uploadClin
      DB[["Clinical"]] <- read.table(inFile$datapath, sep=mySep, header=TRUE, row.names = 1)
    }
    
    return(DB)
  })

  summarize_dataUpload <- eventReactive(input$uploadSummaryButton, {
    smartModal(error=F, title = "Summarizing Uploaded Data", content = "Summarizing uploaded data, please wait for a little while...")
    assign("DB", dataL(), envir = .GlobalEnv)
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
  
  output$dataUploadSummary <- renderTable({
    return(summarize_dataUpload())
  }, include.colnames=FALSE)
  
  sampleVenn_res <- reactive({
    vennData <- list()
    DB <- dataL()
    for(i in setdiff(names(DB), "Clinical")){
      vennData[[i]] <- colnames(DB[[i]])
    }
    if("Clinical" %in% names(DB)){
      vennData[["Clinical"]] <- rownames(DB[["Clinical"]]) 
    }
    venn(vennData)
  })
  output$sampleVenn <- renderPlot({
    sampleVenn_res()   
  }, height = 500, width = 500)
  
  sampleSelection <- observeEvent(input$sampleSelectionButton, {
    smartModal(F,"Selecting","Selecting samples shared by all subjects...")
    allDataTypes <- colnames(DB)
    if ("Clinical" %in% allDataTypes){
      samples <- rownames(DB[["Clinical"]])
      } 
    else if (length(setdiff(allDataTypes, "Clinical")) >=1 ){
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
    removeModal()
  })
  
  CatClin <- reactive({
    input$checkCatClin
  })
  QuanClin <- reactive({
    input$checkQuanClin
  })
  
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
  # Oncoplot
  mutation_OncoPlotInputGenes <- reactive({
    gsub("\\s","", strsplit(input$MutationInputGenes,",")[[1]])
  })
  mutation_OncoPlotClinInfo <- reactive({
    input$OncoPlotClin
  })
  mutation_OncoPlotProtein <- reactive({
    if(input$OncoPlotHasProtein){
      gsub("\\s","", strsplit(input$MutationInputPriteins,",")[[1]])
    }
  })
  mutation_OncoPlotRna <- reactive({
    if(input$OncoPlotHasRna){
      gsub("\\s","", strsplit(input$MutationInputRna,",")[[1]])
    }
  })
  
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
  
  OncoPlot_res <- reactive({
    my_heatmap_mutation(mutation_genes=mutation_OncoPlotInputGenes(), 
                        rna_genes=mutation_OncoPlotRna(), 
                        protein_genes=mutation_OncoPlotProtein(), 
                        clinical_lab=mutation_OncoPlotClinInfo(), 
                        order_by="mutation")
  })
  
  output$OncoPlot <- renderPlot({
    OncoPlot_res()[["plot"]]   
  }, height = heightSize, width = widthSize)
  # download ordered data for heatmap
  output$downloadOncoPlotData <- downloadHandler(
    filename = function() { "data_ordered_by_mutation.csv" },
    content = function(file) {
      write.csv(OncoPlot_res()[["table"]], file, row.names=TRUE)
    }) ###
  

  ########################################################################
  # image panel 
  ########################################################################
  # Imageheat
  Image_ImageheatImage <- reactive({
    gsub("\\s","", strsplit(input$ImageInputFeatures,",")[[1]])
  })
  Image_ImageheatClinInfo <- reactive({
    input$ImageheatClin
    })
  Image_ImageheatMutation <- reactive({
    if(input$ImageheatHasProtein){
      gsub("\\s","", strsplit(input$ImageInputMutations,",")[[1]])
    }
  })
  Image_ImageheatProtein <- reactive({
    if(input$ImageheatHasProtein){
      gsub("\\s","", strsplit(input$ImageInputPriteins,",")[[1]])
    }
  })
  Image_ImageheatRna <- reactive({
    if(input$ImageheatHasRna){
      gsub("\\s","", strsplit(input$ImageInputRna,",")[[1]])
    }
  })
  
  Imageheat_res <- reactive({
    my_heatmap_mutation(mutation_genes=Image_ImageheatMutation(), 
                        rna_genes=Image_ImageheatRna(), 
                        protein_genes=Image_ImageheatProtein(), 
                        clinical_lab=Image_ImageheatClinInfo(), 
                        image_features=Image_ImageheatImage(),
                        order_by="image")
  })
  
  output$Imageheat <- renderPlot({
    Imageheat_res()[["plot"]]   
  }, height = heightSize, width = widthSize)
  # download ordered data for heatmap
  output$downloadImageheatData <- downloadHandler(
    filename = function() { "data_ordered_by_mutation.csv" },
    content = function(file) {
      write.csv(Imageheat_res()[["table"]], file, row.names=TRUE)
    }) ###
  
    
  ########################################################################
  # rna expression panel 
  ########################################################################
  # gene expression clustering and heatmap
  rna_RNAheatInputGenes <- reactive({
    if(input$RNAheatIsInputGenes==0){
      gsub("\\s","", strsplit(input$RNAheatInputGenes,",")[[1]])
    }
  })
  rna_RNAheatGeneCriteria <- reactive({
    if(input$RNAheatIsInputGenes==1){
      strsplit(input$RNAheatGeneCutoff,"and")[[1]]
    }
  })
  rna_RNAheatClustPara <- reactive({
    clust_para <- list()
    if(input$RNAheatIsInputGenes == "1"){
      clust_para[["method"]] <- "hc"
    } else {
      clust_para[["method"]] <- c("hc", "km")[as.numeric(input$RNAheatClustMethod) + 1]
      if(clust_para[["method"]]=="km"){
        clust_para[["k"]] <- as.numeric(input$RNAheatKmeansK)
      }
    }
    clust_para
  })
  rna_RNAheatClinInfo <- reactive({
    input$RNAheatClin
    })
  rna_RNAheatProtein <- reactive({
    if(input$RNAheatHasProtein){
      gsub("\\s","", strsplit(input$RNAheatInputPriteins,",")[[1]])
    }
  })
  rna_RNAheatMutation <- reactive({
    if(input$RNAheatHasMutation){
      gsub("\\s","", strsplit(input$RNAheatInputMutation,",")[[1]])
    }
  })
  
  RNAheat_res <- reactive({
    my_heatmap_rna(mode=input$RNAheatIsInputGenes,
                   clust_para=rna_RNAheatClustPara(),
                   mutation_genes=rna_RNAheatMutation(), 
                   protein_genes=rna_RNAheatProtein(), 
                   clinical_lab=rna_RNAheatClinInfo(),
                   rna_criteria=rna_RNAheatGeneCriteria(), 
                   rna_genes=rna_RNAheatInputGenes())
  })
  output$RNAheat <- renderPlot({
    RNAheat_res()[["plot"]]
  }, height = heightSize, width = widthSize)
  # gene clustering dendrogram
  output$RNAdendro <- renderPlot({
    if(input$RNAheatClustMethod==0){
      plot(RNAheat_res()[["sample_order_res"]][["hc"]], cex=0.5)
    }
  }, height = heightSize, width = widthSize)
  # download ordered data for heatmap
  output$downloadRNAheatData <- downloadHandler(
    filename = function() { "data_ordered_by_rna.csv" },
    content = function(file) {
      write.csv(RNAheat_res()[["table"]], 
                file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # protein panel 
  ########################################################################
  # plot heatmap order by protein
  protein_ProteinheatInputGenes <- reactive({
    if(input$ProteinheatIsInputGenes==0){
      gsub("\\s","", strsplit(input$ProteinheatInputGenes,",")[[1]])
    }
  })
  protein_ProteinheatGeneCriteria <- reactive({
    if(input$ProteinheatIsInputGenes==1){
      strsplit(input$ProteinheatGeneCutoff,"and")[[1]]
    }
  })
  protein_ProteinheatClinInfo <- reactive({
    input$ProteinheatClin
    })
  protein_ProteinheatRNA <- reactive({
    if(input$ProteinheatHasRNA){
      gsub("\\s","", strsplit(input$ProteinheatInputRNA,",")[[1]])
    }
  })
  protein_ProteinheatMutation <- reactive({
    if(input$ProteinheatHasMutation){
      gsub("\\s","", strsplit(input$ProteinheatInputMutation,",")[[1]])
    }
  })
  
  Proteinheat_res <- reactive({
    my_heatmap_protein(mode=input$ProteinheatIsInputGenes,
                       mutation_genes=protein_ProteinheatMutation(), 
                       rna_genes=protein_ProteinheatRNA(), 
                       clinical_lab=protein_ProteinheatClinInfo(),
                       protein_criteria=protein_ProteinheatGeneCriteria(), 
                       protein_genes=protein_ProteinheatInputGenes())
  })
  output$Proteinheat <- renderPlot({
    Proteinheat_res()[["plot"]]
  }, height = heightSize, width = widthSize)
  # gene clustering dendrogram
  output$Proteindendro <- renderPlot({
    plot(Proteinheat_res()[["sample_order_res"]][["hc"]], cex=0.5)
  }, height = heightSize, width = widthSize)
  # download ordered data for heatmap
  output$downloadProteinheatData <- downloadHandler(
    filename = function() { "data_ordered_by_protein.csv" },
    content = function(file) {
      write.csv(Proteinheat_res()[["table"]], 
                file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # Clinical panel 
  ########################################################################
  # plot heatmap order by clinical feature
  clinical_ClinheatClinInfo <- reactive({
    input$ClinheatClin
  })
  clinical_ClinheatRNA <- reactive({
    if(input$ClinheatHasRNA){
      gsub("\\s","", strsplit(input$ClinheatInputRNA,",")[[1]])
    }
  })
  clinical_ClinheatProtein <- reactive({
    if(input$ClinheatHasProtein){
      gsub("\\s","", strsplit(input$ClinheatInputProtein,",")[[1]])
    }
  })
  clinical_ClinheatMutation <- reactive({
    if(input$ClinheatHasMutation){
      gsub("\\s","", strsplit(input$ClinheatInputMutation,",")[[1]])
    }
  })
  clinical_ClinheatOrderFeature <- reactive({
    input$ClinheatSelectOrderFeature
  })
  
  Clinheat_res <- reactive({
    my_heatmap_mutation(mutation_genes=clinical_ClinheatMutation(), 
                        rna_genes=clinical_ClinheatRNA(), 
                        protein_genes=clinical_ClinheatProtein(), 
                        clinical_lab=clinical_ClinheatClinInfo(),
                        order_by = "clinical", 
                        order_clin_feature=clinical_ClinheatOrderFeature())
  })
  output$Clinheat <- renderPlot({
    Clinheat_res()[["plot"]]
  }, height = heightSize, width = widthSize)
  # download ordered data for heatmap
  output$downloadClinheatData <- downloadHandler(
    filename = function() { "data_ordered_by_clinical.csv" },
    content = function(file) {
      write.csv( Clinheat_res()[["table"]], file, row.names=TRUE)
    }) ###
  
  
  ########################################################################
  # analysis panel 
  ########################################################################
  get_patient_groups <- reactive({
    if(input$patientGroups == "") {return(NULL)} 
    tmp<-matrix(strsplit(input$patientGroups, "\n")[[1]])
    myColnames<-strsplit(tmp[1], "\t")[[1]]
    data <- ldply(tmp[2:length(tmp)], function(x){
      gsub(" ", "", unlist(strsplit(x, "\t")))
    })
    colnames(data) <- myColnames
    all_patients <- as.vector(as.matrix(data))
    all_patients <- all_patients[all_patients != ""]
    timesOfPat <- table(all_patients)
    if(any(timesOfPat > 1)){
      return(paste(names(timesOfPat)[timesOfPat>1], "showed up multiple times"))
    }
    data
  })
  ## *** output user inputted patient groups ***
  output$inputtedPatientGroups <- renderTable({ 
    if(input$patientGroups != ""){
      out <- t(get_patient_groups())
      out <- cbind(apply(out, 1, function(x){paste(sum(x!=""), "patients")}), out)
      out
    }
  })
  get_analysis_res <- eventReactive(input$goAnalysisButton, {
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    progress$set(message = "Calculating... This may take a while", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    dataTypes <- list("0"= "mutation", "1"="rna", "2"="protein", "3"="clinical")
    run_analysis(dataTypes[input$AnalysisDataType], get_patient_groups())
  })
  # output analysis result
  output$analysisResTable <- renderTable({ 
    res <- get_analysis_res()
    if(nrow(res)>10){
      res[1:10, ]
    } else {
      res
    }
  }, include.rownames=TRUE)
  output$dowloadAnalysisRes <- downloadHandler(
    filename = function() { "full_significant_genes.csv" },
    content = function(file) {
      write.csv(get_analysis_res(),file, row.names=TRUE)
    }) 
  
  # disease free survival
  get_df_survival <- eventReactive(input$goAnalysisButton, {
    PatList <- as.list(get_patient_groups())
    PatList <- lapply(PatList, function(x){setdiff(x, "")})
    clin_d <- data.frame(time=as.numeric(DB[["Clinical"]][unlist(PatList),"DiseaseFreeMonths"]),
                         event=DB[["Clinical"]][unlist(PatList),"DiseaseFreeStatus"]=="WITHTUMOR",
                         group = c(rep(names(PatList)[1], length(unlist(PatList[1]))), rep(names(PatList)[2], length(unlist(PatList[2]))))
    )
    clin_d <- clin_d[apply(clin_d, 1, function(x){!any(is.na(x))}),]
    survd <- survdiff(Surv(time, event, type="right") ~ group, data = clin_d)
    survf <- survfit(Surv(time,event) ~ group, data = clin_d)
    print(ggsurv(survf) + labs(title=paste("pvalue:", pchisq(survd$chisq, 1)),
                               x='Time (Month)', y='Disease free survival'))
  })
  output[["DFSurvivalPlot"]] <- renderPlot({ get_df_survival()
  }, height = 500, width = 700)
  # survival
  get_survival <- eventReactive(input$goAnalysisButton, {
    PatList <- as.list(get_patient_groups())
    PatList <- lapply(PatList, function(x){setdiff(x, "")})
    clin_d <- data.frame(time=as.numeric(DB[["Clinical"]][unlist(PatList),"OverallSurvivalMonths"]),
                         event=DB[["Clinical"]][unlist(PatList),"OverallSurvivalStatus"]=="DECEASED",
                         group = c(rep(names(PatList)[1], length(unlist(PatList[1]))), rep(names(PatList)[2], length(unlist(PatList[2]))))
    )
    clin_d <- clin_d[apply(clin_d, 1, function(x){!any(is.na(x))}),]
    survd <- survdiff(Surv(time, event, type="right") ~ group, data = clin_d)
    survf <- survfit(Surv(time,event) ~ group, data = clin_d)
    print(ggsurv(survf) + labs(title=paste("pvalue:", pchisq(survd$chisq, 1)),
                               x='Time (Month)', y='Survival'))
  })
  output[["SurvivalPlot"]] <- renderPlot({ get_survival()
  }, height = 500, width = 700)
  
  
  ########################################################################
  # Data navigator panel 
  ########################################################################
  # output gene mutation test result
  output$geneMutationTestResTable <- renderTable({
    assign("DB", dataL(), envir = .GlobalEnv)
    print(head(DB[["Mutation_gene"]]))
    DB[["mutation_gene_test"]] <- run_gene_mutation_association(DB[["Mutation_gene"]])
    print(DB[["mutation_gene_test"]])
    d <- head(DB[["mutation_gene_test"]], n=min(10, nrow(DB[["mutation_gene_test"]])))
    d[,"oddsRatio"] <- format(as.numeric(d[,"oddsRatio"]),nsmall=2, digits=2)
    d[,"pvalue"] <- format(as.numeric(d[,"pvalue"]),scientific=TRUE, nsmall=2,digits=2)
    d[,"adj_pvalue"] <- format(as.numeric(d[,"adj_pvalue"]),scientific=TRUE, nsmall=2,digits=2)
    d
  },include.rownames=FALSE)
  output$dowloadGeneMutationTestRes <- downloadHandler(
    filename = function() { "mutation_gene_association_test.csv" },
    content = function(file) {
      write.csv(DB[["mutation_gene_test"]],file, row.names=FALSE)
    }) 
# output selected mutations
  output$selectedGeneMutationsTable <- renderTable({
    DB[["Mutation_gene"]][unlist(strsplit(gsub(" ", "", input$genesToPullMutation), ",", fixed = TRUE)),]
  },include.rownames=TRUE)
  
  
  # output dot plot for pair-wise expression profile
  dotplot_RNA_genes <- reactive({
    if(input$DataNavStartDataType == 1){
      gsub("\\s","", strsplit(input$RNAPairToDoPlot,",")[[1]])
    }
  })
  dotplot_protein_genes <- reactive({
    if(input$DataNavStartDataType == 2){
      gsub("\\s","", strsplit(input$ProteinPairToDoPlot,",")[[1]])
    }
  })
  output$RNADotPlot <- renderPlot({
    genes <- dotplot_RNA_genes()
    d <- DB[["RNA"]][genes,,drop=FALSE]
    d <- d[,!apply(d, 2, function(z){any(is.na(z))})]
    x <- d[genes[1],,drop=TRUE]
    y <- d[genes[2],,drop=TRUE]
    corr <- cor(x, y)
    plot(x, y, xlab=genes[1], ylab=genes[2], main="RNA expression", sub=paste("correlation:", corr))
  }, height = 500, width = 500)
  
  output[["ProteinDotPlot1"]] <- renderPlot({
    genes <- dotplot_protein_genes()
    d <- DB[["Protein"]][genes,,drop=FALSE]
    d <- d[,!apply(d, 2, function(z){any(is.na(z))})]
    d <- apply(d, c(1,2), as.numeric)
    x <- d[genes[1],,drop=TRUE]
    y <- d[genes[2],,drop=TRUE]
    corr <- cor(x, y)
    plot(x, y, xlab=genes[1], ylab=genes[2], main="Protein expression", sub=paste("correlation:", corr))
  }, height = 500, width = 500)
  
  output$ClinicalInfoTable <- DT::renderDataTable({
    DT::datatable(DB[["Clinical"]], extensions = 'Responsive', escape=F, selection = 'none', rownames = T)
    
  })
  
  get_all_clin <- reactive({
    DB <- dataL()
    ClinList <- list()
    for(i in colnames(DB[["Clinical"]])){
      ClinList[[i]] <- i
    }
    ClinList
  })
  output$checkCatClinUI <- renderUI({
    if (is.null(input$uploadClin))
      return()
    checkboxGroupInput('checkCatClin', 'Select categorical clinical feature', get_all_clin(), 
                       selected = "") 
  })
  output$checkQuanClinUI <- renderUI({
    if (is.null(input$uploadClin))
      return()
    checkboxGroupInput('checkQuanClin', 'Select quantitative clinical feature', get_all_clin(), 
                       selected = "") 
  })
  
})

