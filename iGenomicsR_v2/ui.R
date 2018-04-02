

shinyUI(pageWithSidebar(
  
  
  headerPanel(h3("iGenomicsR: an integrative platform to explore, visualize, and analyze multidimensional genomics data for disease"),
              tags$head(tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
                        tags$style(type="text/css", "select { max-width: 200px; }"),
                        tags$style(type="text/css", "textarea { max-width: 185px; }"),
                        tags$style(type="text/css", ".jslider { max-width: 200px; }"),
                        tags$style(type='text/css', ".well { max-width: 330px; }"),
                        #tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                        tags$style(type='text/css', ".span4 { max-width: 330px; }")
              )
  ),
  
  sidebarPanel(
    
    conditionalPanel(condition="input.tabs1=='Data upload'",
                     radioButtons("fileSepDF", "Delimiter:", list("Tab"=2,"Comma"=1,"Semicolon"=3)),
                     fileInput("uploadMutation", "Load mutation table", multiple = FALSE),
                     fileInput("uploadRNA", "Load RNA expression table", multiple = FALSE),
                     fileInput("uploadProtein", "load protein expression table", multiple = FALSE),
                     fileInput("uploadClin", "load clinical data", multiple = FALSE),
                     actionButton("uploadSummaryButton", "Summarize uploaded data")),
    conditionalPanel(condition="input.tabs1=='Data navigator'",
                     h4("Select data type"),
                     radioButtons("DataNavStartDataType", "", list("Mutation"=0, "RNA expression"=1, "Protein expression"=2, "Clinical data"=3)),
                     conditionalPanel(condition="input.DataNavStartDataType==0",
                                      checkboxInput("isPullMutation", "Pull out mutations for genes", FALSE),
                                      conditionalPanel(condition="input.isPullMutation",
                                                       textInput("genesToPullMutation", "Paste genes here", "PTEN, TP53, POLE"))
                     ),
                     conditionalPanel(condition="input.DataNavStartDataType==1",
                                      textInput("RNAPairToDoPlot", "Input two genes to do dot plot", "PTEN, TP53")
                     ),
                     conditionalPanel(condition="input.DataNavStartDataType==2",
                                      textInput("ProteinPairToDoPlot", "Input genes to do dot plot", "PTEN, ARID1B")
                     ),
                     conditionalPanel(condition="input.DataNavStartDataType==3")
    ),
    
    conditionalPanel(condition="input.tabs1=='Data integration'",
                     h4("Select data type to start"),
                     radioButtons("DataIntStartDataType", "", list("Mutation"=0, "Image"=4, "RNA expression"=1, "Protein expression"=2, "Clinical data"=3)),
                     checkboxInput("plotSize", "Adjust plot size", FALSE),
                     conditionalPanel(condition="input.plotSize",
                                      numericInput("myHeight", "Plot height:", value=800),
                                      numericInput("myWidth", "Plot width:", value=1200)),
                     conditionalPanel(condition="input.DataIntStartDataType==0",
                                      h4("de novo identification of associated genes or plot selected genes"),
                                      radioButtons("MutationIsInputGenes", "", list("Input genes"=0, "de nov"=1)),
                                      conditionalPanel(condition="input.MutationIsInputGenes==0",
                                                       textInput("MutationInputGenes", "Paste genes here", "PTEN, TP53")),
                                      h4("Add more features to heatmap"),
                                      checkboxInput("OncoPlotHasProtein", "Input genes to plot protein level", FALSE),
                                      conditionalPanel(condition="input.OncoPlotHasProtein==1",
                                                       textInput("MutationInputPriteins", "Paste genes here", "PTEN")),
                                      checkboxInput("OncoPlotHasRna", "Input genes to plot RNA level", FALSE),
                                      conditionalPanel(condition="input.OncoPlotHasRna==1",
                                                       textInput("MutationInputRna", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("OncoPlotHasClin", "Select clinical data to plot", FALSE),
                                      conditionalPanel(condition="input.OncoPlotHasClin",
                                                       uiOutput("OncoPlotClinUI")
                                                       )
                                      
                     ),
                     conditionalPanel(condition="input.DataIntStartDataType==1",
                                      h4("de novo clustering of whole transcriptome or selected genes"),
                                      radioButtons("RNAheatIsInputGenes", "", list("Input genes"=0, "de nov"=1), 1),
                                      conditionalPanel(condition="input.RNAheatIsInputGenes==0",
                                                       textInput("RNAheatInputGenes", "Paste genes here", "PTEN, TP53"),
                                                       radioButtons("RNAheatClustMethod", "", list("hierarchical clustering"=0, "kmeans clustering"=1), 0),
                                                       conditionalPanel(condition="input.RNAheatClustMethod==1",
                                                                        textInput("RNAheatKmeansK", "Number of clusters", "2"))),
                                      conditionalPanel(condition="input.RNAheatIsInputGenes==1",
                                                       textInput("RNAheatGeneCutoff", "Gene filter criteria", "var > 0.95"),
                                                       p("eg. maxExp > 0.5 and var > 0.8 and cv > 0.5\nFilter order: maxExp, var, cv. Check FAQ for more detail.")),
                                      h4("Add more features to heatmap"),
                                      checkboxInput("RNAheatHasProtein", "Input genes to plot protein level", FALSE),
                                      conditionalPanel(condition="input.RNAheatHasProtein==1",
                                                       textInput("RNAheatInputPriteins", "Paste genes here", "PTEN")),
                                      checkboxInput("RNAheatHasMutation", "Input genes to plot mutation", FALSE),
                                      conditionalPanel(condition="input.RNAheatHasMutation==1",
                                                       textInput("RNAheatInputMutation", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("RNAheatHasClin", "Select clinical data to plot", FALSE),
                                      conditionalPanel(condition="input.RNAheatHasClin",
                                                       checkboxGroupInput('RNAheatClin', '', c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]]), 
                                                                          selected = "")
                                                       )
                                      
                     ),
                     conditionalPanel(condition="input.DataIntStartDataType==2",
                                      h4("de novo clustering of whole proteome or selected genes"),
                                      radioButtons("ProteinheatIsInputGenes", "", list("Input genes"=0, "de nov"=1), 0),
                                      conditionalPanel(condition="input.ProteinheatIsInputGenes==0",
                                                       textInput("ProteinheatInputGenes", "Paste genes here", "PTEN, AATF")),
                                      conditionalPanel(condition="input.ProteinheatIsInputGenes==1",
                                                       textInput("ProteinheatGeneCutoff", "Gene filter criteria", "var > 0.95"),
                                                       p("eg. maxExp > 0.5 and var > 0.8 and cv > 0.5\nFilter order: maxExp, var, cv. Check FAQ for more detail.")),
                                      h4("Add more features to heatmap"),
                                      checkboxInput("ProteinheatHasRNA", "Input genes to plot RNA level", FALSE),
                                      conditionalPanel(condition="input.ProteinheatHasRNA==1",
                                                       textInput("ProteinheatInputRNA", "Paste genes here", "PTEN")),
                                      checkboxInput("ProteinheatHasMutation", "Input genes to plot mutation", FALSE),
                                      conditionalPanel(condition="input.ProteinheatHasMutation==1",
                                                       textInput("ProteinheatInputMutation", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("ProteinheatHasClin", "Select clinical data to plot", FALSE),
                                      conditionalPanel(condition="input.ProteinheatHasClin",
                                                       checkboxGroupInput('ProteinheatClin', '', c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]]), 
                                                                          selected = "")
                                                       )
                                      
                     ),
                     conditionalPanel(condition="input.DataIntStartDataType==3",
                                      h4("Select clinical features to plot"),
                                      checkboxGroupInput('ClinheatClin', '', c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]]), 
                                                         selected = c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]])),
                                      
                                      selectInput('ClinheatSelectOrderFeature', "Order samples by", c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]]), 
                                                         selected = c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]])[1]),
                                      h4("Add more features to heatmap"),
                                      checkboxInput("ClinheatHasRNA", "Input genes to plot RNA level", FALSE),
                                      conditionalPanel(condition="input.ClinheatHasRNA==1",
                                                       textInput("ClinheatInputRNA", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("ClinheatHasMutation", "Input genes to plot mutation", FALSE),
                                      conditionalPanel(condition="input.ClinheatHasMutation==1",
                                                       textInput("ClinheatInputMutation", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("ClinheatHasProtein", "Input genes to plot protein", FALSE),
                                      conditionalPanel(condition="input.ClinheatHasProtein==1",
                                                       textInput("ClinheatInputProtein", "Paste genes here", "PTEN"))
                                      
                     ),
                     conditionalPanel(condition="input.DataIntStartDataType==4",
                                      textInput("ImageInputFeatures", "Paste features here", "Fraction_Fiber, Nuclei_Density, Nuclei_Area"),
                                      h4("Add more features to heatmap"),
                                      checkboxInput("ImageheatHasMutation", "Input genes to plot mutation profile", FALSE),
                                      conditionalPanel(condition="input.ImageheatHasMutation==1",
                                                       textInput("ImageInputMutations", "Paste genes here", "TP53")),
                                      checkboxInput("ImageheatHasProtein", "Input genes to plot protein level", FALSE),
                                      conditionalPanel(condition="input.ImageheatHasProtein==1",
                                                       textInput("ImageInputPriteins", "Paste genes here", "TP53")),
                                      checkboxInput("ImageheatHasRna", "Input genes to plot RNA level", FALSE),
                                      conditionalPanel(condition="input.ImageheatHasRna==1",
                                                       textInput("ImageInputRna", "Paste genes here", "PTEN, TP53")),
                                      checkboxInput("ImageheatHasClin", "Select clinical data to plot", FALSE),
                                      conditionalPanel(condition="input.ImageheatHasClin",
                                                       checkboxGroupInput('ImageheatClin', '', c(DB[["Clinical_cat_lab"]], DB[["Clinical_quan_lab"]]), 
                                                                          selected = "")
                                                       )
                                      
                     )
                     
    ),
    
    conditionalPanel(condition="input.tabs1=='Data analysis'",
                     h4("Select analysis module"),
                     radioButtons("AnalysisDataType", "", list("Mutation"=0, "RNA expression"=1, "Protein expression"=2, "Clinical data"=3, "Survival"=4)),
                     h4("Define patient groups"),
                     tags$textarea(id="patientGroups", ncol=2),
                     h5(""),
                     actionButton("goAnalysisButton", "Run Analysis!"),
                     p("Click the button to start analysis.")
    )
  ),
  
  mainPanel(
    tabsetPanel(id="tabs1",
                # Welcome tab
                tabPanel("About",
                         HTML('<p> <br> This application was developed to help biologists to do integrative analysis of multi-dimentsional genomics data. We provide plots and statistical test to identify patterns from their data and also test their hypothesis. We would like to thank everyone who has made constructive suggestions so far. We will document the addition of new features in the News tab.</p>'),
                         h4("Software references"),
                         HTML('<p>R Development Core Team. <i><a href="http://www.r-project.org/">R</a>:  A Language and Environment for Statistical Computing.</i> R Foundation for Statistical Computing, Vienna (2013) <br>
                              RStudio and Inc. <i><a href="http://www.rstudio.com/shiny/">shiny</a>: Web Application Framework for R.</i> R package version 0.5.0 (2013) <br> 
                              Hadley Wickham. <i><a href="http://docs.ggplot2.org/current/">ggplot2</a>: a plotting system for R.</i> R package version 1.0.1 <br>
                              Erich Neuwirth. <i><a href="http://cran.r-project.org/web/packages/RColorBrewer/index.html">RColorBrewer</a>: ColorBrewer palettes.</i> R package version 1.1-2<br>
                              Hadley Wickham. <i><a href="https://cran.r-project.org/web/packages/reshape2/index.html">reshape2</a>: Flexibly Reshape Data.</i> R package version 1.4.1 <br>
                              David B. Dahl. <i><a href="https://cran.r-project.org/web/packages/xtable/index.html">xtable</a>: Export tables to LaTeX or HTML.</i> R package version 1.7-4 <br>
                              </p>')
                         ),
                
                # Data upload tab
                tabPanel("Data upload",
                         tableOutput("dataUploadSummary"),
                         conditionalPanel(condition="input.uploadClin!=''", 
                                          uiOutput("checkCatClinUI"),
                                          uiOutput("checkQuanClinUI")),
                         h4("venn diagram for sample overlaps between uploaded data types"),
                         plotOutput("sampleVenn", height='100%', width='100%'),
                         h5("Select samples shared by each of the data type is highly recommended"),
                         actionButton("sampleSelectionButton", "Select samples shared by all data types"),
                         h5("download data after sample selection"),
                         downloadButton("downloadSelectedMutation", "Download selected mutation data in .csv format"),
                         downloadButton("downloadSelectedRNA", "Download selected RNA data in .csv format"),
                         downloadButton("downloadSelectedProtein", "Download selected protein data in .csv format"),
                         downloadButton("downloadSelectedClin", "Download selected clinical data in .csv format")
                ),
                # Data navigator tab
                tabPanel("Data navigator",
                         conditionalPanel(condition="input.DataNavStartDataType==0", 
                                          h4("Genes mutation profile association test results"),
                                          tableOutput("geneMutationTestResTable"),
                                          downloadButton("dowloadGeneMutationTestRes", "Download full table of significant associations as .CSV file"),
                                          h4("Mutations for genes you selected"),
                                          tableOutput("selectedGeneMutationsTable")
                         ),
                         conditionalPanel(condition="input.DataNavStartDataType==1", 
                                          h4("RNA expression dotplot"),
                                          plotOutput("RNADotPlot", height='100%', width='100%')
                         ),
                         conditionalPanel(condition="input.DataNavStartDataType==2", 
                                          h4("Protein expression dotplot"),
                                          plotOutput("ProteinDotPlot1", height='100%', width='100%')
                         ),
                         conditionalPanel(condition="input.DataNavStartDataType==3", 
                                          h4("Clinical information"),
                                          DT::dataTableOutput("ClinicalInfoTable")
                         )
                ),
                
                # Boxplot tab
                tabPanel("Data integration",
                         
                         conditionalPanel(condition="input.DataIntStartDataType==0", 
                                          plotOutput("OncoPlot", height='100%', width='100%'),
                                          h4("Download data"),
                                          downloadButton("downloadOncoPlotData", "Download data for heatmap in .csv format")
                         ),
                         conditionalPanel(condition="input.DataIntStartDataType==4", 
                                          plotOutput("Imageheat", height='100%', width='100%'),
                                          h4("Download data"),
                                          downloadButton("downloadImageheatData", "Download data for heatmap in .csv format")
                         ),
                         conditionalPanel(condition="input.DataIntStartDataType==1", 
                                          plotOutput("RNAheat", height='100%', width='100%'),
                                          plotOutput("RNAdendro", height='100%', width='100%'),
                                          h4("Download data"),
                                          downloadButton("downloadRNAheatData", "Download data for heatmap in .csv format")
                         ),
                         conditionalPanel(condition="input.DataIntStartDataType==2", 
                                          plotOutput("Proteinheat", height='100%', width='100%'),
                                          plotOutput("Proteindendro", height='100%', width='100%'),
                                          h4("Download data"),
                                          downloadButton("downloadProteinheatData", "Download data for heatmap in .csv format")
                         ),
                         conditionalPanel(condition="input.DataIntStartDataType==3", 
                                          plotOutput("Clinheat", height='100%', width='100%'),
                                          h4("Download data"),
                                          downloadButton("downloadClinheatData", "Download data for heatmap in .csv format")
                         )
                ), 
                
                # Data analysis tab
                tabPanel("Data analysis",
                         h4("Patients you inputted"),
                         tableOutput("inputtedPatientGroups"),
                         conditionalPanel(condition="input.AnalysisDataType==4",
                                          plotOutput("SurvivalPlot", height='100%', width='100%'),
                                          plotOutput("DFSurvivalPlot", height='100%', width='100%')),
                         conditionalPanel(condition="input.AnalysisDataType!=4",
                                          h4("Associated genes or clinical features"),
                                          tableOutput("analysisResTable"),
                                          downloadButton("dowloadAnalysisRes", "Download full table of significant genes as .CSV file")
                         )
                ),
                
                # News
                tabPanel("News",
                         h5("April 28, 2016"), 
                         p("This is a test version")
                ),			
                
                # FAQ 
                tabPanel("FAQ",
                         "FAQ"
                )
    ),
    h6("This application was created by the ", a("Xing Tang", href="https://scholar.google.com/citations?user=F9arEIMAAAAJ&hl=en"), 
       " from The Ohio State University. Please send bugs and feature requests to Xing Tang (tangx1986(at)gmail.com) or (tang.811(at)osu.edu). This application uses the ", 
       a("shiny package from RStudio", href="http://www.rstudio.com/shiny/"), ".")
  )
  )
  )







