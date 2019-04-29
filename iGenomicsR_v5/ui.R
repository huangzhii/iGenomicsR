# Zhi Huang 06/09/2018
library(shinyWidgets)
library(markdown)
library(shinyjs)

navbarPage(title=div(a(img(src="images/iGenomicsR_logo2.png",
                           height = 28,
                           style = "margin:0px 0px; padding-bottom: 5px"), href=""),escape=F),
           tabPanel("Welcome",
                    sidebarLayout(
                      position = "left",
                      sidebarPanel(width = 3,
                                   h4("File Uploader", style="color: STEELBLUE"),
                                   actionButton("action_load_example", "Load Example"),
                                   useShinyjs(),
                                   div(
                                     id = "upload_panel",
                                     fileInput("csvfile_mutation", "Mutation Table",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv", ".xlsx", ".xls")),
                                     fileInput("csvfile_image", "Image Profile",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv", ".xlsx", ".xls")),
                                     fileInput("csvfile_mRNA", "RNA Expression Table",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv", ".xlsx", ".xls")),
                                     fileInput("csvfile_protein", "Protein Expression Table",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv", ".xlsx", ".xls")),
                                     fileInput("csvfile_clinical", "Clinical Profile (* Required)",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv", ".xlsx", ".xls")),
                                     
                                     # Include clarifying text ----
                                     helpText("Note: Maximum file size allowed for uploading is 300MB."),
                                     
                                     # Input: Checkbox if file has header ----
                                     checkboxInput("header", "Header", TRUE),
                                     
                                     fluidRow(
                                       # Input: Select separator ----
                                       column(6, radioButtons("sep", "Separator",
                                                              choices = c(Comma = ",",
                                                                          Semicolon = ";",
                                                                          Tab = "\t",
                                                                          Space = " "),
                                                              selected = ",")),
                                       # Input: Select quotes ----
                                       column(6, radioButtons("quote", "Quote",
                                                              choices = c(None = "",
                                                                          "Double Quote" = '"',
                                                                          "Single Quote" = "'"),
                                                              selected = '"'))
                                     ),
                                     # Horizontal line ----
                                     tags$hr(),
                                     p('If you want a sample .csv file to upload,',
                                       'you can first download the sample',
                                       a(href = 'data/Protein.csv', 'Protein.csv'), ', ',
                                       a(href = 'data/mutation.csv', 'mutation.csv'), ', ',
                                       a(href = 'data/RNA.csv', 'RNA.csv'), ', ',
                                       a(href = 'data/Clinical.csv', 'Clinical.csv'), ' and ',
                                       a(href = 'data/Image_Features_Ass_General_CPTAC_merged_by_mean.csv', 'image_features.csv' ),
                                       'files, and then try uploading them.')
                                     # actionButton("action1", "Confirm when Complete")
                                   )
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("About",
                                   h2("iGenomicsR:", style="color: STEELBLUE; font-size: 22px"),
                                   h2("A Multi-omics Platform for Integrative Genomics Data Analysis and Visualization", style="color: STEELBLUE; font-size: 20px; margin: 0px"),
                                   HTML('<p> <br> This application was developed to help biologists to do integrative analysis of multi-dimentsional genomics data. We provide plots and statistical test to identify patterns from their data and also test their hypothesis. We would like to thank everyone who has made constructive suggestions so far. We will document the addition of new features in the News tab.</p>'),
                                   img(src="images/iGenomicsR_banner.png",
                                       height = 600,
                                       style = "margin:0px 0px; padding-bottom: 5px")
                                   ),
                          tabPanel("Dataset",
                                   h4("Data File required", style="color: STEELBLUE; padding-top: 10px"),
                                   
                                   fluidRow(
                                     column(2, "Mutation Profile", htmlOutput("check1")),
                                     column(2, "Image Feature", htmlOutput("check5")),
                                     column(2, "RNA Expression", htmlOutput("check2")),
                                     column(2, "Protein Expression", htmlOutput("check3")),
                                     column(2, "Clinical Profile", htmlOutput("check4")),
                                     style="text-align: center"
                                   ),
                                   conditionalPanel(condition="output.data_ready_flag == 1",
                                                    actionButton("uploadSummaryButton", "Summarize Data"),
                                                    tableOutput("dataUploadSummary"),
                                                    h4("Verify Categorical and Numerical Clinical Profile", style="color: STEELBLUE; padding-top: 10px"),
                                                    fluidRow(
                                                      column(6, uiOutput("checkCatClinUI"),
                                                             helpText("Please check if they are categorical Data.")),
                                                      column(6, uiOutput("checkQuanClinUI"),
                                                             helpText("Please check if they are numerical Data."))
                                                    ),
                                                    h4("Venn diagram for sample overlaps between uploaded data types", style="color: STEELBLUE; padding-top: 10px"),
                                                    
                                                    tags$div(
                                                      plotOutput("sampleVenn", height='100%', width='100%'),
                                                      style="text-align: center"
                                                    ),
                                                    h5("Select samples shared by each of the data type is highly recommended", style="color: STEELBLUE; padding-top: 10px"),
                                                    actionButton("sampleSelectionButton", "Select samples shared by all data types"),
                                                    
                                                    h5("download data after sample selection"),
                                                    downloadButton("downloadSelectedMutation", "Download selected mutation data in .csv format"),
                                                    br(),
                                                    downloadButton("downloadSelectedImageFeature", "Download selected image feature data in .csv format"),
                                                    br(),
                                                    downloadButton("downloadSelectedRNA", "Download selected RNA data in .csv format"),
                                                    br(),
                                                    downloadButton("downloadSelectedProtein", "Download selected protein data in .csv format"),
                                                    br(),
                                                    downloadButton("downloadSelectedClin", "Download selected clinical data in .csv format"),
                                                    tags$hr(),
                                                    actionButton("action2", "Proceed to Data Navigator",style="color: WHITE; background-color: DODGERBLUE")
                                   ),
                                   tags$head(
                                     tags$script(HTML("(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
                                                      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
                                                      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
                                                      
                                                      ga('create', 'UA-113406500-2', 'auto');
                                                      ga('send', 'pageview');"))
                                     ),
                                   tags$head(tags$script('Shiny.addCustomMessageHandler("buttonCallbackHandler",
                                                         function(typeMessage) {console.log(typeMessage)
                                                         if(typeMessage == "tab1"){
                                                         console.log("got here");
                                                         $("a:contains(Dataset)").click();
                                                         }
                                                         if(typeMessage == "tab2"){
                                                         $("a:contains(Data Navigator)").click();
                                                         }
                                                         });
                                                         
                                                         // disable download at startup.
                                                         $(document).ready(function() {
                                                         $("#downloadSelectedMutation").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedImageFeature").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedRNA").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedProtein").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedClin").attr("disabled", "true").attr("onclick", "return false;");
                                                         
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready_mutation", function(message) {
                                                         $("#downloadSelectedMutation").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready_image", function(message) {
                                                         $("#downloadSelectedImageFeature").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready_rna", function(message) {
                                                         $("#downloadSelectedRNA").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready_protein", function(message) {
                                                         $("#downloadSelectedProtein").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready_clinical", function(message) {
                                                         $("#downloadSelectedClin").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         })
                                                         '))
                                   ) # end of mainPanel
                                   )
                                   )
                                   ) # end of sidebarLayout
                    
                                   ),
           tabPanel("Data Navigator",
                    navlistPanel(widths = c(2, 10),
                                 id = "data.navigator",
                                 tabPanel("Mutation",
                                    conditionalPanel(condition="output.hasMutationData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h4("Mutations for Genes List", style="color: STEELBLUE"),
                                              textAreaInput(inputId="MutationInputGenes",label="Paste genes here:",
                                                            value = "PTEN, TP53, POLE", height = 200),
                                              actionButton("action.navigator.mutation", "Run",style="color: WHITE; background-color: DODGERBLUE")
                                            ),
                                            mainPanel(
                                              h4("Genes mutation profile association test results", style="color: STEELBLUE"),
                                              DT::dataTableOutput("geneMutationTestResTable"),
                                              downloadButton("dowloadGeneMutationTestRes", "Download full table of significant associations as .CSV file"),
                                              h4("Mutations for genes you selected", style="color: STEELBLUE; padding-top: 10px"),
                                              DT::dataTableOutput("selectedGeneMutationsTable")
                                            )
                                          )
                                 )
                                 ),
                             tabPanel("Image features",
                                  conditionalPanel(condition="output.hasImageData == 1",
                                          h4("Image Features Table", style="color: STEELBLUE"),
                                          DT::dataTableOutput("ImageFeaturesTable"),
                                          h4("Image Features Heatmap", style="color: STEELBLUE; padding-top: 10px"),
                                          plotOutput("ImageFeaturesNavigatorPlot", height='600px', width='1000px')
                                 )),
                             tabPanel("RNA expression",
                                  conditionalPanel(condition="output.hasRNAData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h4("Input two genes to perform dot plot", style="color: STEELBLUE"),
                                              fluidRow(
                                                column(6, textInput(inputId="navigator.RNA.expression.gene.1",label="Gene 1:",
                                                                    value = "PTEN")),
                                                column(6, textInput(inputId="navigator.RNA.expression.gene.2",label="Gene 2:",
                                                                    value = "TP53"))
                                              ),
                                              actionButton("action.navigator.RNA", "Plot",style="color: WHITE; background-color: DODGERBLUE")
                                            ),
                                            mainPanel(
                                              h4("RNA expression dotplot", style="color: STEELBLUE"),
                                              plotOutput("RNADotPlot", height='100%', width='100%')
                                            )
                                          )
                                 )),
                             tabPanel("Protein expression",
                                  conditionalPanel(condition="output.hasProteinData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h4("Input two genes to perform dot plot", style="color: STEELBLUE"),
                                              fluidRow(
                                                column(6, textInput(inputId="navigator.protein.expression.gene.1",label="Gene 1:",
                                                                    value = "PTEN")),
                                                column(6, textInput(inputId="navigator.protein.expression.gene.2",label="Gene 2:",
                                                                    value = "ARID1B"))
                                              ),
                                              actionButton("action.navigator.protein", "Plot",style="color: WHITE; background-color: DODGERBLUE")
                                            ),
                                            mainPanel(
                                              h4("Protein expression dotplot", style="color: STEELBLUE"),
                                              plotOutput("ProteinDotPlot1", height='100%', width='100%')
                                            )
                                          )
                                 )),
                                 tabPanel("Clinical data",
                                          h4("Clinical Information", style="color: STEELBLUE"),
                                          DT::dataTableOutput("ClinicalInfoTable")
                                 )
                    )
           ),
           tabPanel("Data Integration",
                    navlistPanel(widths = c(2, 10),
                                 id = "data.integration",
                                 tabPanel("Mutation",
                                      conditionalPanel(condition="output.hasMutationData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              # h3("Choose Method", style="color: STEELBLUE"),
                                              # h4("De novo identification of associated genes", style="color: STEELBLUE"),
                                              # actionButton("action.integration.mutation.denovo", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              
                                              
                                              h4("Plot selected genes", style="color: STEELBLUE"),
                                              textAreaInput(inputId="genesToPullMutation",label="Paste genes here:",
                                                            value = "PTEN, TP53", height = 200),
                                              materialSwitch(inputId = "do_hclust_mutation", label = "Sort mutation genes decreasingly based on its score (number of mutates)", status = "primary"),
                                              tags$hr(),
                                              h4("Add more features to heatmap", style="color: STEELBLUE"),
                                              conditionalPanel(condition="output.hasImageData == 1",
                                                prettyCheckbox(inputId = "OncoPlotHasImage", label = "Select features to plot image data", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.OncoPlotHasImage==1",
                                                               uiOutput("ImageInputFeaturesSubUI_Mutation"),
                                                               materialSwitch(inputId = "do_hclust_mutation_image", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasProteinData == 1",
                                                prettyCheckbox(inputId = "OncoPlotHasProtein", label = "Input genes to plot protein level", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.OncoPlotHasProtein==1",
                                                               textInput("MutationInputProteins", "Paste genes here", "PTEN"),
                                                               materialSwitch(inputId = "do_hclust_mutation_protein", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasRNAData == 1",
                                                prettyCheckbox(inputId = "OncoPlotHasRna", label = "Input genes to plot RNA level", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.OncoPlotHasRna==1",
                                                               textInput("MutationInputRna", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_mutation_rna", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              prettyCheckbox(inputId = "OncoPlotHasClin", label = "Select clinical data to plot", icon = icon("check")),
                                              conditionalPanel(condition="input.OncoPlotHasClin",
                                                               uiOutput("OncoPlotClinUI")),
                                              tags$hr(),
                                              sliderInput(inputId="myHeight1", label="Plot height:", min=200, max=2000, value=800),
                                              sliderInput(inputId="myWidth1", label="Plot width:", min=200, max=2000, value=800),
                                              actionButton("action.integration.mutation.inputgenes", "Run",style="color: WHITE; background-color: DODGERBLUE")
                                            ),
                                            mainPanel(
                                              h4("Onco plot", style="color: STEELBLUE"),
                                              plotOutput("OncoPlot", height='100%', width='100%'),
                                              h4("Download onco plot data", style="color: STEELBLUE"),
                                              downloadButton("downloadOncoPlotData", "Download data for heatmap in .csv format")
                                            )
                                          )
                                      )
                                 ),
                                 tabPanel("Image features",
                                      conditionalPanel(condition="output.hasImageData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h4("Input Image Features", style="color: STEELBLUE"),
                                              uiOutput("ImageInputFeaturesUI"),
                                              materialSwitch(inputId = "do_hclust_image", label = "Sort image features by hierarchical clustering (average distance)", status = "primary"),
                                              h4("Add more features to heatmap", style="color: STEELBLUE"),
                                              conditionalPanel(condition="output.hasMutationData == 1",
                                                prettyCheckbox("ImageheatHasMutation", "Input genes to plot mutation profile", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ImageheatHasMutation==1",
                                                               textInput("ImageInputMutations", "Paste genes here", "TP53"),
                                                               materialSwitch(inputId = "do_hclust_image_mutation",
                                                                              label = "Sort genes decreasingly based on its score (number of mutates)",
                                                                              status = "primary")),
                                              conditionalPanel(condition="output.hasProteinData == 1",
                                                prettyCheckbox("ImageheatHasProtein", "Input genes to plot protein level", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ImageheatHasProtein==1",
                                                               textInput("ImageInputProteins", "Paste genes here", "PTEN"),
                                                               materialSwitch(inputId = "do_hclust_image_protein", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasRNAData == 1",
                                                prettyCheckbox("ImageheatHasRna", "Input genes to plot RNA level", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ImageheatHasRna==1",
                                                               textInput("ImageInputRna", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_image_rna", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              prettyCheckbox("ImageheatHasClin", "Select clinical data to plot", FALSE, icon = icon("check")),
                                              conditionalPanel(condition="input.ImageheatHasClin",
                                                               uiOutput("ImageheatClinUI")),
                                              tags$hr(),
                                              sliderInput(inputId="myHeight2", label="Plot height:", min=200, max=2000, value=800),
                                              sliderInput(inputId="myWidth2", label="Plot width:", min=200, max=2000, value=800),
                                              actionButton("action.integration.image.inputgenes", "Run",style="color: WHITE; background-color: DODGERBLUE")
                                            ),
                                            mainPanel(
                                              h4("Image heatmap", style="color: STEELBLUE"),
                                              plotOutput("Imageheat", height='100%', width='100%'),
                                              h4("Download image heatmap data", style="color: STEELBLUE"),
                                              downloadButton("downloadImageheatData", "Download data for heatmap in .csv format")
                                            )
                                          )
                                      )
                                 ),
                                 tabPanel("RNA expression",
                                      conditionalPanel(condition="output.hasRNAData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h3("Choose Method", style="color: STEELBLUE"),
                                              h4("De novo clustering of whole transcriptome", style="color: STEELBLUE"),
                                              textInput("RNAheatGeneCutoff", "Gene filter criteria", "var > 0.95"),
                                              helpText("eg. maxExp > 0.5 and var > 0.8 and cv > 0.5\nFilter order: maxExp, var, cv. Check FAQ for more detail."),
                                              prettyCheckbox(inputId = "show.RNA.name.1", label = "Show RNA names", FALSE, icon = icon("check")),
                                              actionButton("action.integration.RNA.denovo", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              
                                              h4("Clustering on selected genes", style="color: STEELBLUE"),
                                              textAreaInput(inputId="RNAheatInputGenes",label="Paste genes here:",
                                                            value = "PTEN, TP53", height = 200),
                                              awesomeRadio("RNAheatClustMethod", "", list("hierarchical clustering"=0, "kmeans clustering"=1), 0),
                                              conditionalPanel(condition="input.RNAheatClustMethod==1",
                                                               numericInput("RNAheatKmeansK", "Number of clusters", value = 2, min=2)),
                                              prettyCheckbox(inputId = "show.RNA.name.2", label = "Show RNA names", FALSE, icon = icon("check")),
                                              actionButton("action.integration.RNA.inputgenes", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              tags$hr(),
                                              h4("Add more features to heatmap", style="color: STEELBLUE"),
                                              conditionalPanel(condition="output.hasImageData == 1",
                                                prettyCheckbox(inputId = "RNAheatHasImage", label = "Select features to plot image data", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.RNAheatHasImage==1",
                                                               uiOutput("ImageInputFeaturesSubUI_RNA"),
                                                               materialSwitch(inputId = "do_hclust_rna_image", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasProteinData == 1",
                                                prettyCheckbox(inputId = "RNAheatHasProtein", label = "Input genes to plot protein level", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.RNAheatHasProtein==1",
                                                               textInput("RNAheatInputProteins", "Paste genes here", "PTEN"),
                                                               materialSwitch(inputId = "do_hclust_rna_protein", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasMutationData == 1",
                                                prettyCheckbox(inputId = "RNAheatHasMutation", label = "Input genes to plot mutation", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.RNAheatHasMutation==1",
                                                               textInput("RNAheatInputMutation", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_rna_mutation", label = "Sort genes decreasingly based on its score (number of mutates)", status = "primary")),
                                              prettyCheckbox(inputId = "RNAheatHasClin", label = "Select clinical data to plot", icon = icon("check")),
                                              conditionalPanel(condition="input.RNAheatHasClin",
                                                               uiOutput("RNAheatClinUI")),
                                              tags$hr(),
                                              sliderInput(inputId="myHeight3", label="Plot height:", min=200, max=2000, value=800),
                                              sliderInput(inputId="myWidth3", label="Plot width:", min=200, max=2000, value=800)
                                            ),
                                            mainPanel(
                                              h4("RNA heatmap", style="color: STEELBLUE"),
                                              plotOutput("RNAheat", height='100%', width='100%'),
                                              h4("RNA dendrogram", style="color: STEELBLUE"),
                                              plotOutput("RNAdendro", height='100%', width='100%'),
                                              h4("Download RNA heatmap data", style="color: STEELBLUE"),
                                              downloadButton("downloadRNAheatData", "Download data for heatmap in .csv format")
                                            )
                                          )
                                      )
                                 ),
                                 tabPanel("Protein expression",
                                      conditionalPanel(condition="output.hasProteinData == 1",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h3("Choose Method", style="color: STEELBLUE"),
                                              h4("De novo clustering of whole proteome", style="color: STEELBLUE"),
                                              textInput("ProteinheatGeneCutoff", "Gene filter criteria", "var > 0.95"),
                                              helpText("eg. maxExp > 0.5 and var > 0.8 and cv > 0.5\nFilter order: maxExp, var, cv. Check FAQ for more detail."),
                                              prettyCheckbox(inputId = "show.protein.name", label = "Show protein names", FALSE, icon = icon("check")),
                                              actionButton("action.integration.protein.denovo", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              h4("Clustering of selected genes", style="color: STEELBLUE"),
                                              h5("Hierarchical Clustering (Average)", style="color: STEELBLUE"),
                                              textAreaInput(inputId="ProteinheatInputGenes",label="Paste genes here:",
                                                            value = "PTEN, AATF", height = 200),
                                              actionButton("action.integration.protein.inputgenes", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              
                                              tags$hr(),
                                              h4("Add more features to heatmap", style="color: STEELBLUE"),
                                              conditionalPanel(condition="output.hasImageData == 1",
                                                prettyCheckbox(inputId = "ProteinheatHasImage", label = "Select features to plot image data", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ProteinheatHasImage==1",
                                                               uiOutput("ImageInputFeaturesSubUI_Protein"),
                                                               materialSwitch(inputId = "do_hclust_protein_image", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasRNAData == 1",
                                                prettyCheckbox("ProteinheatHasRNA", "Input genes to plot RNA level", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ProteinheatHasRNA==1",
                                                               textInput("ProteinheatInputRNA", "Paste genes here", "PTEN"),
                                                               materialSwitch(inputId = "do_hclust_protein_rna", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasMutationData == 1",
                                                prettyCheckbox("ProteinheatHasMutation", "Input genes to plot mutation", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ProteinheatHasMutation==1",
                                                               textInput("ProteinheatInputMutation", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_protein_mutation", label = "Sort genes decreasingly based on its score (number of mutates)", status = "primary")),
                                              prettyCheckbox("ProteinheatHasClin", "Select clinical data to plot", FALSE, icon = icon("check")),
                                              conditionalPanel(condition="input.ProteinheatHasClin",
                                                               uiOutput("ProteinheatClinUI")),
                                              tags$hr(),
                                              sliderInput(inputId="myHeight4", label="Plot height:", min=200, max=2000, value=800),
                                              sliderInput(inputId="myWidth4", label="Plot width:", min=200, max=2000, value=800)
                                            ),
                                            mainPanel(
                                              h4("Protein Heatmap", style="color: STEELBLUE"),
                                              plotOutput("Proteinheat", height='100%', width='100%'),
                                              h4("Protein Dentrogram", style="color: STEELBLUE"),
                                              plotOutput("Proteindendro", height='100%', width='100%'),
                                              h4("Download protein heatmap data", style="color: STEELBLUE"),
                                              downloadButton("downloadProteinheatData", "Download data for heatmap in .csv format")
                                            )
                                          )
                                      )
                                 ),
                                 tabPanel("Clinical data",
                                          sidebarLayout(
                                            position = "left",
                                            sidebarPanel(
                                              width=3,
                                              h3("Select clinical features to plot", style="color: STEELBLUE"),
                                              h4("Please check following selctions:", style="color: STEELBLUE"),
                                              uiOutput("ClinheatClinUI"),
                                              uiOutput("ClinheatSelectOrderFeatureUI"),
                                              tags$hr(),
                                              h4("Add more features to heatmap", style="color: STEELBLUE"),
                                              conditionalPanel(condition="output.hasImageData == 1",
                                                prettyCheckbox(inputId = "ClinheatHasImage", label = "Select features to plot image data", icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ClinheatHasImage==1",
                                                               uiOutput("ImageInputFeaturesSubUI_Clinical"),
                                                               materialSwitch(inputId = "do_hclust_clin_image", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasRNAData == 1",
                                                prettyCheckbox("ClinheatHasRNA", "Input genes to plot RNA level", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ClinheatHasRNA==1",
                                                               textInput("ClinheatInputRNA", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_clin_rna", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              conditionalPanel(condition="output.hasMutationData == 1",
                                                prettyCheckbox("ClinheatHasMutation", "Input genes to plot mutation", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ClinheatHasMutation==1",
                                                               textInput("ClinheatInputMutation", "Paste genes here", "PTEN, TP53"),
                                                               materialSwitch(inputId = "do_hclust_clin_mutation", label = "Sort genes decreasingly based on its score (number of mutates)", status = "primary")),
                                              conditionalPanel(condition="output.hasProteinData == 1",
                                                prettyCheckbox("ClinheatHasProtein", "Input genes to plot protein", FALSE, icon = icon("check"))
                                              ),
                                              conditionalPanel(condition="input.ClinheatHasProtein==1",
                                                               textInput("ClinheatInputProtein", "Paste genes here", "PTEN"),
                                                               materialSwitch(inputId = "do_hclust_clin_protein", label = "Sort genes by hierarchical clustering (average distance)", status = "primary")),
                                              actionButton("action.integration.clinical", "Run",style="color: WHITE; background-color: DODGERBLUE"),
                                              tags$hr(),
                                              sliderInput(inputId="myHeight5", label="Plot height:", min=200, max=2000, value=800),
                                              sliderInput(inputId="myWidth5", label="Plot width:", min=200, max=2000, value=800),
                                              helpText("Note: Mutation genes are being sorted decreasingly based on its score (number of mutates).")
                                              
                                            ),
                                            mainPanel(
                                              h4("Clinical heatmap", style="color: STEELBLUE"),
                                              plotOutput("Clinheat", height='100%', width='100%'),
                                              h4("Download clinical heatmap data", style="color: STEELBLUE"),
                                              downloadButton("downloadClinheatData", "Download data for heatmap in .csv format")
                                            )
                                          )
                                 )
                    )
           ),
           tabPanel("Data Analysis",
                    sidebarLayout(
                      position = "left",
                      sidebarPanel(
                        width=3,
                        h4("Select analysis module", style="color: STEELBLUE"),
                        uiOutput("AnalysisDataTypeUI"),
                        h4("Define patient groups", style="color: STEELBLUE"),
                        fluidRow(
                          column(6, uiOutput("PatientGroupsInputUI1")),
                          column(6, uiOutput("PatientGroupsInputUI2"))

                          ),
                        
                        tags$br(),
                        conditionalPanel(condition="input.AnalysisDataType!=5",
                                         actionButton("goAnalysisButton", "Run Analysis"),
                                         helpText("Click the button to start analysis.")
                        )
                          ),
                      mainPanel(
                        width=9,
                        uiOutput("analysisTitle"),
                        h4("Patients you inputted"),
                        DT::dataTableOutput("inputtedPatientGroups"),
                        conditionalPanel(condition="input.AnalysisDataType==5",
                                         uiOutput("SurvivalStandardTerminologyUI1"),
                                         uiOutput("SurvivalStandardTerminologyUI2"),
                                         plotOutput("SurvivalPlot", height='100%', width='100%'),
                                         plotOutput("DFSurvivalPlot", height='100%', width='100%')),
                        conditionalPanel(condition="input.AnalysisDataType!=5",
                                         h4("Associated genes or clinical features"),
                                         DT::dataTableOutput("analysisResTable"),
                                         downloadButton("dowloadAnalysisRes", "Download full table of significant genes as .CSV file")
                        )
                        
                      )
                          )
                          ),
           tabPanel("News",
                    h3("News", style="color: STEELBLUE; padding-bottom: 20px"),
                    h4("May 05, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
                    tags$ul(
                      tags$li("Integrated image features."),
                      tags$li("Minor bugs fixed.")
                    ),
                    h4("April 18, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
                    tags$ul(
                      tags$li("Modifications finished.")
                    )
           ),
           tabPanel("Tutorial",
                    includeMarkdown("www/README.md")
           ),
           tabPanel("About",
                    h3("About Us", style="color: STEELBLUE; padding-bottom: 20px"),
                    "iGenomicsR",
                    tags$div(
                      tags$img(src='images/IUSM2.png',
                               height="100",
                               alt="IUSM", class="center", style="padding: 30px"),
                      tags$img(src='images/regenstrief.png',
                               height="100",
                               alt="Regenstrief", class="center", style="padding: 30px"),
                      style="text-align: center; padding: 20px"
                    ),
                    h4("Our Other Softwares", style="color: STEELBLUE; padding-bottom: 20px"),
                    tags$div(
                      a(tags$img(src='images/tsunami_logo.png',
                                 height="45",
                                 alt="TSUNAMI", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/tsunami/", target="_blank"),
                      br(),a("TSUNAMI: Translational Bioinformatics Tool SUite for Network Analysis and MIning",
                             href="https://apps.medgen.iupui.edu/rsc/tsunami/", target="_blank"),
                      br(),br(),
                      a(tags$img(src='images/lmQCM_logo.png',
                                 height="60",
                                 alt="lmQCM", class="center", style="padding: 5px"), href="https://CRAN.R-project.org/package=lmQCM", target="_blank"),
                      br(),a("R package: lmQCM", href="https://CRAN.R-project.org/package=lmQCM", target="_blank"),
                      br(),br(),
                      a(tags$img(src='images/annoPeak_logo.png',
                                 height="40",
                                 alt="annoPeak", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/content/19/", target="_blank"),
                      br(),a("annoPeakR: a web-tool to annotate, visualize and compare peak sets from ChIP-seq/ChIP-exo", href="https://apps.medgen.iupui.edu/rsc/content/19/", target="_blank"),
                      style="text-align: center; padding: 5px"
                    ),
                    br(),
                    tags$div(
                      a(tags$img(src='images/iGPSeplus_logo.png',
                                 height="40",
                                 alt="iGenomicsR", class="center", style="padding: 5px"),
                        href="https://apps.medgen.iupui.edu/rsc/content/23/", target="_blank"),
                      br(),a("iGPSe Plus: Integrative Genomic based Canser Patient Stratification", href="https://apps.medgen.iupui.edu/rsc/content/23/", target="_blank"),
                      style="text-align: center; padding: 5px"
                    ),
                    h4("Development Team", style="color: STEELBLUE; padding-bottom: 20px"),
                    h5("Prof. Kun Huang's Laboratory", style="color: STEELBLUE"),
                    h4("Publications", style="color: STEELBLUE; padding-bottom: 20px"),
                    tags$ul(
                      tags$li("-")
                    ),
                    h4("Funding for the iGPSe Plus is or has been provided by:", style="color: STEELBLUE; padding-bottom: 20px"),
                    tags$ul(
                      tags$li("Partially supported by IUSM startup fund, the NCI ITCR U01 (CA188547)."),
                      tags$li("Data Science and Bioinformatics Program for Precision Health Initiative, Indiana University.")
                    )
           ),
           tags$head(tags$script(HTML("document.title = 'iGenomicsR';"))), # rename the title by JS
           tags$div(
             p(a("iGenomicsR", href=""), "Version v1.2 | ", a("IUSM",href="https://medicine.iu.edu/", target="_blank"), " | ", a("RI",href="http://www.regenstrief.org/", target="_blank"), style="color: grey; font-size: 12px"), 
             p("Questions and feedback:  | ", a("Report Issue", href="https://github.com/huangzhii/iGenomicsR/issues", target="_blank"), " | ", a("Github", href="https://github.com/huangzhii/iGenomicsR/", target="_blank"), style="color: grey; font-size: 12px"),
             style="text-align: center; padding-top: 40px"
           )
)