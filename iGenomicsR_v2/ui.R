# Zhi Huang 04/02/2018
library(shinyWidgets)
navbarPage(title=div(a(img(src="images/iGenomicsR_logo2.png",
                           height = 28,
                           style = "margin:0px 0px; padding-bottom: 5px"), href=""),escape=F),
           tabPanel("Welcome",
                    sidebarLayout(
                      position = "left",
                      sidebarPanel(width = 3,
                        h4("File Uploader", style="color: STEELBLUE"),
                        fileInput("csvfile_mutation", "Mutation Table",
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
                        fileInput("csvfile_clinical", "Clinical Profile",
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
                          a(href =  'data/mRNA.sample.csv', 'mRNA.sample.csv'), ', ',
                          a(href = 'data/miRNA.sample.csv', 'miRNA.sample.csv'), ' and ',
                          a(href = 'data/time.cencer.csv', 'time.cencer.csv' ),
                          'files, and then try uploading them.'),
                        actionButton("action1", "Confirm when Complete")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("About",
                                   h2("iGenomicsR:", style="color: STEELBLUE; font-size: 22px"),
                                   h2("An integrative platform to explore, visualize, and analyze multidimensional genomics data for disease", style="color: STEELBLUE; font-size: 20px; margin: 0px"),
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
                          tabPanel("Dataset",
                                   h4("Data File required", style="color: STEELBLUE; padding-top: 10px"),
                                   
                                   fluidRow(
                                     column(3, "Mutation Profile", htmlOutput("check1")),
                                     column(3, "RNA Expression Profile", htmlOutput("check2")),
                                     column(3, "Protein Expression Profile", htmlOutput("check3")),
                                     column(3, "Clinical Profile", htmlOutput("check4")),
                                     style="text-align: center"
                                   ),
                                   tableOutput("dataUploadSummary"),
                                   conditionalPanel(condition="output.data_ready_flag == 1",
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
                                                   downloadButton("downloadSelectedMutation", "Download selected mutation data in .csv format"),br(),
                                                   downloadButton("downloadSelectedRNA", "Download selected RNA data in .csv format"),br(),
                                                   downloadButton("downloadSelectedProtein", "Download selected protein data in .csv format"),br(),
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
                                                         if(typeMessage == "tab3"){
                                                         $("a:contains(Survival Analysis)").click();
                                                         }
                                                         });
                                                         
                                                         // disable download at startup.
                                                         $(document).ready(function() {
                                                         $("#downloadSelectedMutation").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedRNA").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedProtein").attr("disabled", "true").attr("onclick", "return false;");
                                                         $("#downloadSelectedClin").attr("disabled", "true").attr("onclick", "return false;");
                                                         
                                                         Shiny.addCustomMessageHandler("download_seleted_data_ready", function(message) {
                                                         $("#downloadSelectedMutation").removeAttr("disabled").removeAttr("onclick");
                                                         $("#downloadSelectedRNA").removeAttr("disabled").removeAttr("onclick");
                                                         $("#downloadSelectedProtein").removeAttr("disabled").removeAttr("onclick");
                                                         $("#downloadSelectedClin").removeAttr("disabled").removeAttr("onclick");
                                                         });
                                                         })
                                                         '))
                                               ) # end of mainPanel
                          )
                        )
                    ) # end of sidebarLayout
                    
           ),
           tabPanel("Data Navigator"
                    
           ),
           tabPanel("Data Integration"
                    
           ),
           tabPanel("Data Analysis"
                    
           ),
           tabPanel("News"
                    
           ),
           tabPanel("FAQ"
                    
           ),
           tabPanel("About",
                    h3("About Us", style="color: STEELBLUE; padding-bottom: 20px"),
                    "iGenomicsR",
                    tags$div(
                      tags$img(src='images/IUSM2.png',
                               height="100",
                               alt="TSUNAMI", class="center", style="padding: 30px"),
                      tags$img(src='images/regenstrief.png',
                               height="100",
                               alt="TSUNAMI", class="center", style="padding: 30px"),
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
                      tags$img(src='images/iGPSeplus_logo.png',
                               height="40",
                               alt="iGenomicsR", class="center", style="padding: 5px"),
                      br(),"Coming Soon",
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
             p(a("iGenomicsR", href=""), "Version v1.1 | ", a("IUSM",href="https://medicine.iu.edu/", target="_blank"), " | ", a("RI",href="http://www.regenstrief.org/", target="_blank"), style="color: grey; font-size: 12px"), 
             p("Questions and feedback:  | ", a("Report Issue", href="https://github.com/huangzhii/iGenomicsR/issues", target="_blank"), " | ", a("Github", href="https://github.com/huangzhii/iGenomicsR/", target="_blank"), style="color: grey; font-size: 12px"),
             style="text-align: center; padding-top: 40px"
           )
)