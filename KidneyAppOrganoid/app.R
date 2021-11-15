library(shiny)
library(shinydashboard)
library(shinyWidgets)
#library(shinyjs)
#library(rintrojs)
library(DT)
library(tidyverse)
library(RSQLite)
library(plotly)
library(heatmaply)
library(zip)
con <- dbConnect(RSQLite::SQLite(), "dbKidneyExplorer_final.db")

##########################################################################
##########################################################################
# This Shiny app has been developed using R version 3.6.3 and packages :
# heatmaply_1.1.1      viridis_0.5.1        viridisLite_0.3.0   
# plotly_4.9.2.1       RSQLite_2.2.1        forcats_0.5.0       
# stringr_1.4.0        dplyr_1.0.2          purrr_0.3.4         
# readr_1.4.0          tidyr_1.1.2          tibble_3.0.4        
# ggplot2_3.3.2        tidyverse_1.3.0      DT_0.16             
# rintrojs_0.2.2       shinyWidgets_0.5.4   shinydashboard_0.7.1
# shiny_1.5.0 
##########################################################################
##########################################################################


############################ Attention ###################################
# We need to have the list of table by species, this cannot be got directly as
# we can only get the name of all tables,
# So if you add a new dataset to the database, you need to add the name in the 
# list of the species it belongs to
###########################################################################

list_tables_human <- c('hPodocult_PAN_33', 
                       'hPodocult_PAN_37',
                       'hGlom_FSGS_1',
                       'hGlom_FSGS_2',
                       'hGlom_NephrinMutation_1',
                       'hGlom_NephrinMutation_2',
                       'hPUC_ACTN4Mutation',
                       'hPodocult_Actn4Mut_vsActn4WT',
                       'hPodocult_Actn4Mut_vsGFP',
                       'hOrganoid_d25_d21',
                       'hOrganoid_d29_d21',
                       'hOrganoid_24h_TNFa',
                       'hOrganoid_48h_TNFa')

list_tables_mouse <- c('mPodo_LPS',
                       'mPodo_Doxorubicin')

list_tables_rat <- c("rGlom_PAN_d2",
                     "rGlom_PAN_d4")

list_keywords <- sort(as.vector(dbGetQuery(con, "SELECT * FROM list_keywords")$keyword))

list_BiologicalProcess <- sort(as.vector(dbGetQuery(con, "SELECT * FROM list_BiologicalProcess")$BiologicalProcess))

list_CellularComponent <- sort(as.vector(dbGetQuery(con, "SELECT * FROM list_CellularComponent")$CellularComponent))

list_MolecularFunctions <- sort(as.vector(dbGetQuery(con, "SELECT * FROM list_MolecularFunction")$MolecularFunction))

##### Definition of a new button to download as SVG in addition to PNG ####
icon_svg_path = "M15.608,6.262h-2.338v0.935h2.338c0.516,0,0.934,0.418,0.934,0.935v8.879c0,0.517-0.418,0.935-0.934,0.935H4.392c-0.516,0-0.935-0.418-0.935-0.935V8.131c0-0.516,0.419-0.935,0.935-0.935h2.336V6.262H4.392c-1.032,0-1.869,0.837-1.869,1.869v8.879c0,1.031,0.837,1.869,1.869,1.869h11.216c1.031,0,1.869-0.838,1.869-1.869V8.131C17.478,7.099,16.64,6.262,15.608,6.262z M9.513,11.973c0.017,0.082,0.047,0.162,0.109,0.226c0.104,0.106,0.243,0.143,0.378,0.126c0.135,0.017,0.274-0.02,0.377-0.126c0.064-0.065,0.097-0.147,0.115-0.231l1.708-1.751c0.178-0.183,0.178-0.479,0-0.662c-0.178-0.182-0.467-0.182-0.645,0l-1.101,1.129V1.588c0-0.258-0.204-0.467-0.456-0.467c-0.252,0-0.456,0.209-0.456,0.467v9.094L8.443,9.553c-0.178-0.182-0.467-0.182-0.645,0c-0.178,0.184-0.178,0.479,0,0.662L9.513,11.973z"

dl_button <- list(
    name = "Download as SVG",
    icon = list(
        path = icon_svg_path,
        transform = "scale(0.84) translate(-1, -1)"
    ),
    click = htmlwidgets::JS("
          function (gd) {
        Plotly.downloadImage(gd, {
        filename: 'svg_download',
        format: 'svg',
        width: gd._fullLayout.width,
        height: gd._fullLayout.height
    })
  }
   "))

########################################################################
####################### UI Function ####################################
########################################################################

ui <- dashboardPage(
    dashboardHeader(title = "Kidney diseases explorer", 
                    titleWidth = 300,
                    tags$li(class = "dropdown",
                        
                            # Information button
                            dropMenu(
                                dropdownButton("Info", icon = icon('info')),
                                h3(strong('Description of datasets')),
                                p(source("about.R")),
                                placement = "left",
                                arrow = TRUE
                            )),
                    tags$li(
                        a(
                            strong("HELP"),
                            height = 40,
                            href = "https://github.com/LenaBonin73/KidneyDiseasesExplorer/blob/main/README.md",
                            title = "",
                            target = "_blank"
                        ),
                        class = "dropdown"
                    )),
    
    # Sidebar ------------------------------------------------------------
    dashboardSidebar(
        # Selection of proteins and datasets
        width = 300,
        tags$br(),
        submitButton("Click here to apply changes", width = "100%"),
        menuItem("Selection of proteins",
                 inputId = "idBigMenu",
                 tabName = "selection_proteins",
                 startExpanded = TRUE,
                 menuItem( "By Gene name",
                           tabName = "GeneName",
                            selectizeInput(
                               inputId = "nameProt",
                               label = "Select gene name(s)",
                               choices = sort(unique((dbGetQuery(con, "SELECT Gene_names
                                                         FROM mergeTable ")$Gene_names))),
                               selected = NULL,
                               multiple = TRUE
                           ),
                           textAreaInput(
                                      inputId = "GeneList",
                                      label = "Insert a list of gene names (separated by a space)",
                                      value = "",
                                      width = '100%',
                                      height = '1000%',
                                      placeholder = "NPHS1 NPHS2 PDGFRB PECAM1"
                                    )),
                 
                 menuItem("Proteins annotated with a certain keyword",
                          tabName = "protKeyword",
                          selectizeInput(
                              inputId = "keywordSelected",
                              label = "Protein's keyword",
                              choices = list_keywords,
                              selected = NULL,
                              multiple = TRUE
                          )
                 ),
                 menuItem("Proteins annotated with a certain Gene Ontology",
                          tabName = "ProtGO",
                          menuItem("Biological Process",
                                   tabName = "GO_BiologicalProcess",
                                   startExpanded = TRUE,
                                   selectizeInput(
                                       inputId = "BiologicalProcessSelected",
                                       label = "",
                                       choices = list_BiologicalProcess,
                                       selected = NULL,
                                       multiple = TRUE
                                   )
                          ),
                          menuItem("Cellular Component",
                                   tabName = "GO_CellularComponent",
                                   startExpanded = TRUE,
                                   selectizeInput(
                                       inputId = "cellularComponentSelected",
                                       label = "",
                                       choices = list_CellularComponent,
                                       selected = NULL,
                                       multiple = TRUE
                                   )
                          ),
                          menuItem("Molecular Function",
                                   tabName = "GO_MolecularFunction",
                                   startExpanded = TRUE,
                                   selectizeInput(
                                       inputId = "MolecularfunctionSelected",
                                       label = "",
                                       choices = list_MolecularFunctions,
                                       selected = NULL,
                                       multiple = TRUE
                                   )
                          )
                 ),
                 menuItem("Proteins must be present in at least ... datasets",
                          tabName = "protNumber",
                          sliderInput(
                              inputId = "sliderNbProt",
                              label = "Number of datasets",
                              min = 1L,
                              max = length(list_tables_human)+length(list_tables_mouse)+length(list_tables_rat),
                              value = 1L,
                              step = 1L
                          )
                 )),
        
        menuItem("Confidence level for t-test",
                 tabName = "ttest",
                 selectizeInput(
                   inputId = "confidence",
                   label = "",
                   choices = c("Display all", "95%", "99%"),
                   selected = "Display all",
                   multiple = FALSE
                 )
                 ),
        
        menuItem("Show only one copy per gene ?",
                 tabName = "copy",
                 selectizeInput(
                   inputId = "showCopy",
                   label = "",
                   choices = c("No", "Yes"),
                   selected = "No",
                   multiple = FALSE
                 )
        ),
        
        menuItem("Selection of dataset for the heatmap",
                 tabName = "Dataselect",
                 checkboxInput('bar', 'All/None', value = TRUE),
                 submitButton("Apply"),
                 checkboxGroupInput(
                     inputId = "dataInputHuman",
                     label = "Human data",
                     choices = list_tables_human,
                     selected = list_tables_human
                 ),
                 checkboxGroupInput(
                     inputId = "dataInputMouse",
                     label = "Mouse data",
                     choices = list_tables_mouse,
                     selected = list_tables_mouse
                 ),
                 checkboxGroupInput(
                     inputId = "dataInputRat",
                     label = "Rat data",
                     choices = list_tables_rat,
                     selected = list_tables_rat
                 ))
    ),
    
    # Body ---------------------------------------------------------------
    dashboardBody(
        
        fluidRow(
            div(id = "scatter_panel",
                column(12,
                       tabBox(
                           title = "",
                           id = "tabset2", height = 480,
                           tabPanel("log2(Fold-changes)", textOutput("text_scatter_fold"),
                                    tabPanel("Fold-change2",plotlyOutput("plot_ab"))),
                           tabPanel("-log10(P-values)", textOutput("text_scatter_p"),
                                    tabPanel("P-values2",plotlyOutput("pvalue_scatter")))
                       ),
                       
                       box(
                           title = "Choose parameters for the scatter plot", solidHeader = TRUE,
                           selectInput("abscissa", "Abscissa (x-axis)", choices = c(list_tables_human, list_tables_mouse, list_tables_rat), selected = c(list_tables_human, list_tables_mouse, list_tables_rat)[1]),
                           textOutput('descriptionAbscissa'),
                           br(),
                           selectInput("ordinates", "Ordinates (y-axis)", choices = c(list_tables_human, list_tables_mouse, list_tables_rat), selected = c(list_tables_human, list_tables_mouse, list_tables_rat)[2]),
                           textOutput('descriptionOrdinates'),
                           br(),
                           selectInput("density", "Display density", choices = c("Yes","No"), selected = "No" ),
                           selectInput("genes1", "Proteins to display", choices = c("All (color the selected ones)", "Selected"), selected = "All (color the selected ones)"),
                           selectInput("annots", "Proteins to label", choices = c("None", "Selected", "Only selected by gene name"), selected = "Only selected by gene name"), 
                           submitButton("Display")
                       )))
        ),
        
        fluidRow(
            div(
                id = "heatmap_panel",
                tabBox(
                    title = "",
                    id = "tabset1", width = 12, #,height=500,
                    tabPanel("log2(Fold-changes)", 
                             textOutput("text_fold_heatmap"),
                             tabPanel("log2(Fold-changes)", plotlyOutput("heatmap", height = "600%")
                             )),
                    tabPanel("-log10(P-values)", textOutput("text_pvalue_heatmap"),
                             tabPanel("heat_graph",plotlyOutput("pvalue_heatmap", height = "600%")))
                )
            )
        ),
        
        fluidRow(
            div(
                id = "DataTables_panel",
                tabBox(
                    title = "Data",
                    id = "tabset3", width = 12,
                    tabPanel("Scatter Fold-change",DT::dataTableOutput("datatable_SF")),
                    tabPanel("scatter P-value", DT::dataTableOutput("datatable_SP")),
                    tabPanel("Heatmap Fold-change", DT::dataTableOutput("datatable_HF")),
                    tabPanel("Heatmap P-value", DT::dataTableOutput("datatable_HP"))
                )
            )
        ),
        fluidRow(
          div(
            id = "OriginalDataDiv",
            tabBox(
              title = "Original datasets",
              id = "OriginalData",
              width = 12,
                tabPanel(
                  title = "Download Orginal Datasets",
                    selectizeInput(inputId = "datasetDownload",
                            label = "",
                            choices = c(list_tables_human, list_tables_mouse, list_tables_rat),
                            selected = c(list_tables_human, list_tables_mouse, list_tables_rat)[1],
                            multiple = F),
                  submitButton("Click to apply changes"),
                  tags$br(),
                  downloadButton("downloadData", "Download as CSV"),
                  downloadButton('downloadAll', "Download all datasets"),
                  dataTableOutput("firstRowsOriginal")
                    ),
              tabPanel("Description of datasets", DT::dataTableOutput("dataset_desc"))
                )
            )
        )
    )
)

##########################################################################
########################### Server function ##############################
##########################################################################

server <- function(input, output, session) { 
  
    observe({
      updateCheckboxGroupInput(
        session, "dataInputHuman", choices = list_tables_human,
        selected = if (input$bar) list_tables_human
      )
      updateCheckboxGroupInput(
        session, "dataInputMouse", choices = list_tables_mouse,
        selected = if (input$bar) list_tables_mouse
      )
      updateCheckboxGroupInput(
        session, "dataInputRat", choices = list_tables_rat,
        selected = if (input$bar) list_tables_rat
      )
    })
  
    # Create a R list if a list of genes is provided in textAreaInput
    # The list contains gene names
    ListOfGenes <- reactive({
      if (str_length(input$GeneList)>0){
        gene_names <- strsplit(input$GeneList, " ")
        gene_names <- gene_names[[1]]
      }
      else{
        gene_names <- c()
      }
      gene_names
    })
  
    ########################### Scatter plot #############################
    
    ##### List of selected proteins (using "selection of proteins")
    data_scatter <- reactive({
        ### Genes under study ###
        list_prot <- input$nameProt
        keywords <- input$keywordSelected
        BP <- input$BiologicalProcessSelected
        CC <- input$cellularComponentSelected
        MF <- input$MolecularfunctionSelected
        nbDataset <- input$sliderNbProt
        
        # First we look at keywords and GO
        if(length(keywords)+length(BP)+length(CC)+length(MF)+nbDataset>1){
            genes_studied <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
            if (length(keywords)+length(BP)+length(CC)+length(MF)>0){
                
                request <- "SELECT Gene_names FROM Keyword_annotations_tidy WHERE"
                
                WHERE <- ""
                if (length(keywords>0)){ # keywords
                    for (i in 1:length(keywords)){
                        WHERE <- paste(WHERE, "Keywords LIKE '%", keywords[i], "%' AND ", sep = "" )
                    }
                }
                
                if (length(BP)>0){ # GO Biological Process
                    for (i in 1:length(BP)){
                        WHERE <- paste(WHERE, "GO_BiologicalProcess LIKE '%", BP[i], "%' AND ", sep = "" )
                    }
                }
                if (length(CC)>0){ # GO Cellular Component
                    for (i in 1:length(CC)){
                        WHERE <- paste(WHERE, "GO_CellularComponent LIKE '%", CC[i], "%' AND ", sep = "" )
                    }
                }
                if (length(MF)>0){ # GO Molecular Function
                    for (i in 1:length(MF)){
                        WHERE <- paste(WHERE, "GO_molecular_function LIKE '%", MF[i], "%' AND ", sep = "" )
                    }
                }
                
                WHERE <- substr(WHERE, 1, nchar(WHERE)-5)
                
                request <- paste (request, WHERE)
                request
                
                genes_studied <- dbGetQuery(con, request) # We get a dataset
            }
            
            # Second we look at how many datasets proteins must be in => a second dataset
            all_genes_nb <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
            if (nbDataset >1){
                request <- "SELECT Gene_names FROM NbTableForGene WHERE Nb_tables >= ?"
                all_genes_nb <- distinct(dbGetQuery(con, request, params = nbDataset))
            }
            
            # we merge the 2 dataset in order to keep only the genes that meet both requirements
            if(nrow(all_genes_nb)>0 & nrow(genes_studied)>0){
                data_genes <- distinct(merge(genes_studied, all_genes_nb))
            }
            else if (nrow(genes_studied)>0){
                data_genes <- genes_studied
            }
            else if (nrow(all_genes_nb)>0) {
                data_genes <- all_genes_nb
            }
            else {data_genes <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))}
            
            list_studied_genes <- as.vector(data_genes$Gene_names)
        }
        # If no selection criterion is selected, we arbitrary set list_studied_genes to "()"
        else{list_studied_genes <- "()"}
        
        # If no selection is selected except by Gene names : 
        # We want the gene names to be considered as a selection
        if (list_studied_genes == "()"){
          if(length(ListOfGenes())>0){
          list_studied_genes <- ListOfGenes()
          }
          else if(length(input$nameProt)>0){
            list_studied_genes <- input$nameProt
          }
        }
        #else if (list_studied_genes == "()" & length(input$nameProt)>0){
         # list_studied_genes <- input$nameProt
        #}
        
        list_studied_genes
    })
    
    #########Scatter plot of fold changes ###########
    
    ###### dataset for the fold-change scatter plot, scatter text and fold-change scatter datatable
    data_scatter_fold <- reactive({
        list_studied_genes <- data_scatter()
        
        # selection of data 
        matrix1 <- input$abscissa
        matrix2 <- input$ordinates
        
        # If a confidence level is selected
        if (input$confidence == "95%"){
          where <- "WHERE ttest_pvalue1 >= -log10(0.05)"
        }
        else if (input$confidence == "99%"){
          where <- "WHERE ttest_pvalue1 >= -log10(0.01)"
        }
        else {
          where <-""
        }
        
        # We need to match the gene_names for both matrix if at least one is not a human dataset
        if (matrix1 %in% list_tables_mouse){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix1,
                             ".ttest_difference1 AS ttest_matrix1 FROM mergeTable JOIN ",
                             matrix1, " ON mergeTable.Symbol_Mouse = ", matrix1, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix1, ".Gene_names", sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix1 %in% list_tables_rat){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix1,
                             ".ttest_difference1 AS ttest_matrix1 FROM mergeTable JOIN ",
                             matrix1, " ON mergeTable.Symbol_Rat = ", matrix1, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix1, ".Gene_names", sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix1 %in% list_tables_human){
            request <- paste("SELECT ", matrix1, ".Gene_names AS Gene_names, ", matrix1,".ttest_difference1 AS ttest_matrix1 FROM ",matrix1, sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        
        # idem for matrix2
        if (matrix2 %in% list_tables_mouse){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix2,
                             ".ttest_difference1 AS ttest_matrix2 FROM mergeTable JOIN ",
                             matrix2, " ON mergeTable.Symbol_Mouse = ", matrix2, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix2, ".Gene_names", sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix2 %in% list_tables_rat){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix2,
                             ".ttest_difference1 AS ttest_matrix2 FROM mergeTable JOIN ",
                             matrix2, " ON mergeTable.Symbol_Rat = ", matrix2, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix2, ".Gene_names", sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix2 %in% list_tables_human){
            request <- paste("SELECT ", matrix2, ".Gene_names AS Gene_names, ", matrix2,".ttest_difference1 AS ttest_matrix2",
                             " FROM ", matrix2, sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        
        # End of matching
        # Creation of the dataset containing both statistics
        data <- merge(data_matrix1, data_matrix2)
        data <- distinct(data, .keep_all = T)
        
       
        # If we want to display only selected proteins
        if (input$genes1 == "Selected"){
          data <- filter(data, Gene_names %in% list_studied_genes)
        }
        
        # If we want to display only one line per gene
        if(input$showCopy == "Yes"){
          doublons <- which(duplicated(data[,1], fromLast=TRUE))
          if(length(doublons)>0){
            data <- data[-doublons,]
          }
        }
        
        # If the variables are considered as factor, transform them to numeric
        if(is.numeric(data$ttest_matrix1)==F){
          data$ttest_matrix1 <- as.numeric(gsub(",",".",data$ttest_matrix1,fixed=TRUE))
        }
        if(is.numeric(data$ttest_matrix2)==F){
          data$ttest_matrix2 <- as.numeric(gsub(",",".",data$ttest_matrix2,fixed=TRUE))
        }
  
        # Return data
        data
    })
    
    ##### Display the scatter plot
    output$plot_ab <- renderPlotly({
        data<- data_scatter_fold()
        list_studied_genes <- data_scatter()
        
        if(nrow(data)>0){ 
            # Because if there is no selected protein there is nothing to display
            
            #Genes selected by gene name
            genelist <- ListOfGenes()
            if(length(genelist)>0){
              selected_genes <- genelist
            }
            else if (length(input$nameProt)>0){
              selected_genes <- c()
              for (i in 1:length(input$nameProt)){
                selected_genes<- c(selected_genes, input$nameProt[i])
              }
            }
            else {selected_genes<- c()}
            
            # All selected genes
            all_selected <- c(selected_genes, data_scatter())
            all_selected <- unique(all_selected)
            
            # Defining colors
            if(input$genes1=="Selected"){
              data$color <- "#0370EC"
            }
            else{
              data$color<- ifelse(data$Gene_names%in%all_selected, "red", "#0370EC")
            }
            
            # Defining the scatter plot t015EC7  004FA8 
            p <- ggplot(data, aes(x=ttest_matrix1, y=ttest_matrix2))+
              geom_point(aes(text =Gene_names), color = data$color)+
              #geom_text(aes(label = ifelse(Gene_names %in% selected_genes, Gene_names, ""), fontface = "bold", size = 5),color = "black")+
              xlab(paste("Log2(fold change)", input$abscissa))+
              ylab(paste("Log2(fold change)", input$ordinates))+
              NULL
            
            # if density = yes, we plot the density
            if(input$density == "Yes"){
              p <- p + geom_density2d(color = "navyblue")
            }
            
            fig <- ggplotly(p, tooltip = c("text")) %>% 
              config(
                modeBarButtonsToRemove = list("hoverCompareCartesian", "hoverClosestCartesian", "lasso2d", "select2d", "autoScale2d"),
                modeBarButtonsToAdd = list(dl_button)
              )
            # Annotations :
            if (input$annots!="None"){
              if (input$annots == "Selected"){
                m<- filter(data, Gene_names %in% all_selected)
              }
              else{
                m<- filter(data, Gene_names %in% selected_genes)
              }
            if (nrow(m)>0){
              a <- list(
              x = m$ttest_matrix1,
              y = m$ttest_matrix2,
              text = m$Gene_names,
              xref = "x",
              yref = "y",
              showarrow = TRUE,
              arrowhead = 5,
              arrowsize = 0.5,
              ax = 10,
              ay = -30,
              font = list(color = 'black',
                          family = 'sans serif',
                          size = 14)
            )
            fig <- fig %>% add_markers()
            fig <- fig %>% layout(annotations = a)
            }
            }
            # returning the plot
            fig
        }
    })

    
    ############ Scatter plot of p-values ############
    
    ##### Dataset for the p-value scatter plot and p-value scatter datatable
    data_scatter_pvalue <- reactive({
        list_studied_genes <- data_scatter()
        # selection of data 
        matrix1 <- input$abscissa
        matrix2 <- input$ordinates
        
        # If a confidence level is selected
        if (input$confidence == "95%"){
          where <- "WHERE ttest_pvalue1 >= -log10(0.05)"
        }
        else if (input$confidence == "99%"){
          where <- "WHERE ttest_pvalue1 >= -log10(0.01)"
        }
        else {
          where <-""
        }
        
        #We need to match the gene_names for both matrix if at least one is not a human dataset
        if (matrix1 %in% list_tables_mouse){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix1,
                             ".ttest_pvalue1 AS ttest_matrix1 FROM mergeTable JOIN ",
                             matrix1, " ON mergeTable.Symbol_Mouse = ", matrix1, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix1, ".Gene_names", sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix1 %in% list_tables_rat){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix1,
                             ".ttest_pvalue1 AS ttest_matrix1 FROM mergeTable JOIN ",
                             matrix1, " ON mergeTable.Symbol_Rat = ", matrix1, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix1, ".Gene_names", sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix1 %in% list_tables_human){
            request <- paste("SELECT ", matrix1, ".Gene_names AS Gene_names, ", matrix1,".ttest_pvalue1 AS ttest_matrix1 FROM ",matrix1, sep = "")
            data_matrix1 <- dbGetQuery(con, paste(request, where))
        }
        # idem for matrix2
        if (matrix2 %in% list_tables_mouse){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix2,
                             ".ttest_pvalue1 AS ttest_matrix2 FROM mergeTable JOIN ",
                             matrix2, " ON mergeTable.Symbol_Mouse = ", matrix2, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix2, ".Gene_names", sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix2 %in% list_tables_rat){
            request <- paste("SELECT mergeTable.Gene_names AS Gene_names, ", matrix2,
                             ".ttest_pvalue1 AS ttest_matrix2 FROM mergeTable JOIN ",
                             matrix2, " ON mergeTable.Symbol_Rat = ", matrix2, ".Gene_names",
                             " OR mergeTable.Gene_names = ", matrix2, ".Gene_names", sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        if (matrix2 %in% list_tables_human){
            request <- paste("SELECT ", matrix2, ".Gene_names AS Gene_names, ", matrix2,".ttest_pvalue1 AS ttest_matrix2",
                             " FROM ", matrix2, sep = "")
            data_matrix2 <- dbGetQuery(con, paste(request, where))
        }
        # End of matching
        # creation of a dataset containing both statistics
        data <- merge(data_matrix1, data_matrix2)
        data <- distinct(data, .keep_all = T)
        
        # If we have a list of genes to study (selected by GO, keywords or nb of datasets)
        #if (list_studied_genes != "()"){
        #    data <- filter(data, Gene_names %in% list_studied_genes)
        #}
        # If we want to display only selected proteins
        if (input$genes1 == "Selected"){
          data <- filter(data, Gene_names %in% list_studied_genes)
        }
        
        # If we want to display only one line per gene
        if(input$showCopy == "Yes"){
          doublons <- which(duplicated(data[,1], fromLast=TRUE))
          if(length(doublons)>0){
            data <- data[-doublons,]
          }
        }
        
        # If the variables are considered as factor, transform them to numeric
        if(is.numeric(data$ttest_matrix1)==F){
          data$ttest_matrix1 <- as.numeric(gsub(",",".",data$ttest_matrix1,fixed=TRUE))
        }
        if(is.numeric(data$ttest_matrix2)==F){
          data$ttest_matrix2 <- as.numeric(gsub(",",".",data$ttest_matrix2,fixed=TRUE))
        }
        
        # return the dataset
        data
    })
    
    ##### Display the scatter plot
    output$pvalue_scatter <- renderPlotly({
      data<- data_scatter_pvalue()
      list_studied_genes <- data_scatter()
        
        if(nrow(data)>0){ 
          # Because if there is no selected protein there is nothing to display
          
          #Genes selected by gene name
          genelist <- ListOfGenes()
          if(length(genelist)>0){
            selected_genes <- genelist
          }
          else if (length(input$nameProt)>0){
            selected_genes <- c()
            for (i in 1:length(input$nameProt)){
              selected_genes<- c(selected_genes, input$nameProt[i])
            }
          }
          else {selected_genes<- c()}
          
          # All selected genes
          all_selected <- c(selected_genes, data_scatter())
          all_selected <- unique(all_selected)
          
          # Defining colors
          if(input$genes1=="Selected"){
            data$color <- "#0370EC"
          }
          else{
            data$color<- ifelse(data$Gene_names%in%all_selected, "red", "#0370EC")
          }
          
          # Defining the scatter plot t015EC7  004FA8 
          p <- ggplot(data, aes(x=ttest_matrix1, y=ttest_matrix2))+
                geom_point(aes(text =Gene_names), color = data$color)+
                #geom_text(aes(label = ifelse(Gene_names %in% selected_genes, Gene_names, ""), fontface = "bold", size = 5),color = "black")+
                xlab(input$abscissa)+
                ylab(input$ordinates)+
                NULL
            
            # if density = yes, we plot the density
            if(input$density == "Yes"){
                p <- p + geom_density2d(color = "navyblue")
            }
            
            # Creating the plot
            fig <- ggplotly(p, tooltip = c("text"))%>% 
                config(
                  modeBarButtonsToRemove = list("hoverCompareCartesian", "hoverClosestCartesian", "lasso2d", "select2d", "autoScale2d"),
                  modeBarButtonsToAdd = list(dl_button)
                )
            # Annotations :
            if (input$annots!="None"){
              if (input$annots == "Selected"){
                m<- filter(data, Gene_names %in% all_selected)
              }
              else{
                m<- filter(data, Gene_names %in% selected_genes)
              }
              if (nrow(m)>0){
                a <- list(
                  x = m$ttest_matrix1,
                  y = m$ttest_matrix2,
                  text = m$Gene_names,
                  xref = "x",
                  yref = "y",
                  showarrow = TRUE,
                  arrowhead = 5,
                  arrowsize = 0.5,
                  ax = 10,
                  ay = -30,
                  font = list(color = 'black',
                              family = 'sans serif',
                              size = 14)
                )
                fig <- fig %>% add_markers()
                fig <- fig %>% layout(annotations = a)
              }
            }
            # returning the plot
            fig
        }
    })

    
    ############################## Heatmap #################################
    
    ##### Protein selection 
    ##### Creation of a SQL list containing proteins selected by gene name, keyword, GO, dataset number
    data_heatmap <- reactive({
        # We do not need to do that is no dataset is selected because then there is no heatmap to display
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            # Inputs
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
            list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # If a selection criterion is selected (excepted selection by gene name)
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
                # Creation of an empty dataset in case this criterion is not selected (to make us able to use the merge function later) 
                genes_studied <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
                
                # First we look at keywords and GO 
                if (length(keywords)+length(BP)+length(CC)+length(MF)>0){
                    
                    request <- "SELECT DISTINCT Gene_names FROM Keyword_annotations_tidy WHERE"
                    
                    WHERE <- ""
                    if (length(keywords>0)){ # Keywords
                        for (i in 1:length(keywords)){
                            WHERE <- paste(WHERE, "Keywords LIKE '%", keywords[i], "%' AND ", sep = "" )
                        }
                    }
                    
                    if (length(BP)>0){ # GO Biological Process
                        for (i in 1:length(BP)){
                            WHERE <- paste(WHERE, "GO_BiologicalProcess LIKE '%", BP[i], "%' AND ", sep = "" )
                        }
                    }
                    if (length(CC)>0){ # GO Cellular COmponent
                        for (i in 1:length(CC)){
                            WHERE <- paste(WHERE, "GO_CellularComponent LIKE '%", CC[i], "%' AND ", sep = "" )
                        }
                    }
                    if (length(MF)>0){ # GO Molecular Function
                        for (i in 1:length(MF)){
                            WHERE <- paste(WHERE, "GO_molecular_function LIKE '%", MF[i], "%' AND ", sep = "" )
                        }
                    }
                    
                    WHERE <- substr(WHERE, 1, nchar(WHERE)-5)
                    
                    request <- paste (request, WHERE)
                    
                    # We get a dataset containing proteins under study according to these criterions
                    # (It replace the empty dataset we previoulsy created)
                    genes_studied <- dbGetQuery(con, request)
                }
                
                # Second : if protein must be present in more than 1 dataset
                
                # Creation of an empty dataset in case this criterion is not selected (to make us able to use the merge function later)
                all_genes_nb <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
                if (nbDataset >1){
                    request <- "SELECT Gene_names FROM NbTableForGene WHERE Nb_tables >= ?"
                    all_genes_nb <- distinct(dbGetQuery(con, request, params = nbDataset))
                }
                
                # Third : we merge the 2 datasets in order to keep only the genes that meet both requirements
                if(nrow(all_genes_nb)>0 & nrow(genes_studied)>0){
                    data_genes <- distinct(merge(genes_studied, all_genes_nb))
                }
                else if (nrow(genes_studied)>0){
                    data_genes <- genes_studied
                }
                else if (nrow(all_genes_nb)>0) {
                    data_genes <- all_genes_nb
                }
                else {data_genes <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))}
                
                # Fourth :  If genes are selected by their name
                if(length(list_prot)>=1){
                    if (nrow(data_genes)>0){
                        # some proteins are selected by previous criterions
                        vector_prot <- c()
                        for (i in 1:length(list_prot)){
                            vector_prot <- c(vector_prot, list_prot[i])
                        }
                        data_genes_final <- filter(data_genes, Gene_names %in% vector_prot)
                    }
                    else{
                        # They is no protein selected by previous criterions 
                        prot <- paste("('", list_prot[1], "'", sep = "")
                        if(length(list_prot)>1){
                            for (i in 2:length(list_prot)){
                                prot <- paste(prot, ",'", list_prot[i], "'", sep = "")
                            }
                        }
                        prot <- paste(prot, ")", sep = "")
                        request <- paste("SELECT Gene_names FROM mergeTable WHERE Gene_names IN ", prot, sep = "")
                        data_genes_final <- distinct(dbGetQuery(con, request))
                    }
                } # end of if "nameProt"
                else {
                    # No gene is selected by its name
                    data_genes_final <- data_genes
                }
                
                # If there are more than 3000 genes to study, we do not try to look for proteins in the datasets
                if(nrow(data_genes_final)<=2000){
                # We have a dataset with all the genes under study, we need to make a SQL list
                list_studied_genes <- "("
                for (i in 1:nrow(data_genes_final)){
                    list_studied_genes <- paste(list_studied_genes, "'", data_genes_final$Gene_names[i], "', ", sep = "")
                }
                list_studied_genes <- substr(list_studied_genes, 1, nchar(list_studied_genes)-2)
                list_studied_genes <- paste(list_studied_genes, ")", sep = "")
                }
                else{
                  list_studied_genes <- "too many"
                }
          } # close if at least 1 selection criterion
        } # close if at least 1 dataset
      
            # If there is no criterion selected we arbitraly set list_studied_genes to "()"
            else{list_studied_genes <- "()"} 
        
            # At this point we have the list of the genes we have to display
            list_studied_genes
    })
    
    
    ############ Heatmap of p-values ############
    
    ##### Creation of the dataset containing Gene names and selected ttest p-values
    data_pvalue_heat <- reactive({
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
              list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # Step 1 : finding proteins selected via keywords, GO, gene name or dataset numbers (what we did in data_heatmap)
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
                list_studied_genes <- data_heatmap()
                
                #if there are less than 1000 proteins
                if (list_studied_genes != "too many"){
                
                # Step 2 : finding the t-stat for the genes we selected
                
                # We have to treat each species separately because for non-human species 
                #we need to match the selected gene names with the gene names of that species
                datasetMouse <- input$dataInputMouse
                datasetRat <- input$dataInputRat
                datasetHuman <- input$dataInputHuman
                
                # Confidence level 
                if(input$confidence == "95%"){
                  conf <- "-log10(0.05)"
                }
                else if(input$confidence == "99%"){
                  conf <- "-log10(0.01)"
                }
                else {
                  conf <- "-10000000000"
                }
                
                
                # Retrieving values
                selection <- "SELECT mergeTable.Gene_names"
                from <- "FROM mergeTable"
                if (list_studied_genes != "()"){
                  where <- paste("WHERE mergeTable.Gene_names IN ", list_studied_genes, " AND ", sep = "")
                }
                else {
                  where <- "WHERE "
                }
                data <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
                
                if (length(datasetMouse)>0){
                    for (i in 1:length(datasetMouse)){
                        selectionI<- paste(selection, ", ", datasetMouse[i],'.ttest_pvalue1 AS MinusLogPvalue_', datasetMouse[i], sep = "")
                        fromI <- paste(from, " LEFT JOIN ", datasetMouse[i], ' ON mergeTable.Symbol_Mouse = ', datasetMouse[i], ".Gene_names OR mergeTable.Gene_names = ", datasetMouse[i], ".Gene_names", sep="")
                        whereI <- paste (where, datasetMouse[i], '.ttest_pvalue1 >= ', conf, sep = "")
                        dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
                        data <- merge(data, dataI, all = T)
                    }}
                if(length(datasetRat)>0){
                    for (i in 1:length(datasetRat)){
                        selectionI<- paste(selection, ", ", datasetRat[i],'.ttest_pvalue1 AS MinusLogPvalue_', datasetRat[i], sep = "")
                        fromI <- paste(from, " LEFT JOIN ", datasetRat[i], ' ON mergeTable.Symbol_Rat = ', datasetRat[i], ".Gene_names OR mergeTable.Gene_names = ", datasetRat[i], ".Gene_names", sep="")
                        whereI <- paste (where, datasetRat[i], '.ttest_pvalue1 >= ', conf, sep = "")
                        dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
                        data <- merge(data, dataI, all = T)
                    }}
                if (length(datasetHuman)>0){
                    for (i in 1:length(datasetHuman)){
                        selectionI<- paste(selection, ", ", datasetHuman[i],'.ttest_pvalue1 AS MinusLogPvalue_',datasetHuman[i], sep = "")
                        fromI <- paste(from, " LEFT JOIN ", datasetHuman[i], ' ON mergeTable.Gene_names = ', datasetHuman[i], ".Gene_names", sep="")
                        whereI <- paste (where, datasetHuman[i], '.ttest_pvalue1 >= ', conf, sep = "")
                        dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
                        data <- merge(data, dataI, all = T)
                    }}
                
                # If we want to display only one line per gene
                if(input$showCopy == "Yes"){
                  doublons <- which(duplicated(data[,1], fromLast=TRUE))
                  if(length(doublons)>0){
                    data <- data[-doublons,]
                  }
                }
                
                # Returning the data table under study
                data
                }
            }
        }
        
    })
    
    ##### Displaying the heatmap of p-values
    output$pvalue_heatmap<- renderPlotly({
        # If datasets are selected we potentially can display a heatmap
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            # Inputs
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
            list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # If selection criterions are selected we potentially can display a heatmap (if there are more than one protein that meet the requirements)
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
              # If there are more than 1000 selected proteins
              if(data_heatmap()!="too many"){
                # Getting the dataset of selected genes and ttest p-values
                test <- data_pvalue_heat()
            
                # Function to treat the problem of NA in heatmapply clustering
                dist_no_na <- function(mat1) {
                    mat <- mat1
                    for (i in 1:nrow(mat)){
                        for (j in 1:ncol(mat)){
                            if (is.na(mat[i,j])){
                                mat[i,j]<- 100 # we arbitraly decide to assign 100 because it is much higher than every existed number
                            }
                        }
                    }
                    edist <- dist(mat)
                    return(edist)
                }
                
                # If there are between 2 and 1000 selected protein, we can display a heatmap
                if(nrow(test)<=1000 & nrow(test)>1){
                    # To display the heatmap, we need at least one value of p-value (not only NA)
                    if (sum(is.na(test[,-1]))!=nrow(test)*ncol(test[,-1])){
                        
                        # Treatment of the duplicates
                      if(TRUE %in% duplicated(test$Gene_names)){
                        test[1,1]<-paste(test[1,1], "_1", sep = "")
                        i<-2
                        j<-1
                        while(i<=nrow(test)){
                            if (test[i,1]==substr(test[i-j,1], 1, nchar(test[i-j,1])-2)){
                                j<-j+1
                                test[i,1]<- paste(test[i,1],"_", j, sep = "")
                                i<-i+1
                            }
                            else {
                                j<-1
                                test[i,1]<- paste(test[i,1], '_1', sep = "")
                                i<-i+1
                            }
                        }
                      }
                        
                        col_names <- c()
                        for (i in 2:ncol(test)){
                          cn <- colnames(test)[i]
                          col_names <- c(col_names, substr(cn, 16, nchar(cn)))
                        }
                        
                        # heatmap
                        heatmaply(test[,-1],
                                  distfun = dist_no_na,
                                  dendrogram = "row",
                                  ylab = "",
                                  xlab = "",
                                  main = "",
                                  #scale = "column",
                                  margins = c(30,100,0,0),
                                  #grid_color = "white",
                                  grid_width = 0.00001,
                                  #limits = c(-2,2),
                                  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                                      low = "blue", 
                                      high = "red", 
                                      mid = 'grey92',
                                      midpoint = 1.3,
                                      limits = c(min(test[,-1], na.rm = TRUE) , max(test[,-1], na.rm = TRUE) )
                                  ),
                                  titleX = FALSE,
                                  hide_colorbar = FALSE,
                                  branches_lwd = 0.1,
                                  labRow = test[,1],
                                  labCol = col_names,
                                  heatmap_layers = theme(axis.line=element_blank()) )%>% 
                            config(
                                modeBarButtonsToAdd = list(dl_button)
                            )
                      }
                    }
                }
            }
        } # close if "there is at least one dataset selected
    })
    
    
    ############ Heatmap of fold-changes ############ 
    
    ##### Creation of the dataset containing Gene names and selected ttest difference values
    ##### This dataset is also used to know if we have a print a message (if it has one value or zero) and for the datatable
    data_fold_heat <- reactive({
        # To display a heatmap we need to have selected datasets
      if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
        if(length(ListOfGenes())>0){
          list_prot <- ListOfGenes()
        }
        else{
          list_prot <- input$nameProt}
        keywords <- input$keywordSelected
        BP <- input$BiologicalProcessSelected
        CC <- input$cellularComponentSelected
        MF <- input$MolecularfunctionSelected
        nbDataset <- input$sliderNbProt
        
        # Step 1 : finding proteins selected via keywords, GO, gene name or dataset numbers (what we did in data_heatmap)
        if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
          list_studied_genes <- data_heatmap()
          
          if (list_studied_genes!="too many"){
          
          # Step 2 : finding the t-stat for the genes we selected
          
          # We have to treat each species separately because for non-human species 
          #we need to match the selected gene names with the gene names of that species
          datasetMouse <- input$dataInputMouse
          datasetRat <- input$dataInputRat
          datasetHuman <- input$dataInputHuman
          
          # Confidence level 
          if(input$confidence == "95%"){
            conf <- "-log10(0.05)"
          }
          else if(input$confidence == "99%"){
            conf <- "-log10(0.01)"
          }
          else {
            conf <- "-10000000000"
          }
          
          
          # Retrieving values
          selection <- "SELECT mergeTable.Gene_names"
          from <- "FROM mergeTable"
          if (list_studied_genes != "()"){
            where <- paste("WHERE mergeTable.Gene_names IN ", list_studied_genes, " AND ", sep = "")
          }
          else {
            where <- "WHERE "
          }
          data <- as.data.frame(set_names(replicate(1,numeric(0), simplify = F), c('Gene_names')))
          
          if (length(datasetMouse)>0){
            for (i in 1:length(datasetMouse)){
              selectionI<- paste(selection, ", ", datasetMouse[i],'.ttest_difference1 AS ttest_difference_', datasetMouse[i], sep = "")
              fromI <- paste(from, " LEFT JOIN ", datasetMouse[i], ' ON mergeTable.Symbol_Mouse = ', datasetMouse[i], ".Gene_names OR mergeTable.Gene_names = ", datasetMouse[i], ".Gene_names", sep="")
              whereI <- paste (where, datasetMouse[i], '.ttest_pvalue1 >= ', conf, sep = "")
              dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
              data <- merge(data, dataI, all = T)
            }}
  
          if(length(datasetRat)>0){
            for (i in 1:length(datasetRat)){
              selectionI<- paste(selection, ", ", datasetRat[i],'.ttest_difference1 AS ttest_difference_', datasetRat[i], sep = "")
              fromI <- paste(from, " LEFT JOIN ", datasetRat[i], ' ON mergeTable.Symbol_Rat = ', datasetRat[i], ".Gene_names OR mergeTable.Gene_names = ", datasetRat[i], ".Gene_names", sep="")
              whereI <- paste (where, datasetRat[i], '.ttest_pvalue1 >= ', conf, sep = "")
              dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
              data <- merge(data, dataI, all = T)
            }}
          
          if (length(datasetHuman)>0){
            for (i in 1:length(datasetHuman)){
              selectionI<- paste(selection, ", ", datasetHuman[i],'.ttest_difference1 AS ttest_difference_',datasetHuman[i], sep = "")
              fromI <- paste(from, " LEFT JOIN ", datasetHuman[i], ' ON mergeTable.Gene_names = ', datasetHuman[i], ".Gene_names", sep="")
              whereI <- paste (where, datasetHuman[i], '.ttest_pvalue1 >= ', conf, sep = "")
              dataI <- distinct(dbGetQuery(con, paste(selectionI, fromI, whereI)))
              data <- merge(data, dataI, all = T)
            }}
          
          # If we want to display only one line per gene
          if(input$showCopy == "Yes"){
            doublons <- which(duplicated(data[,1], fromLast=TRUE))
            if(length(doublons)>0){
              data <- data[-doublons,]
            }
          }
          
          # Returning the data table under study
          data
          } # close else
        } # close a selection criterion is selected
      } # close at least a dataset is selected 
    })
    
    
    ### Displaying heatmap of fold-changes
    output$heatmap<- renderPlotly({
        # If datasets are selected, we can potentially display a heatmap
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            # Inputs
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
              list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # If selection criterions are selected we potentially can display a heatmap (if there are more than one protein that meet the requirements)
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
              
              if (data_heatmap()!="too many"){
                # Getting the dataset of selected genes and fold-changes
                test <- data_fold_heat()
        
        
        #Function to treat the problem of NA in heatmapply clustering
        dist_no_na <- function(mat1) {
            mat <- mat1
            for (i in 1:nrow(mat)){
                for (j in 1:ncol(mat)){
                    if (is.na(mat[i,j])){
                        mat[i,j]<- 100 # we arbitraly decide to assign 100 because it is much higher than every existed number
                    }
                }
            }
            edist <- dist(mat)
            return(edist)
        }
        
        # If there is more than 1 selected protein, we can potentially display a heatmap
        if(nrow(test)>1 & nrow(test)<=200){
            # To display the heatmap, we need at least one value of fold-change (not only NA)
        if (sum(is.na(test[,-1]))!=nrow(test)*ncol(test[,-1])){
            
            # Treatment of the duplicates
          if(TRUE %in% duplicated(test$Gene_names)){
            test[1,1]<-paste(test[1,1], "_1", sep = "")
            i<-2
            j<-1
            while(i<=nrow(test)){
                if (test[i,1]==substr(test[i-j,1], 1, nchar(test[i-j,1])-2)){
                    j<-j+1
                    test[i,1]<- paste(test[i,1],"_", j, sep = "")
                    i<-i+1
                }
                else {
                    j<-1
                    test[i,1]<- paste(test[i,1], '_1', sep = "")
                    i<-i+1
                }
            }
          }
            
            col_names <- c()
            for (i in 2:ncol(test)){
              cn <- colnames(test)[i]
              col_names <- c(col_names, substr(cn, 18, nchar(cn)))
            }
            
            # heatmap
            heatmaply(test[,-1],
                      distfun = dist_no_na,
                      dendrogram = "row",
                      xlab = "", ylab = "", 
                      main = "",
                      #scale = "column",
                      margins = c(30,100,0,0),
                      #grid_color = "white",
                      grid_width = 0.00001,
                      #limits = c(-2,2),
                      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                          low = "blue", 
                          high = "red", 
                          mid = 'grey92', #ivory2
                          midpoint = 0,
                          limits = c(min(test[,-1], na.rm = TRUE) , max(test[,-1], na.rm = TRUE) )
                      ),
                      titleX = FALSE,
                      hide_colorbar = FALSE,
                      branches_lwd = 0.1,
                      labRow = test[,1],
                      labCol = col_names,
                      heatmap_layers = theme(axis.line=element_blank())
                      ) %>% 
                config(
                    modeBarButtonsToAdd = list(dl_button)
                )
                  } # close if not all na
                } # close if nrow(test)>1
              } # close if less than 1500 proteins
            } # close if "at least a selection criterion"
       } # close if "there is at least one dataset selected
    })

    ############################### Messages ################################
    
    ##### Message for the heatmap 
    output$text_fold_heatmap <- output$text_pvalue_heatmap <- renderText({
        # Inputs
        if(length(ListOfGenes())>0){
          list_prot <- ListOfGenes()
        }
        else{
          list_prot <- input$nameProt}
        keywords <- input$keywordSelected
        BP <- input$BiologicalProcessSelected
        CC <- input$cellularComponentSelected
        MF <- input$MolecularfunctionSelected
        nbDataset <- input$sliderNbProt
        text_to_print <- ""
        
        # If criterions are selected but less than 1 proteins meet the requirements 
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
              if(data_heatmap()!="too many"){
                test <- data_fold_heat()
                # Zero protein
                if(nrow(test)<1){
                    text_to_print <- "There is no data which meets all the selected requirements"
                }
                # 1 protein
                else if(nrow(test)==1){
                    text_to_print <- paste("There is only one data which meets all the selected requirements. It is the one associated with the gene name",test[1]$Gene_names)
                }
                # More than 200 proteins
                else if (nrow(test)>200){
                  text_to_print <- "There are more than 200 data which meet the selected requirements in the datasets, the heatmap cannot be displayed.\n
                  You may consider ploting only one copy per gene name (if not already done)"
                }
                # Otherwise we display a heatmap, so there is no message to print
                else {text_to_print <- ""}
              }
              else{
                text_to_print <- "There are more than 2000 selected proteins to be looked for:\n
                You may consider adding filters to the selection"
              }
            }
        }
        # If no dataset is selected or no criterion is selected, we cannot display a heatmap and so a message is printed
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))< 1
           | (length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset)<=1){
            text_to_print <- "Select proteins and datasets to display their heatmap"
        }
        text_to_print
    })
    
    ##### Message for the scatter plot 
    output$text_scatter_fold <- output$text_scatter_p <-  renderText({
        text_to_print <- ""
        if(nrow(data_scatter_fold())==0){
            text_to_print <- "There is no data which meets all the selected requirements"
        }
        text_to_print
    })
    
    ##### Description of selected datasets for the scatter plot
    
    # Open the dataset that contains description
    descr_table <- reactive({
      data <- read.csv("dataset_descriptions.csv", sep = "\t")
      # Transformation of the interesting variables into characters (instead of factors)
      data$short_name <- as.character(data$short_name) 
      data$short_description <- as.character(data$short_description)
      data
    })
    
    # Description of the selected dataset for abscissa
    output$descriptionAbscissa <- renderText({
      data <- descr_table()
      data2 <- filter(data, short_name == input$abscissa)
      data2[1,3]
    })
    
    # Description of the selected dataset for ordinates
    output$descriptionOrdinates <- renderText({
      data <- descr_table()
      data2 <- filter(data, short_name == input$ordinates)
      paste(data2[1,3], "\n")
    })
    
    
    ###################### Table with selected data #######################
    
    ##### Data table scatter plot fold-change
    output$datatable_SF <- DT::renderDataTable({
        data <- data_scatter_fold()
        colnames(data) <- c("Gene_names", paste("Fold-change", input$abscissa), paste("Fold-change", input$ordinates))
        datatable(
            data,
            rownames = FALSE,
            extensions = "Buttons",
            options = list(
                lengthMenu = c(10, 25, 50, 100),
                #dom = 'lBSfrtpi',
                dom = '<"row"<"col-sm-4"l><"col-sm-4 text-center"<B>><"col-sm-4"f>>tipS',
                buttons = c('csv', 'excel', 'pdf'),
                style = "bootstrap"
            )
        ) %>% 
            formatRound(columns = names(data)[-1], digits=3)
    }, server = FALSE)
    
    ##### Data table scatter plot p-value
    output$datatable_SP <- DT::renderDataTable({
        data <- data_scatter_pvalue()
        colnames(data) <- c("Gene_names", paste("P-value", input$abscissa), paste("P-value", input$ordinates))
        datatable(
            data,
            rownames = FALSE,
            extensions = "Buttons",
            options = list(
              dom = '<"row"<"col-sm-4"l><"col-sm-4 text-center"<B>><"col-sm-4"f>>tipS',
                buttons = c('csv', 'excel', 'pdf'),
                style = "bootstrap",
                lengthMenu = c(10, 25, 50, 100)
            )
        ) %>% 
            formatRound(columns = names(data)[-1], digits=3)
    }, server = FALSE)
    
    ##### Data table heatmap fold-change
    output$datatable_HF <- DT::renderDataTable({
        
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
            list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # Getting the dataset of selected genes and fold-changes
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
                data <- data_fold_heat()
        
            datatable(
                data,
                rownames = FALSE,
                extensions = "Buttons",
                options = list(
                    dom = '<"row"<"col-sm-4"l><"col-sm-4 text-center"<B>><"col-sm-4"f>>tipS',
                    buttons = c('csv', 'excel', 'pdf'),
                    style = "bootstrap",
                    lengthMenu = c(10, 25, 50, 100),
                    scrollX = TRUE
                )) %>% 
                formatRound(columns = names(data)[-1], digits=3)
    }}}, server = FALSE)
    
    ##### Data table heatmap p-value
    output$datatable_HP <- DT::renderDataTable({
        
        if((length(input$dataInputHuman)+length(input$dataInputMouse)+length(input$dataInputRat))>1){
            if(length(ListOfGenes())>0){
              list_prot <- ListOfGenes()
            }
            else{
              list_prot <- input$nameProt}
            keywords <- input$keywordSelected
            BP <- input$BiologicalProcessSelected
            CC <- input$cellularComponentSelected
            MF <- input$MolecularfunctionSelected
            nbDataset <- input$sliderNbProt
            
            # Getting the dataset of selected genes and p-values
            if(length(keywords)+length(BP)+length(CC)+length(MF)+length(list_prot)+nbDataset>1){
                data <- data_pvalue_heat()
                
                datatable(
                    data,
                    rownames = FALSE,
                    extensions = "Buttons",
                    options = list(
                      dom = '<"row"<"col-sm-4"l><"col-sm-4 text-center"<B>><"col-sm-4"f>>tipS',
                        buttons = c('csv', 'excel', 'pdf'),
                        style = "bootstrap",
                        lengthMenu = c(10, 25, 50, 100),
                        scrollX = TRUE
                    )
                ) %>% 
                    formatRound(columns = names(data)[-1], digits=3)
            }}}, server = FALSE)
    
    
    output$dataset_desc <- DT::renderDataTable({
        data <- read.csv("dataset_descriptions.csv", sep = "\t")
        datatable(
            data,
            rownames = FALSE,
            extensions = "Buttons",
            options = list(
                dom = '<"row"<"col-sm-4"l><"col-sm-4 text-center"<B>><"col-sm-4"f>>tipS',
                pageLength = 13,
                buttons = c('csv', 'excel', 'pdf'),
                style = "bootstrap",
                scrollX = TRUE
            ),
        )
    }, server = FALSE)
    
    ##### Original datasets
    # Datatable
    output$firstRowsOriginal <- renderDataTable({
      request <- paste("SELECT HomoloGeneID,", 
                       input$datasetDownload, ".Gene_names,",
                       input$datasetDownload, ".ttest_difference1,",
                       input$datasetDownload, ".ttest_pvalue1",
                       " FROM mergeTable, ", input$datasetDownload,
                       " WHERE mergeTable.Gene_names = ", input$datasetDownload,".Gene_names", 
                       " LIMIT 10", sep = "") 
      datatable(
      dbGetQuery(con, request),
      caption = "First ten rows of the dataset",
      options = list(
        dom = '<"top">rt',
        scrollX=TRUE, 
        scrollCollapse=TRUE))
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$datasetDownload, ".csv", sep = "")
      },
      content = function(file) {
        request <- paste("SELECT * FROM ", input$datasetDownload)
        data <- dbGetQuery(con, request)
        write.csv(data, file, row.names = FALSE)
      }
    )
    
    output$downloadAll <- downloadHandler(
      filename = function(){
        "origninal_datasets.zip"
      },
      content = function(file){
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- c()
        
        AllData <- c(list_tables_human, list_tables_mouse, list_tables_rat)
        #loop through the sheets
        for (i in 1:length(AllData)){
          #write each sheet to a csv file, save the name
          fileName <- paste(AllData[i],".csv",sep = "")
          request <- paste("SELECT * FROM ", AllData[i])
          data <- dbGetQuery(con, request)
          write.csv(data, fileName, row.names = FALSE)
          files <- c(files, fileName)
        }
        #create the zip file
        zip(file,files)
      }
    )
   
} # Close server function

#####################################################
# Launching the application 
#####################################################
options(servr.port = 4322)
shinyApp(ui = ui, server = server, options = "port")