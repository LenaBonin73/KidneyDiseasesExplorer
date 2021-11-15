########################################################################
####################### UI Function ####################################
########################################################################
source("global.R")

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

