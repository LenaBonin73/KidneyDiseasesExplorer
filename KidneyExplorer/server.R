########################################################################
##################### Server Function ##################################
########################################################################
source("global.R")
server <- function(input, output, session) { 
    
    # update dataset checkbox selection in the sidebar
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
    
    ########################### Scatter plots #############################
    
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
    
    source("ScatterFoldChange.R", local = T)
    
    ############ Scatter plot of p-values ############
    
    source("ScatterPvalues.R", local = T)
    
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
                        prot <- paste("( '", list_prot[1], "'", sep = "")
                        if(length(list_prot)>1){
                            for (i in 2:length(list_prot)){
                                prot <- paste(prot, ", '", list_prot[i], "'", sep = "")
                            }
                        }
                        prot <- paste(prot, ")", sep = "")
                        request <- paste("SELECT Gene_names FROM mergeTable WHERE Gene_names IN ", prot, sep = "")
                        data_genes_final <- dbGetQuery(con, request)
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
    
    source("HeatmapPvalues.R", local = T)
    
    ########## Heatmap of fold-changes ##########
    
    source("HeatmapFoldChange.R", local = T)
    
    ############################### Messages ################################
    
    source("WarningMessages.R", local = T)
    
    ######## Description of selected datasets for the scatter plot ##########
    
    # Open the CSV file that contains descriptions
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
    
    source("DataTables.R", local = T)
    
    ########################## Original datasets ##########################
    
    source("OriginalDatasets.R", local = T)
    
} # Close server function