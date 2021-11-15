############# Code for fold-change scatter plot #####################

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