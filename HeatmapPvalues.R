###################### Code for p-value heatmap #####################

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
