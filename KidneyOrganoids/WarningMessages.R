########################### Warning messages #########################

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