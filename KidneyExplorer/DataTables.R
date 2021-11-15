#################### Code to print selected data #########################

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