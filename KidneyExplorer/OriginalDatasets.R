############## Code to print and download original datasets ###################

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

# Download selected file
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

# Download all files in a zip
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