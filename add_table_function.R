library(tidyverse)
library(RSQLite)
library(DBI)

add_table <- function(db_path, table_path, table_name, dec_char, separator, Species, ttest_difference_name, ttest_pvalue_name,
                      ttest_difference_name2 = "", ttest_pvalue_name2 = ""){
  ################################### Parameters : ######################################
  ## path to the database, path to the new table file, name of the new table,
  ## decimal character ("," or "."), separator character, 
  ## Species ('Human', 'Mouse' or 'Rat'), 
  ## name of the log2 fold change column, name of the -log10 p-value column
  ## If there are two tests in the table, then name of the second log2 fold-change column 
  ## name of the second -log10 p-value column
  #######################################################################################
  
  ## First, open the table in R
  # CSV file
  if (str_sub(table_path, str_length(table_path)-3, str_length(table_path)) == ".csv"){
  newTable <-read.csv(table_path, header = T, dec = dec_char, sep = separator) # dec = dec_char
  }
  # Text file
  else if (str_sub(table_path, str_length(table_path)-3, str_length(table_path)) == ".txt"){
    newTable <- read.table(table_path, header = T, sep=separator, dec = dec_char) # sep="\t"
  }
  # Other formats not accepted
  else {
    print ("The table must be in csv or txt format")
    return("")
  }
  # So that the user can have an idea of whether the table is read correctly
  print(head(newTable))
  
  ### Then, connect the database
  con <- dbConnect(RSQLite::SQLite(), db_path)
  
  ### If there are points in variable names, change it to underscores
  variables <- names(newTable)
  variables<- gsub("\\.","_",variables)
  names(newTable)<- variables
  
  ### Write the table in the database
  dbWriteTable(con,table_name, newTable)
  
  # If there are two tests in this table, then duplicate it in the database
  if (ttest_difference_name2 != ""){
  request <- paste("CREATE TABLE ", table_name, "2 AS 
          SELECT * FROM ", table_name, sep = "")
  dbExecute(con, request)
  }
  
  ### Rename the ttest-difference column so that the name is the same as in the other tables 
  request <- paste ("ALTER TABLE ", table_name,
                    " RENAME ", ttest_difference_name,
                    " TO ttest_difference1")
  dbExecute(con, request)
  
  # If you have created a second table, do that but we the other column
  if (ttest_difference_name2 != ""){
  ttest_diff_col2 <- "write the second name here"
  request <- paste ("ALTER TABLE ", table_name, "2",
                    " RENAME ", ttest_difference_name2,
                    " TO ttest_difference1", sep = "")
  dbExecute(con, request)
  }
  
  ### Rename the p-value column
  request <- paste ("ALTER TABLE ", table_name,
                    " RENAME ", ttest_pvalue_name,
                    " TO ttest_pvalue1")
  dbExecute(con, request)
  
  # If you have created a second table, do that but we the other column
  if (ttest_difference_name2 != ""){
    request <- paste ("ALTER TABLE ", table_name, "2",
                      " RENAME ", ttest_pvalue_name2,
                      " TO ttest_pvalue1", sep = "")
    dbExecute(con, request)
  }
  
  ### Make tidy dataset
  request <- paste("CREATE TABLE ", table_name, "_tidy AS
                   SELECT Gene_names, * FROM ", table_name, sep = "")
  dbExecute(con, request)
  dbExecute(con, paste("DELETE FROM ", table_name, "_tidy", sep = ""))
  
  data <- dbGetQuery(con, paste("SELECT * FROM ",table_name))
  for (i in 1:nrow(data)){
    gene_name <- strsplit(data$Gene_names[i],';')
    if (length(gene_name[[1]])!=0){
      for (j in 1:length(gene_name[[1]])){
        request <- paste("INSERT INTO ", table_name, "_tidy VALUES (", str_dup("?,", times=ncol(data)), "?)", sep = "")
        param <- c()
        for (h in 1:ncol(data)){
          param <- c(param, data[i,h])
        }
        dbExecute(con,request, params=c(gene_name[[1]][j],param))
      }
    }
  }
  
  # If a duplicated table was created : 
  if (ttest_difference_name2 != ""){
    request <- paste("CREATE TABLE ", table_name, "2_tidy AS
                   SELECT Gene_names, * FROM ", table_name, "2", sep = "")
    dbExecute(con, request)
    dbExecute(con, paste("DELETE FROM ", table_name, "2_tidy", sep = ""))
    
    data <- dbGetQuery(con, paste("SELECT * FROM ",table_name, "2", sep = ""))
    for (i in 1:nrow(data)){
      gene_name <- strsplit(data$Gene_names[i],';')
      if (length(gene_name[[1]])!=0){
        for (j in 1:length(gene_name[[1]])){
          request <- paste("INSERT INTO ", table_name, "2_tidy VALUES (", str_dup("?,", times=ncol(data)), "?)", sep = "")
          param <- c()
          for (h in 1:ncol(data)){
            param <- c(param, data[i,h])
          }
          dbExecute(con,request, params=c(gene_name[[1]][j],param))
        }
      }
    }
  }
  
  ### Modify the merge table
  mergeTable <-dbGetQuery(con, "SELECT * FROM mergeTable")
  # if there a duplicated table was created 
  if (ttest_difference_name2 != ""){
    request <- paste("SELECT matchSymbols.Symbol_Human AS Gene_names, ",
                     table_name, "_tidy.Gene_names AS Gene_names", table_name, ", ",
                     table_name, "2_tidy.Gene_names AS Gene_names", table_name, "2",
                     " FROM ", table_name, "_tidy 
      LEFT JOIN ", table_name, "2_tidy ON ", table_name, "_tidy.Gene_names = ", table_name, "2_tidy.Gene_names 
      LEFT JOIN matchSymbols ON ", table_name, "_tidy.Gene_names = matchSymbols.Symbol_", Species, sep = "")
    
  }
  else {
    request <- paste("SELECT matchSymbols.Symbol_Human AS Gene_names, ", table_name, "_tidy.Gene_names AS Gene_names", table_name,
                     " FROM ", table_name, "_tidy 
      LEFT JOIN matchSymbols ON ", table_name, "_tidy.Gene_names = matchSymbols.Symbol_", Species, sep = "")
  }
  newTable <- dbGetQuery(con, request )
  
  for (i in 1:nrow(newTable)){
    if (is.na(newTable[i,1])){
      newTable[i,1]<- newTable[i,2]
    }
  }
  
  data <- merge(mergeTable, newTable, all = T )
  dbExecute(con, "DROP TABLE mergeTable")
  dbWriteTable(con, "mergeTable", data)
  
  ### Table containing for each gene the number of table we can found it in
  dbExecute(con, 'DROP TABLE NbTableForGene')
  dbExecute(con, "CREATE TABLE NbTableForGene (
                HomoloGeneID INTEGER,
                Gene_names TEXT,
                Nb_tables INTEGER)")
  
  f <- function(n){
    for (i in 0:(n-1)){
      genes<- dbGetQuery(con, "SELECT *
           FROM mergeTable
           LIMIT 1 OFFSET ? ", params=i)
      nb<- ncol(genes)-5
      for (j in 6:ncol(genes)){
        if (is.na(genes[1,j])){
          nb <- nb-1
        }
      }
      test <- dbGetQuery(con,"SELECT * 
                       FROM NbTableForGene
                       Where Gene_names = ?", params = genes$Gene_names)
      if (nrow(test)>0){
        dbExecute(con, "DELETE FROM NbTableForGene
                      WHERE Gene_names = ? ", params = genes$Gene_names)
      }
      
      dbExecute(con, "INSERT INTO NbTableForGene VALUES (?,?, ?)", params = c(genes$HomoloGeneID, genes$Gene_names,nb))
    }
  }
  
  n<- dbGetQuery(con, "SELECT COUNT(*) FROM mergeTable")
  f(n[1,1])
  
  ### Delete the original table
  dbExecute(con, paste("DROP TABLE ", table_name))
  
  ### Rename the tidy table
  dbExecute(con, paste("ALTER TABLE ", table_name, "_tidy RENAME TO ", table_name, sep = ""))
  
  # Do the same if a duplicated table was created
  if (ttest_difference_name2 != ""){
    dbExecute(con, paste("DROP TABLE ", table_name, "2", sep = ""))
    dbExecute(con, paste("ALTER TABLE ", table_name, "2_tidy RENAME TO ", table_name, "2", sep = ""))
  }
  
  ### Close the connection
  dbDisconnect(con)
}

# Example of use of the function
# add_table("dbKidneyExplorer_final.db", "organoid_tnfa.csv", "hOrganoid_tnfa", ",", "", "Human", "log2_FCs_24h_TNFa_24h_VC", "neglog10_pvalue_24h_TNFa_24h_VC", "log2_FCs_48h_TNFa_48h_VC", "neglog10_pvalue_48h_TNFa_48h_VC")

# If one wants to rename one of the table, one can do : 
# dbExecute(con, "ALTER TABLE <current name> RENAME TO <new name>")

# If one needs to delete a table from the database 
# dbExecute(con, "DROP TABLE <table name>")

############################### Note ###################################
# If one encounter problem, one can try to delete the sep and dec arguments in read.csv 
# or read.table, so that the default parameters will be used

