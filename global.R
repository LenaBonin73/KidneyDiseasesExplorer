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
                       'hPodocult_Actn4Mut_vsGFP')

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
