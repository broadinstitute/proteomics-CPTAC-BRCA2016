#################################################################
## Filename: ui.r
## Created: February 10, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2 Breast
##          Cancer dataset (Mertins et al. 2016)
## This file defines the user interface of the Shiny-app.
#################################################################
library(shiny)

#########################################################
## Define UI
#########################################################
shinyUI(fluidPage(

  # Application title
  titlePanel(HTML(TITLESTRING), windowTitle = "Mertins et al. 2016. Nature"),

  ##################################
  ## side bar with:
  ## - text input field
  ## - Submit button
  ## - download buttons for pdf/xlsx
  sidebarLayout(

      sidebarPanel(
          HTML('<br>'),
          ## text input
          textInput('genes', label=paste('Paste a list of gene names (max. ', GENEMAX,')', sep=''), value=GENESSTART),
          ## submit button
          fluidRow(
              column(6, submitButton('GO')),
              column(6, checkboxInput('zscore', 'Apply (row) Z-score'), value=FALSE)
              ##column(6, HTML('<a href=\"help.html\" target=\"blank\">Help</a>'))
          ),
          HTML('<br><br>'),

          ## download buttons
          fluidRow(
            column(6, downloadButton('downloadHM', 'Download pdf')),
            column(6, downloadButton('downloadTab', 'Download Excel'))
          ),

          HTML('<br><br>'),

          HTML('<p><b>Getting started</b></p>'),
          helpText('Simply enter or paste your gene names of interest (official gene symbols, e.g. ERBB2) into the text field. The text field accepts lists of up to 20 gene symbols in either comma-, semicolon-, or space-separated form. The dataset provides quantitative data on 16,826 genes, however, not every data type will be available for every gene. If enabled Z-scoring will be applied on expression data and not on (discrete) CNA data.'),
          HTML('<p>For more details please see our publication <a href="https://www.nature.com/articles/nature18003" target="_blank_">Mertins et al. Nature. 2016</a></p>')   
           ),
    ################################
    ## main panel: heatmap
    mainPanel(
      plotOutput("plot")
    )
   # HTML('<p align="right"><a href="https://www.broadinstitute.org/proteomics" target="_blank_"><b>Proteomics Platform@Broad</b></a></p>')
  )
))
