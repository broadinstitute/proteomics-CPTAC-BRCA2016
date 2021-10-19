#################################################################
## Filename: global.r
## Created: February 10, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2 Breast
##          Cancer dataset (Mertins et al. 2016)
## This defines the server logic of the Shiny-app.
#################################################################
library(shiny)


########################################################
## Define server logic
########################################################
shinyServer( function(input, output, session) {

    global <- reactiveValues()

    ##############################
    ## update list of input genes
    observeEvent(input$genes, {
        global$genes.input <- extractGenes(input$genes)
    })

    ##############################
    ## generate the heatmap
    output$plot <- renderPlot({

        genes.vec <- extractGenes( input$genes )
        if(length(genes.vec)==0) return()

        hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, zscore=as.logical(input$zscore))
        global$expr.select <- hm
    },
    width = function(){ width=1200},
    height= function(){ height=dynamicHeightHM(length( global$genes.input ))}
    )

    #############################
    ## download heatmap
    output$downloadHM <- downloadHandler(

        ##filename = paste( FILENAMESTRING, ifelse(input$zscore, 'Zscore', ''),'.pdf', sep=''),
        filename = paste( FILENAMESTRING,  '-',  gsub(' |\\:','-', Sys.time()), '.pdf', sep=''),
        content = function(file){
            genes.vec <- extractGenes( input$genes )
            if(length(genes.vec)==0) return()
            hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, filename=file, main=TITLESTRING.HEATMAP, height=ifelse(length(genes.vec) < 4, 6.5, NA), as.logical(input$zscore))
            }
    )
    #############################
    ## download Excel
    output$downloadTab <- downloadHandler(
        filename = function(){paste( FILENAMESTRING, '-',  gsub(' |\\:','-', Sys.time()) ,'.xlsx', sep='') },
        content = function(file){
            tab=as.data.frame(global$expr.select)
            WriteXLS('tab', ExcelFileName=file, SheetNames=FILENAMESTRING, FreezeCol=6, FreezeRow=5, row.names=T)
            }
    )
})


