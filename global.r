#################################################################
## Filename: global.r
## Created: February 10, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2 Breast
##          Cancer dataset (Mertins et al. 2016)
## This file imports the underlying data and contains global functions
## and variables used by 'ui.R' and 'server.R'
#################################################################

source('pheatmap.r')
library(scales)
library(gtable)

## import the data
load('data_bc2016.RData')

## global parameters
GENESSTART <<- 'TP53 ERBB2 PIK3CA GATA3 ESR1 PGR'
GENEMAX <<- 20
#TITLESTRING <<- 'Supplementary information: <a href="https://www.nature.com/articles/nature18003" target-"_blank_">Mertins et al. 2016</a>'
#TITLESTRING <<- 'Supplementary information to <i>Proteogenomics connects somatic mutations to signalling in breast cancer.</i><a href="https://www.nature.com/articles/nature18003" target-"_blank_"> (Mertins <i>et al</i>. 2016)</a>'
#TITLESTRING <<- '<a href="https://www.nature.com/articles/nature18003" target="_blank_">Proteogenomics connects somatic mutations to signalling in breast cancer</a>. Mertins <i>et al.</i> (2016). Nature.<br><br>Supplementary information'
TITLESTRING <<- '<font size="5" face="times"><i><b>"Proteogenomics connects somatic mutations to signalling in breast cancer"</b></i> (<a href="https://www.nature.com/articles/nature18003" target="_blank_">Mertins <i>et al.</i> Nature. 2016</a>)</font><br>'


FILENAMESTRING <<- 'CPTAC2_BC2016'

##library(pheatmap)
library(RColorBrewer)
library(gplots)
library(WriteXLS)
library(grid)

##################################################################
## function to extract gene names from a string of
## character
extractGenes <- function(genes.char){

    gene.max=GENEMAX

    if( nchar(genes.char) == 0 ){
        return(NULL)
    }
    ## extract genes
    genes.vec= unlist(strsplit(genes.char, ','))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ' '))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ';'))

    ## unique gene names
    genes.vec <- unique(genes.vec)

    ## limit to 'gene.max' genes
    if(length(genes.vec) > gene.max){
        warning(paste('more than', gene.max,'gene ids submitted! Showing results for the first 20 genes in the list.\n'))
        genes.vec <- genes.vec[1:gene.max]
    }
    return(genes.vec)
}

##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n){
    if( n < 5)
        height=500
    if( n >= 5 & n < 10)
        height=500
    if( n >= 10)
        height=500+n*20
    return(height)
}

#################################################################
## draw the actual heatmap
##
#################################################################
makeHM <- function(gene, filename=NA, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, zscore=F, ...){

    min.val=-3
    max.val=3
    n.bins=12

    ## remove spaces
    gene <- gsub(' ', '', gene )

    ## reverse the order
    gene <- gene[length(gene):1]

    ## make unique
    gene <- unique(gene)

    ## check whether the genes are present in the dataset
    gene.idx <-  grep( paste(paste('^', gene,'$', sep=''), collapse='|'), row.anno[, 'Gene name'] )
    if( length(gene.idx) == 0 ){
        stop('None of the gene ids you have entered could be found in the dataset!\n')
    }
    #################################
    ## extract sample ids
    sampleIDs <- colnames(expr)

    #################################
    ## extract genes of interest
    expr.select <- expr[gene.idx, ]
    row.anno.select <- row.anno[gene.idx, ]

    #####################################################
    ## labels for the rows in the heatmap
    featureIDs.anno.select <- paste(row.anno.select[ , 'Gene name'],  gsub('6-RPPA_pSTY', 'RPPA pSTY',  gsub('5-RPPA_prot', 'RPPA Protein',  gsub( '4_pSTY', 'MS pSTY', gsub('3_Protein', 'MS Protein', gsub('2_RNAseq', 'RNA-Seq', gsub('1_CNA', 'CNA', row.anno.select[ , 'Data type'])) )))))

    #####################################################
    ## add phosphosite annotation
    #####################################################
    ## RPPA
    rppa.psty.idx <- grep('RPPA pSTY', featureIDs.anno.select)
    if(length(rppa.psty.idx)>0){
        featureIDs.anno.select[ rppa.psty.idx ] <- paste( sub(' pSTY', '', featureIDs.anno.select[ rppa.psty.idx ]), sub('.*_(.*?)\\..*', '\\1', row.anno.select[ rppa.psty.idx, 'ID2']))
    }
    ## MS
    ms.psty.idx <- grep('MS pSTY', featureIDs.anno.select)
    if(length(ms.psty.idx)>0){
        featureIDs.anno.select[ ms.psty.idx ] <- paste( sub(' pSTY', '', featureIDs.anno.select[ ms.psty.idx ]), paste('p', sub('.*_([S|T|Y][0-9]*)[s|t|y].*', '\\1', row.anno.select[ ms.psty.idx, 'ID2']), sep='') )
    }

    #################################
    ## apply zscore
    rownames(expr.select) <- featureIDs.anno.select

    if(zscore){
        ## exclude CNA data from Z-scoreing
        expr.select.zscore.tmp <- lapply( rownames(expr.select), function(xx){x=expr.select[xx,];
            if( length(grep( 'CNA', xx)) == 0)return((x-mean(x, na.rm=T))/sd(x, na.rm=T));
            if( length(grep( 'CNA', xx)) > 0)return(x);
        })
        expr.select.zscore <- matrix(unlist(expr.select.zscore.tmp), ncol=ncol(expr.select), byrow=T, dimnames=list(rownames(expr.select), colnames(expr.select)))
    } else {
        expr.select.zscore <- expr.select
    }

    ## cap at -3/3
    expr.select.zscore[which(expr.select.zscore < min.val)] <- min.val
    expr.select.zscore[which(expr.select.zscore > max.val)] <- max.val

    ##############################
    ## column annotation
    column.anno.fig <- data.frame(
        PAM50=as.character(column.anno[c('PAM50'), ]),
        ER.Status=as.character(column.anno[c('ER.Status'), ]),
        PR.Status=as.character(column.anno[c('PR.Status'), ]),
        HER2.Status=as.character(column.anno[c('HER2.Status'), ])                          )
    rownames(column.anno.fig) <- colnames(column.anno)
    column.anno.fig <- column.anno.fig[, c('HER2.Status', 'PR.Status', 'ER.Status', 'PAM50')]

    ##############################
    ## colors for column annotation
    column.anno.col <- list(
        PAM50=c(Basal='red', Her2='violet', LumA='blue', LumB='cyan', Normal='grey'),
        ER.Status=c(Positive='black', Negative='white', Normal='grey'),
        PR.Status=c(Positive='black', Negative='white', Normal='grey'),
        HER2.Status=c(Positive='black', Negative='white', Normal='grey',  Equivocal='grey90')
    )

    ################################
    ## gaps
    ################################
    ## only possible because matrix is ordered according to PAM50
    gaps.column=cumsum(c(  table(column.anno.fig$PAM50) ))
    gaps.row=cumsum(table(sub(' .*', '', featureIDs.anno.select)))


    ################################
    ## colors misc

    color.breaks = seq(min.val, max.val, length.out=n.bins)
    color.hm =  colorRampPalette( c('blue', 'grey', 'red'))(length(color.breaks))
    color.border = 'white'

    ##legend_labels=c('-3 | CNA deletion', '-2', '-1 | CNA LOH', ' 0 | CNA neutral', '+1 | CNA gain', '+2' ,'+3 | CNA amplification')
    legend_breaks=seq(-3, 3, 1)
    legend_labels=c('-3               ', '-2', '-1', ' 0', '+1', '+2' ,'+3')


    ###############################
    ## heatmap
    cellwidth=10
    cellheight=10
    pheatmap(expr.select.zscore, cluster_row=F, cluster_col=F,  annotation_col=column.anno.fig, annotation_colors=column.anno.col,  scale = "none", labels_row=featureIDs.anno.select, border_color=color.border, gaps_col=gaps.column, gaps_row=gaps.row, color=color.hm, filename=filename, cellwidth=cellwidth, cellheight=cellheight, labels_col=sampleIDs, breaks=color.breaks, legend_breaks=legend_breaks, legend_labels=legend_labels, na_col='white',...)


    #########################################################################################
    ##
    ## - return part of table that is shown in the heatmap and that can be downloaded
    ## - change the CNA values (-3, -1, 0, 1, 3) back to the orignial values (-1, -.3, 0, .3, 1)
    ##
    #########################################################################################
    cna.idx <- grep('CNA', rownames(expr.select))
    if(length(cna.idx) > 0){
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -3 ] <- -1
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -1 ] <- -.3
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 1 ] <- .3
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 3 ] <- 1
    }
    ## add row annotation
    mat.row <- as.data.frame( cbind( row.anno[gene.idx, ], expr.select, deparse.level=0 ), stringsAsFactors=F )
    rownames(mat.row) <- featureIDs.anno.select

    ## column annotation
    mat.col <- as.data.frame( cbind(matrix('', nrow=ncol(column.anno.fig), ncol=ncol(row.anno)), t(column.anno.fig), deparse.level=0), stringsAsFactors=F)
    colnames(mat.col) <- colnames(mat.row)

    ## put everything together
    mat <- rbind( mat.col, mat.row)

    return(mat)
}
