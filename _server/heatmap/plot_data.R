
#----------------------------------------------------------------------
# plot_data.R controls the data plotting actions
#----------------------------------------------------------------------

# define the heatmap colors
getColorPalette <- function(){
    reportProgress('getColorPalette')
    nColorsPerSide <- as.integer(input[['colorsPerSide']])
    nColors <- nColorsPerSide * 2 + 1
    cols <- colorRampPalette(brewer.pal(n=11, name=input[['rcbPalette']]))(nColors)
    if(input[['invertPalette']]) cols <- rev(cols)
    list(cols=cols, nColorsPerSide=nColorsPerSide, nColors=nColors)
}
dataToColors <- function(d=NULL, usePos=TRUE){
    
    # get configuration
    maxD <- as.numeric(input[['zmax']])    
    paletteGradient <- input[['paletteGradient']]
    colorPalette <- getColorPalette()    
    if(is.null(d)) d <- matrix(seq(-maxD, maxD, length=colorPalette$nColors), ncol=1)

    # scale the values for color selection
    if(paletteGradient == 'linear'){
        fracMax <- d / maxD
    } else if(paletteGradient == 'squared'){
        fracMax <- sqrt(ifelse(d < 0, -d, d)) / sqrt(maxD)
        fracMax <- ifelse(d < 0, -fracMax, fracMax)
    } else if(paletteGradient == 'sqrt'){
        fracMax <- d**2 / maxD**2
        fracMax <- ifelse(d < 0, -fracMax, fracMax)
    }
    fracMax <- ifelse(d > 0, pmin(1,fracMax), pmax(-1,fracMax))
    
    # assign colors per grid rectangle
    colorI <- round(colorPalette$nColorsPerSide * fracMax, 0) + colorPalette$nColorsPerSide + 1
    cols <- matrix(colorPalette$cols[colorI], nrow=nrow(colorI))
    cols[!usePos,] <- 'grey'
    list(cols=cols, d=d)
}

# plot variables
replicates <- list()
allPos  <- integer()
allPosK <- integer()
nAllPos <- 0
timepoints  <- character()
nTimepoints <- 0

# generic function to plot a custom heat map
plotHeatMap <- function(d, main, tall, file){
    reportProgress('plotHeatMap')
    
    # set the output colors
    cols <- dataToColors(d, allPos %notin% plotMask$dsb)
    
    # convert the X axis to # bp of resection, i.e. distance from DSB
    # negative values are on the OTHER side of the DSB, only decectable when there is no DSB
    resPos <- allPos - targetData[[input$dsbEnd]]$dsb$dsb
    
    # make the heat map plot
    posLab <- "resection position (bp)"
    tpLab  <- "time (min)"
    posRange <- range(resPos)
    tpRange  <- c(1,nTimepoints+1)
    if(tall){
        if(!is.null(file)){
            png(filename = file,
                width = as.numeric(input$saveWidth), height = as.numeric(input$saveHeight),
                pointsize = as.numeric(input$saveFontSize),
                res = 600, units = "in", type = "cairo")            
            par(mar=c(4.1,4.1,2.1,0.1))
        }
        plot(0, 0, typ="n",
             xlim=tpRange, ylim=posRange, xaxt="n",
             xlab=tpLab,   ylab=posLab,   yaxt="n",
             xaxs="i", yaxs="i")
        ybot  <- rep(resPos, nTimepoints) - 0.5
        xleft <- as.vector(sapply(1:nTimepoints, function(y) rep(y, nAllPos)))    
        rect(1, posRange[1], nTimepoints+1, posRange[2], col="grey", border=NA)
        rect(xleft, ybot, xleft + 1, ybot + 1, col=cols$cols, border=NA)
        tpLab <- sub("T", "", timepoints)
        axis(side=1, at=1:nTimepoints + 0.5, labels=tpLab, las=2, tick=FALSE, cex=7/8)
        axis(side=2, at=seq(-1000,1000,25), cex=7/8)
        mtext(main, side=3, line=0.5, font=2)
        if(!is.null(file)){
            graphics.off()
        } 
    } else {
        plot(0, 0, typ="n", main=main,
             xlim=posRange, ylim=tpRange,
             xlab=posLab,   ylab=tpLab, yaxt="n",
             xaxs="i", yaxs="i")
        xleft <- rep(resPos, nTimepoints) - 0.5
        ybot  <- as.vector(sapply(nTimepoints:1, function(y) rep(y, nAllPos)))
        rect(posRange[1], 1, posRange[2], nTimepoints+1, col="grey", border=NA)
        rect(xleft, ybot, xleft + 1, ybot + 1, col=cols$cols, border=NA)
        tpLab <- sub("T", "", rev(timepoints))
        axis(side=2, at=1:nTimepoints + 0.5, labels=tpLab, las=2, tick=FALSE)        
    }
}
gradientLegend <- function(file=NULL){
    
    # set the output colors
    cols  <- dataToColors()
    cols$d <- cols$d * 100
    nCols <- length(cols$cols)

    # make the legend plot
    xleft <- rep(0, nCols)
    ybot  <- seq(min(cols$d), max(cols$d), length=nCols)
    yinc <- abs(ybot[2] - ybot[1])
    ybot <- ybot - yinc / 2
    ytop <- ybot + yinc
    if(!is.null(file)){
        png(filename = file,
            width = 1.1, height = 1.25,
            pointsize = 8,
            res = 600, units = "in", type = "cairo")            
    }
    par(mar=c(0.2,2.1,0.2,5.1), cex=1, cex.main=1)
    plot(0, 0, typ="n",
         xlim=c(0,1), xlab='', xaxt='n',
         ylim=c(min(ybot),max(ytop)), ylab='', yaxt="n",
         xaxs="i", yaxs="i", las=2)
    rect(xleft, ybot, xleft + 1, ytop, col=cols$cols, border=NA)
    if(!is.null(file)){
        mtext('DSB - Control', # (% of mol.)
              side=2, line=0.5, font=1)
        axis(side=4, las=2)   
        graphics.off()
    } else {
        mtext('DSB - Control (% of mol.)', #
              side=2, line=1, font=2)
        axis(side=4, las=2)        
    }
}

# assemble the data for plotting
getReplicates <- function(){
    reportProgress('getReplicates')
    plotTypes <- c(top='top', bottom='bottom')
    lapply(plotTypes, function(plotType){
        reportProgress(plotType)
        replicateIds <- paste0(plotType, '_', 1:nReplicates)
        replicates   <- sapply(replicateIds, function(replicateId) input[[replicateId]])
        replicateIds <- replicateIds[replicates != '-'] # all non-null samples for a timepoint series     
        list(replicateIds=replicateIds,
             replicates=replicates,
             N=length(replicateIds))
    })  
}
getSeriesData <- function(timepoint, plotType){
    reportProgress('getSeriesData')
    reportProgress(paste(timepoint, plotType))
    nReplicates <- replicates[[plotType]]$N
    stack <- array(NA, dim=c(nStackN, nReplicates, nAllPos), dimnames=list(stackN, NULL, NULL))      
    for(j in 1:nReplicates){
        sampleId <- replicates[[plotType]]$replicateIds[j]
        sample <- input[[sampleId]]
        loadSample(sample, input) # get data from disk on first encounter of a samples
        if(sample %in% badSamples[[timepoint]]){
            stack[stackN$dsb_res,j,] <- NA
        } else {
            sd <- sampleData[[sample]]
            if(timepoint %in% colnames(sd$F_dsb_pos)){
                stack[stackN$dsb_raw,j,] <- sd$F_dsb_pos[allPos,timepoint]                 
                stack[stackN$ctl_bkg,j,] <- predict(sd$fit_ctl[[timepoint]], data.frame(fitPos=allPos))                 
                stack[stackN$dsb_res,j,] <- stack[stackN$dsb_raw,j,] - stack[stackN$ctl_bkg,j,]                    
            } else {
                stack[stackN$dsb_res,j,] <- NA
            }
        }
    }

    # aggregate (i.e. average) all replicates for each timepoint series
    if(nReplicates == 1){
        agg <- stack[stackN$dsb_res,1,]
    } else {
        agg <- apply(stack[stackN$dsb_res,,], 2, mean, na.rm=TRUE)
    }
    agg
}

# function for plotting the individual signal above background
# or difference of that value between two replicate series
makeHeatMap <- function(plotType, tall=FALSE, file=NULL){
    reportProgress('makeHeatMap')
    reportProgress(plotType)
    
    # get the replicate series to aggregate
    replicates <<- getReplicates()
    replicates$difference$N <<- min(replicates$top$N, replicates$bottom$N)
    if(replicates[[plotType]]$N == 0) return(NULL)
    
    # get the list of positions (i.e. PCR product sizes) to retrieve
    allPos  <<- as.integer(input$xmin):as.integer(input$xmax)
    nAllPos <<- length(allPos)    
    allPosK <<- 1:nAllPos
    
    # get the requested timepoints
    timepoints  <<- strsplit(input[['timepoints']], ' +', perl=TRUE)[[1]]
    nTimepoints <<- length(timepoints)

    # get a matrix of aggregated signal above background (i.e. Signal >Control from visualize)
    getSignal <- function(plotType){
        matrix(sapply(timepoints, getSeriesData, plotType),
               ncol=nTimepoints,
               dimnames=list(pos=NULL, timepoint=timepoints))
    }
    main <- ''
    d <- if(plotType == 'difference'){
        if(input[['invertDifference']]) {
            main <- paste(input[['bottom_name']], '-', input[['top_name']])
            getSignal('bottom') - getSignal('top')
        } else {
            main <- paste(input[['top_name']], '-', input[['bottom_name']])
            getSignal('top') - getSignal('bottom')
        }
    } else {
        main <- input[[paste0(plotType, '_name')]]
        getSignal(plotType)
    }
    
    #make plot
    plotHeatMap(d, main, tall, file)
}

