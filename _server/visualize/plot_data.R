
#----------------------------------------------------------------------
# plot_data.R controls the data plotting actions
#----------------------------------------------------------------------

# plot variables
plotPos  <- list()
plotPosK <- list()
plotData <- list()
plotIs   <- integer()
plotMask <- list(dsb=c(), ctl=c())
plotDownload <- list(ctl=data.frame(), dsb=data.frame())

# make one of the ctl or dsb plots
pch.limit   <- 1
pch.valid   <- 19
cex.main    <- 2
cex.text    <- 1.75
cex.legend  <- 1.85
cex.points  <- 0.75
lwd.rule    <- 2
lwd.rule.light <- 1
make_plot <- function(input, ymin, main, ylab, raw, allele, val){
    
    # create the base plot    
    ymax <- as.numeric(input$ymax)
    plot(0, 0, main=main, cex.main=cex.main,
        xlim=range(c(plotPos$dsb, plotPos$ctl)),
        ylim=c(ymin,ymax),
        xlab='Position (bp)',
        ylab=ylab,
        cex.lab=cex.text, cex.axis=cex.text
    )        
    
    # add mask areas underneath plotted points
    for(mask in plotMask[[allele]]) {
        rect(min(mask) - 0.5, ymin,
             max(mask) + 0.5, ymax, col="darkgrey", border=NA)
    }
    if(length(plotIs) > 0){
        
        # add lines underneath plotted points
        cols <- sapply(plotIs, function(i) input[[paste0('color', i)]])
        abline(h=strsplit(input$hLines, ' +')[[1]], lwd=lwd.rule)
        abline(v=strsplit(input$vLines, ' +')[[1]], lwd=lwd.rule.light)
        abline(h=0, lwd=lwd.rule)
        if(raw) for(i in plotIs) lines(plotPos[[allele]], plotData[[i]]$agg[stackN$ctl_bkg,plotPosK[[allele]]], col=cols[i], lwd=lwd.rule)
        
        # add a legend
        lbls <- sapply(plotIs, function(i) plotData[[i]]$lbl )
        legend('topleft', lbls, col=cols, bty="o", pch=pch.valid, cex=cex.legend, bg="white")
        
        # add the requested data series
        plotDownload[[allele]] <<- data.frame(pos=plotPos[[allele]])
        for(i in plotIs){
            y <- pmin(ymax, plotData[[i]]$agg[val,plotPosK[[allele]]])
            if(length(plotPos[[allele]]) == length(y)){
                plotDownload[[allele]][[lbls[i]]] <<- y
                xBump <- as.integer(input[[paste0('xBump', i)]])
                if(is.null(xBump)) xBump <- 0
                lines (plotPos[[allele]] + xBump, y, col=cols[i])
                points(plotPos[[allele]] + xBump, y, col=cols[i],
                       pch=ifelse(y>=ymax, pch.limit, pch.valid), cex=cex.points)                    
            }
        }        
    }
}

# assemble the data array for a set of data series, with possibility for replicates
getSeriesData <- function(input, nPlotPos, timepoint=NULL){
    plotIs <<- c()
    lapply(1:nSeries, function(i){
        if(is.null(timepoint)) timepoint <- input[[paste0('timepoint', i)]]
        replicateIds <- paste0('sample_', 1:nReplicates, i)
        replicates   <- sapply(replicateIds, function(replicateId) input[[replicateId]])
        replicateIds <- replicateIds[replicates != '-'] # all non-null samples for a timepoint series
        nReplicates  <- length(replicateIds)
        if(nReplicates > 0){
            plotIs <<- c(plotIs, i)
            stack <- array(NA, dim=c(nStackN, nReplicates, nPlotPos), dimnames=list(stackN, NULL, NULL))      
            for(j in 1:nReplicates){
                sample <- input[[replicateIds[j]]]
                loadSample(sample, input) # get data from disk on first encounter of a samples
                sd <- sampleData[[sample]]
                stack[stackN$dsb_raw,j,plotPosK$dsb] <- sd$F_dsb_pos[plotPos$dsb,timepoint]                 
                stack[stackN$ctl_raw,j,plotPosK$ctl] <- sd$F_ctl_pos[plotPos$ctl,timepoint]
                stack[stackN$dsb_bkg,j,plotPosK$dsb] <- predict(sd$fit_ctl[[timepoint]], data.frame(fitPos=plotPos$dsb)) # dsb background determined from control allele lm
                stack[stackN$ctl_bkg,j,plotPosK$ctl] <- predict(sd$fit_ctl[[timepoint]], data.frame(fitPos=plotPos$ctl))                
                stack[stackN$dsb_res,j,] <- stack[stackN$dsb_raw,j,] - stack[stackN$dsb_bkg,j,]
                stack[stackN$ctl_res,j,] <- stack[stackN$ctl_raw,j,] - stack[stackN$ctl_bkg,j,]
            }
            
            # apply moving average to each replicate prior to aggregation
            if(input$movingAverage != "-"){
                maN <- as.integer(input$movingAverage)
                for(j in 1:nReplicates){
                    for(val in stackN) stack[val,j,] <- ma(stack[val,j,], maN)
                }
            }
            
            # aggregate (i.e. average) all replicates for each timepoint series
            agg <- array(NA, dim=c(nStackN, nPlotPos), dimnames=list(stackN, NULL))
            if(nReplicates == 1){
                agg = stack[,1,]
            } else {
                for(val in stackN) agg[val,] <- apply(stack[val,,], 2, mean)
            }
            
            # collect the final required information to plot each timepoint series
            list(
                lbl = paste(timepoint, paste(replicates[replicates != '-'], collapse=",")),
                agg = agg
            )
        } else {
            NULL
        }
    })
}

# the master plotting function
make_plots <- function(input){
    
    # get the list of positions (i.e. PCR product sizes) to retrieve
    # omit the positions not meaningful for plotting
    pp <- as.integer(input$xmin):as.integer(input$xmax)
    plotPos <<- list(
        dsb = pp[!(pp %in% plotMask$dsb)],
        ctl = pp[!(pp %in% plotMask$ctl)]
    )
    plotPosK <<- list(
        dsb = 1:length(plotPos$dsb),
        ctl = 1:length(plotPos$ctl)  
    )
    nPlotPos <- max(length(plotPos$dsb), length(plotPos$ctl))
    
    # collect and aggregate the data to plot
    plotData <<- getSeriesData(input, nPlotPos)
    plotDownload <<- list(ctl=data.frame(), dsb=data.frame())
    
    # make stacked dsb and ctl plots
    par(mfrow=c(3,1))
    if(input$plotType == "raw"){
        make_plot(input, 0, 'DSB Raw',     'Fraction of Alleles', TRUE, 'dsb', stackN$dsb_raw)
        make_plot(input, 0, 'Control Raw', 'Fraction of Alleles', TRUE, 'ctl', stackN$ctl_raw)        
    } else if(input$plotType == "minusCtl"){
        mins <- c() # dynamically adjust minY to show all negative data points (maxY set by user)
        for(i in plotIs) mins <- c(mins, min(plotData[[i]]$agg[c(stackN$dsb_res, stackN$ctl_res),], na.rm=TRUE))
        make_plot(input, min(mins), 'DSB Residual',     'Signal Above Control', FALSE, 'dsb', stackN$dsb_res)
        make_plot(input, min(mins), 'Control Residual', 'Signal Above Control', FALSE, 'ctl', stackN$ctl_res)
    } else {
        # collect T0 data and subtract it from a resection profile
        T0data <- getSeriesData(input, nPlotPos, 'T0')
        mins <- c()
        for(i in plotIs) {
            plotData[[i]]$agg[stackN$dsb_rxn,] <<- plotData[[i]]$agg[stackN$dsb_res,] - T0data[[i]]$agg[stackN$dsb_res,]
            plotData[[i]]$agg[stackN$ctl_rxn,] <<- plotData[[i]]$agg[stackN$ctl_res,] - T0data[[i]]$agg[stackN$ctl_res,]
            mins <- c(mins, min(plotData[[i]]$agg[c(stackN$dsb_rxn, stackN$ctl_rxn),], na.rm=TRUE))
        }
        make_plot(input, min(mins), 'DSB Resection',     'Signal Above T0', FALSE, 'dsb', stackN$dsb_rxn)
        make_plot(input, min(mins), 'Control Resection', 'Signal Above T0', FALSE, 'ctl', stackN$ctl_rxn) 
    }
}


