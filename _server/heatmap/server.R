
#----------------------------------------------------------------------
# server.R returns a function that defines reactive server actions
#----------------------------------------------------------------------

# handle user actions and update output
`%notin%` <- Negate(`%in%`)
server <- function(input, output, session){
    sessionEnv <- environment()
    verbose <<- FALSE

    # declare session-specific data resources; only visible within server block
    # see: https://shiny.rstudio.com/articles/scoping.html
    for(script in c(
        "plot_data.R"
    )) source(paste(serverEnv$ACTIONS_PATH, script, sep="/"), local=TRUE)
    
    # load all references and samples again when requested (e.g. when new data is available)
    observeEvent(input$reload, { initalizeApp() })
    
    # clear the current plot selections
    resetSampleSelectors <- function(){
        reportProgress('resetSampleSelectors')
        lapply(1:nReplicates, function(j){
            lapply(c('top','bottom'), function(type){
                sampleId <- paste0(type, '_', j)
                updateSelectInput(session, sampleId, selected='-')
                resetSampleChoices(sampleId, crr='-')                
            })
        })       
    }
    observeEvent(input$reset, {
        reportProgress('input$reset')
        resetSampleSelectors()        
        updateTextInput(session, 'top_name', value='Top')
        updateTextInput(session, 'bottom_name', value='Bottom')
        updateTextInput(session, 'timepoints', value='T0 T35 T60 T75 T90')
        updateTextInput(session, 'xmin', value=300)
        updateTextInput(session, 'xmax', value=532)
        updateTextInput(session, 'zmax', value=0.001)
        updateSelectInput(session, 'rcbPalette', selected='RdBu')
        updateCheckboxInput(session, 'invertPalette', value=TRUE)
        updateCheckboxInput(session, 'invertDifference', value=FALSE)
        updateSelectInput(session, 'paletteGradient', selected='linear')
        updateSelectInput(session, 'colorsPerSide', selected=20)  
    })
    
    # initialize the non-plotted position masks
    observeEvent(input$dsbEnd, {
        plotMask <<- list(
            dsb = targetData[[input$dsbEnd]][['dsb']]$mask,
            ctl = targetData[[input$dsbEnd]][['ctl']]$mask
        )
        resetSampleSelectors()  
    })
    
    # handle cascading input values
    resetSampleChoices <- function(sampleId, crr=NULL){
        reportProgress('resetSampleChoices')
        if(is.null(crr)) crr <- input[[sampleId]] # keep the timepoint the same when changing samples, if possible
        targetI <- targetIs[input$dsbEnd]
        ss1 <- samples[[targetI]][[1]]
        ss2 <- samples[[targetI]][[2]]
        ss <- sort(unique(ss1, ss2))
        selected <- if(crr %in% ss) crr else '-'
        updateSelectInput(session, sampleId, choices=c('-', ss), selected=selected)        
    }
    
    # make the output plot
    output$gradientLegend <- renderPlot({
        gradientLegend()   
    })
    output$topPlot <- renderPlot({
        makeHeatMap('top')
    })
    output$bottomPlot <- renderPlot({
        makeHeatMap('bottom')
    })
    output$differencePlot <- renderPlot({
        makeHeatMap('difference')
    })
    output$topPlotTall <- renderPlot({
        makeHeatMap('top', TRUE)
    })
    output$bottomPlotTall <- renderPlot({
        makeHeatMap('bottom', TRUE)
    })
    output$differencePlotTall <- renderPlot({
        makeHeatMap('difference', TRUE)
    })
    # image download
    output$topPlotTallDownload <- downloadHandler(
        filename = function() {
            paste(input$top_name, "tall.png", sep = ".")
        },
        content = function(file) {
            makeHeatMap('top', tall=TRUE, file=file)
        },
        contentType = "image/png" 
    )
    output$bottomPlotTallDownload <- downloadHandler(
        filename = function() {
            paste(input$bottom_name, "tall.png", sep = ".")
        },
        content = function(file) {
            makeHeatMap('bottom', tall=TRUE, file=file)
        },
        contentType = "image/png" 
    )
    output$differencePlotTallDownload <- downloadHandler(
        filename = function() {
            paste(input$top_name, input$bottom_name, "tall.png", sep = ".")
        },
        content = function(file) {
            makeHeatMap('difference', tall=TRUE, file=file)
        },
        contentType = "image/png" 
    )
    output$gradientLegendDownload <- downloadHandler(
        filename = function() {
            paste("legend", "png", sep = ".")
        },
        content = function(file) {
            gradientLegend(file=file)   
        },
        contentType = "image/png" 
    )
}

