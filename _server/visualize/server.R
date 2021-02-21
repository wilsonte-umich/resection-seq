
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
        lapply(1:nSeries, function(i){
            timepointId <- paste0('timepoint', i)
            updateSelectInput(session, timepointId, selected='T0')
            xBumpId <- paste0('xBump', i)
            updateSelectInput(session, xBumpId, selected=0)
            lapply(1:nReplicates, function(j){
                sampleId <- paste0('sample_', j, i)
                updateSelectInput(session, sampleId, selected='-')
                resetSampleChoices(sampleId, timepointId)
            })
        })        
    }
    observeEvent(input$reset, {
        resetSampleSelectors()
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
    resetSampleChoices <- function(sampleId, timepointId){
        crr <- input[[sampleId]] # keep the timepoint the same when changing samples, if possible
        targetI <- targetIs[input$dsbEnd]
        timepointI <- timepointIs[input[[timepointId]]]
        ss <- samples[[targetI]][[timepointI]] 
        ss <- ss[ss %notin% badSamples[[input[[timepointId]]]]]
        selected <- if(crr %in% ss) crr else '-'
        updateSelectInput(session, sampleId, choices=c('-', ss), selected=selected)        
    }
    lapply(1:nSeries, function(i){
        # respond to timepoint selection in a series row        
        timepointId <- paste0('timepoint', i)
        observeEvent(input[[timepointId]], {
            lapply(1:nReplicates, function(j){
                sampleId <- paste0('sample_', j, i)
                resetSampleChoices(sampleId, timepointId)
            })
        })                  
    })

    # make the output plot
    output$plotOutput <- renderPlot({
        make_plots(input)
    })
    
    # data download
    getExcelName <- function(allele) {
        paste("Resection", allele, Sys.Date(), "xlsx", sep=".")
    }
    getExcelData <- function(allele, file){
        write_xlsx(plotDownload[[allele]], path=file, col_names=TRUE, format_headers=TRUE)  
    }
    output$downloadDsb <- downloadHandler(
        filename = function() getExcelName('dsb'),
        content  = function(file) getExcelData('dsb', file),
        contentType = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )
    output$downloadCtl <- downloadHandler(
        filename = function() getExcelName('ctl'),
        content  = function(file) getExcelData('ctl', file),
        contentType = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )
}

