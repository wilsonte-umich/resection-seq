
#----------------------------------------------------------------------
# ui.R defines the html page layout and defines user inputs
#----------------------------------------------------------------------

# page layout
ui <- fluidPage(
    style="padding: 10px;",
    tags$head(
        tags$style(HTML( # hack fixes default Shiny checkbox group indenting
        ".checkbox-inline, .radio-inline  { 
            margin-left: 4px;
            margin-right: 4px;
        }
        .checkbox-inline+.checkbox-inline, .radio-inline+.radio-inline {
            margin-left: 4px;
            margin-right: 4px;
        }
        hr {
            margin-top: 0;
            margin-bottom: 10px;
            border: 0;
            border-top: 1px solid #666;
        }
        #downloadExcel {
            color: rgb(0,0,200);
        }
        .col-sm-1 {
            padding-left: 3px;
            padding-right: 3px;
        }
        .col-sm-2 {
            padding-left: 3px;
            padding-right: 3px;
        }
        "
    ))),
    
    sidebarLayout(
        
        # top level inputs that control the base plot and all data series
        sidebarPanel(
            selectInput(inputId = 'dsbEnd',
                        label = 'DSB End',
                        choices = targets),
            radioButtons(inputId = "plotType",
                         label = "Plot Type",
                         inline = TRUE,
                         choiceNames = c('Raw Signal', 'Signal >Control', 'Signal >T0'),
                         choiceValues = c('raw', 'minusCtl', 'minusT0')),            
            selectInput(inputId = "movingAverage",
                        label = "Moving Avg.",
                        choices = c("-",3,5,7,9)),
            textInput(inputId = "xmin",
                      label = "X Min",
                      value = 300),
            textInput(inputId = "xmax",
                      label = "X Max",
                      value = 565),            
            textInput(inputId = "ymax",
                      label = "Y Max",
                      value = 0.001),
            textInput(inputId = "hLines",
                      label = "Horiz. Lines (space delim)",
                      value = ""),
            textInput(inputId = "vLines",
                      label = "Vert. Lines (space delim)",
                      value = ""),
            hr(),
            downloadButton('downloadDsb', label="downloadDsb"),
            br(),br(),
            downloadButton('downloadCtl', label="downloadCtl"),
            br(),br(),
            hr(),
            actionButton(inputId = "reload", "Reload Sample Names"),
            actionButton(inputId = "reset", "Reset Form")   ,
            width=2
        ),
        
        # one row of sample selector inputs for each available data series
        mainPanel(
            fluidRow( # input "table" headers
                style="font-weight: bold;",
                column(2, "Color"),
                column(1, "X_Bump"),
                column(1, "Timepoint"),                
                column(2, "Sample_1"),
                column(2, "Sample_2"),
                column(2, "Sample_3"),
                column(2, "Sample_4")
            ), 
            lapply(1:nSeries, function(i) {
                fluidRow(
                    column(2, selectInput(inputId =  paste0('color', i),
                        label = NULL, choices = plotColors, selected=plotColors[i])),
                    column(1, selectInput(inputId =  paste0('xBump', i),
                        label = NULL, choices = -5:5, selected=0)),
                    column(1, selectInput(inputId =  paste0('timepoint', i),
                        label = NULL, choices = timepoints)),                    
                    column(2, selectInput(inputId =  paste0('sample_1', i),
                        label = NULL, choices = c("-"))),
                    column(2, selectInput(inputId =  paste0('sample_2', i),
                        label = NULL, choices = c("-"))),
                    column(2, selectInput(inputId =  paste0('sample_3', i),
                        label = NULL, choices = c("-"))),
                    column(2, selectInput(inputId =  paste0('sample_4', i),
                        label = NULL, choices = c("-")))
                )                
            }),
            
            # the output plot 
            plotOutput(outputId="plotOutput", height='1200px')
            
        )   
    )
)

