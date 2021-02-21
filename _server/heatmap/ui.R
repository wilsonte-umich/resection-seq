
#----------------------------------------------------------------------
# ui.R defines the html page layout and defines user inputs
#----------------------------------------------------------------------

wideHeight <- '300px'
tallHeight <- '600px'
tallWidth  <- '200px'

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
            textInput(inputId = "timepoints",
                      label = "Timepts (space delim)",
                      value = "T0 T35 T60 T75 T90"),
            hr(),
            textInput(inputId = "xmin",
                      label = "X Min",
                      value = 300),
            textInput(inputId = "xmax",
                      label = "X Max",
                      value = 530),            
            textInput(inputId = "zmax",
                      label = "Z Max",
                      value = 0.001),
            hr(),
            selectInput(inputId = "rcbPalette",
                        label = tags$a(href="https://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/0101-rcolorbrewer-palette-rcolorbrewer-palettes-colorblind-friendly-1.png",
                                       target="RColorBrewer",
                                       "Color Palette"),
                        choices = palettes$RColorBrewer$diverging),
            checkboxInput(inputId = "invertPalette",
                          label = "Invert Palette?",
                          value = TRUE),
            checkboxInput(inputId = "invertDifference",
                          label = "Invert Difference?",
                          value = FALSE),            
            selectInput(inputId = "paletteGradient",
                        label = "Palette Gradient",
                        choices = c('linear','squared','sqrt'),
                        selected = 'linear'),            
            selectInput(inputId = "colorsPerSide",
                        label = "# Colors Per Side",
                        choices = c(5,10,20,50,100),
                        selected = 20),
            hr(),
            textInput(inputId = "saveFontSize",
                      label = "Font (points)",
                      value = 8),
            textInput(inputId = "saveWidth",
                      label = "Width (inches)",
                      value = "1.25"),
            textInput(inputId = "saveHeight",
                      label = "Height (inches)",
                      value = "4"),
            hr(),
            actionButton(inputId = "reload", "Reload Sample Names"),
            actionButton(inputId = "reset", "Reset Form")   ,
            width=2
        ),
        
        # one row of sample selector inputs for each available data series
        mainPanel(
            fluidRow( # input "table" headers
                style="font-weight: bold;",               
                column(2, "Name"),
                column(2, "Top 1"),
                column(2, "Top 2"),
                column(2, "Top 3"),
                column(2, "Top 4")
            ),
            fluidRow(
                column(2, textInput(inputId = "top_name", label = NULL, value = "Top")),
                column(2, selectInput(inputId =  'top_1',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'top_2',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'top_3',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'top_4',
                    label = NULL, choices = c("-")))
            ),
            fluidRow( # input "table" headers
                style="font-weight: bold;",
                column(2, "Name"),
                column(2, "Bottom 1"),
                column(2, "Bottom 2"),
                column(2, "Bottom 3"),
                column(2, "Bottom 4")
            ),
            fluidRow(
                column(2, textInput(inputId = "bottom_name", label = NULL, value = "Bottom")),
                column(2, selectInput(inputId =  'bottom_1',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'bottom_2',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'bottom_3',
                    label = NULL, choices = c("-"))),
                column(2, selectInput(inputId =  'bottom_4',
                    label = NULL, choices = c("-")))
            ),
            
            # the output plot
            fluidRow(
                column(3, plotOutput(outputId="topPlotTall", height=tallHeight, width=tallWidth)),
                column(3, plotOutput(outputId="bottomPlotTall", height=tallHeight, width=tallWidth)),
                column(3, plotOutput(outputId="differencePlotTall", height=tallHeight, width=tallWidth)),
                column(3, plotOutput(outputId="gradientLegend", height='250px',    width='150px'))
            ),
            fluidRow(
                column(3, downloadButton("topPlotTallDownload", "download")),
                column(3, downloadButton("bottomPlotTallDownload", "download")),
                column(3, downloadButton("differencePlotTallDownload", "download")),
                column(3, downloadButton("gradientLegendDownload", "download"))
            ),
            plotOutput(outputId="topPlot",        height=wideHeight),
            plotOutput(outputId="bottomPlot",     height=wideHeight),
            plotOutput(outputId="differencePlot", height=wideHeight)

        )   
    )
)
