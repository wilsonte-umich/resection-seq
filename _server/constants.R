
# data configuration
stackN <- list(
    dsb_raw = 'dsb_raw',
    ctl_raw = 'ctl_raw',
    dsb_bkg = 'dsb_bkg',
    ctl_bkg = 'ctl_bkg',
    dsb_res = 'dsb_res',
    ctl_res = 'ctl_res',
    dsb_rxn  = 'dsb_rxn',
    ctl_rxn  = 'ctl_rxn'
)
nStackN <- length(stackN)

# available plotting colors for data series
# number of sample input rows cascades from here
plotColors <- c(
    blue    =   rgb(0,  0,  225,    maxColorValue=255),
    red     =   rgb(200,0,  0,      maxColorValue=255),
    green   =   rgb(0,  175,0,      maxColorValue=255),
    purple  =   rgb(200,0,  200,    maxColorValue=255)
)
nSeries     <- length(plotColors)
nReplicates <- 4

# function to apply a moving average
ma <- function(x, n){
    filter(x, rep(1 / n, n), sides = 2)
}

# code development utilities for debugging and finding slow steps
verbose <- FALSE
reportProgress <- function(message){
    if(verbose) message(message)
}

# RColorBrewer Palettes
palettes <- list(
    RColorBrewer = list(
        diverging = c(
            "RdBu",            
            "RdYlBu",
            "PuOr",
            "PRGn",
            "PiYG",
            "BrBG"
        )
    )
)

