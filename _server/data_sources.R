
#----------------------------------------------------------------------
# data_sources.R has path and scripts to load resection data
# it is called on startup, and can be called again by action button(s)
#----------------------------------------------------------------------

#======================================================================
# set NovaSeq_resection file paths
#----------------------------------------------------------------------
strainTableFile <- paste(serverEnv$DATA_SCRIPTS_PATH, "ResectionMasterStrainTable.txt", sep="/")
badSampleFile <- paste(serverEnv$DATA_SCRIPTS_PATH, "BadLibraryTable.txt", sep="/")
outputDir    <- (paste(serverEnv$DATA_PATH, "output",   sep="/"))
targetsDir   <- (paste(serverEnv$DATA_PATH, "references", sep="/")) # reference (for read alignment) and target (for plotting) are synonymous
samplePrefix <- function(sample){
    paste(outputDir, "/", sample, "/data/", sample, sep="")
}
#======================================================================

#======================================================================
# initialize the data sources for the overall app
#----------------------------------------------------------------------
# load avialable references, i.e. target, allele information
targetData <- list()    # this list is filled by sourced reference-definition scripts
targets  <- character() # the known targets, according to strain table
targetIs <- integer()
getTargets <- function(){
    for (target in list.dirs(targetsDir, full.names=FALSE, recursive=FALSE)){
        if(target != 'OLD') source(paste(targetsDir, "/", target, "/", target, ".R", sep=""))
    }
}
#----------------------------------------------------------------------
# load the strain table
strainTable <- data.frame()
timepoints  <- character()
timepointIs <- integer()
allSamples <- character()
samples <- list()
badSamples <- list()
sampleTargets <- character()
loadStrainTable <- function(){
    strainTable <<- read.table(strainTableFile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    badSampleTable <- read.table(badSampleFile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    targets  <<- c(sort(unique(strainTable$DSB_End[strainTable$DSB_End!=''])), 'all')    
    targetIs <<- 1:length(targets)
    names(targetIs) <<- targets    
    allSamples <<- sort(unique(strainTable$Lab_ID))
    cols <- colnames(strainTable)     
    timepoints  <<- cols[grepl('^T\\d+$', cols, perl = TRUE)]
    timepointIs <<- 1:length(timepoints)
    names(timepointIs) <<- timepoints
    samples <<- lapply(targets, function(target){ # organized by dsb target, then timepoint
        lapply(timepoints, function(timepoint){
            sort(unique(strainTable[strainTable$DSB_End == target & !is.na(strainTable[[timepoint]]), 'Lab_ID']))
        })
    })
    samples[[targetIs['all']]] <<- lapply(timepoints, function(timepoint){ # construct list of sample independent of dsb end
        timepointI <- timepointIs[timepoint]
        tmpSamples <- c()
        for(targetI in targetIs){
            tmpSamples <- c(tmpSamples, samples[[targetI]][[timepointI]]   )
        }
        tmpSamples
    })
    sampleTargets <<- sapply(allSamples, function(sample){ # create a lookup from sample to its DSB target
        strainTable[strainTable$Lab_ID == sample, 'DSB_End'][1]
    })
    badSamples <<- lapply(timepoints, function(timepoint){ # organized by timepoint
        badSampleTable[badSampleTable$Timepoint==timepoint,'Lab_ID']
    })
    names(badSamples) <<- timepoints
}
#======================================================================

#======================================================================
# load ctl/dsb positional crosstab data for a specific sample
#----------------------------------------------------------------------
sampleData <- list()
readSampleFile <- function(sample, target){ # read the position x timepoint crosstab file for a sample
    file <- paste(samplePrefix(sample), ".crosstab.", target, ".txt", sep="")
    read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
#======================================================================

#======================================================================
# top-level function for loading and normalizing all data for a given sample, i.e. time course
#----------------------------------------------------------------------
loadSample <- function(sample, input){ # assemble the complete set of information for a project+sample
    if(is.null(sampleData[[sample]]) & sample != "-"){       
        tgt <- targetData[[ sampleTargets[[sample]] ]] # the reference allele for this sample        
        fitPos <- tgt$unq:(tgt$ctl$nde-3) # the control allele positions used to fit the background signal lm
        ctlPos <- tgt$unq:(tgt$ctl$nde+2) # the control allele positions used to count the number of control alelles in sample
        N_ctl_pos <- readSampleFile(sample, 'ctl') # data.frame of UMI counts by position for all timepoints
        N_dsb_pos <- readSampleFile(sample, 'dsb')
        Perc_ctl_ddPCR <- strainTable[strainTable$Lab_ID == sample, 'SpikeIn'] # external value = ctl/dsb allele ratio
        F_ctl_ddPCR <- Perc_ctl_ddPCR / 100           
        F_ctl_pos <- N_ctl_pos # initialize control-normalized data at all timepoints
        F_dsb_pos <- N_dsb_pos
        fit_ctl <- list()
        for(timepoint in colnames(N_dsb_pos)[3:ncol(N_dsb_pos)]){
            A_ctl <- sum(N_ctl_pos[ctlPos, timepoint]) # _measured_ number of control alleles in sampled cells
            A_dsb <- A_ctl / F_ctl_ddPCR               # _inferred_ number of DSB alleles in sampled cells            
            F_dsb_pos[,timepoint] <- N_dsb_pos[,timepoint] / A_dsb # fraction of dsb and control alleles at each position           
            F_ctl_pos[,timepoint] <- N_ctl_pos[,timepoint] / A_ctl
            fit_ctl[[timepoint]]  <- lm(ma(F_ctl_pos[fitPos,timepoint],5)~fitPos) # fit to control to account for random fragmentation AND PCR size bias
        }        
        sampleData[[sample]] <<- list(
            F_ctl_pos = F_ctl_pos,
            F_dsb_pos = F_dsb_pos,
            fit_ctl   = fit_ctl
        )
    }
}
#======================================================================

