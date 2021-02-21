
# set target sequence landmarks

targetData[['ILV1-R']] <- list(
    unq = 103, # ~first uniquely mappable base after the sequence common to ctl and dsb
    ctl = list(
        nde=450,  
        mask=c(
            448:451 # RE site plus two proximal and one distal
        )
    ),
    dsb = list(
        dsb=450,
        nde=490, # DUMMY VALUE, is not distal RE site
        mask=c(
            450  # HO dsb itself
        )
    )
)
