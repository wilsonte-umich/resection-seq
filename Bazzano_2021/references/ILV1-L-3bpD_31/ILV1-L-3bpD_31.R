
# set target sequence landmarks

targetData[['ILV1-L-3bpD_31']] <- list(
    unq = 190, # ~first uniquely mappable base after the sequence common to ctl and dsb
    ctl = list(
        nde=529,  # reflects at 527
        mask=c(
            527:530 # NdeI site plus two proximal and one distal
        )
    ),
    dsb = list(
        dsb=527,     
        nde=611, # reflects at 609  
        mask=c(
            527,     # HO dsb itself    
            585:600, # unmappable A stretch between HO and Nde  
            609:612  # NdeI site plus two proximal and one distal
        ),
        del=496 # number of bp to the left of the 3 deleted bases
    )
)
