
# set target sequence landmarks

targetData[['ILV1-L']] <- list(
    unq = 190, # ~first uniquely mappable base after the sequence common to ctl and dsb
    ctl = list(
        nde=529,  # reflects at 527
        mask=c(
            527:530 # NdeI site plus two proximal and one distal
        )
    ),
    dsb = list(
        dsb=530,   
        nde=614, # reflects at 612
        mask=c(
            530,     # HO dsb itself
            588:603, # unmappable A stretch between HO and Nde
            612:615  # NdeI site plus two proximal and one distal
        )
    )
)
