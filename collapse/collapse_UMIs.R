
# this script has been carefully optimized for speed and memory;
# recursive edit distance comparisons over many thousands
# of UMIs can be very computationally expensive!

# skip base positions with more than this many input UMIs
# with too many UMIs at a position the process can be very slow
# and would sature the UMI pool complexity anyway
maxNUmis <- as.integer(Sys.getenv('MAX_N_UMIS'))

# transform a vector of fixed-length kmers into 2-bit base encoding
# only A=0,C=1,G=3,T=2 base characters may be present!
# handles up to 16-mers in a 32 bit integer (obviously)
encode2Bit <- function(kmers){
    k <- nchar(kmers[1])
    n <- length(kmers)
    m <- matrix(sapply(1:n, function(i){ # this makes use of a convenience of ASCII values for A,C,G,T
        bitwAnd(bitwShiftR(utf8ToInt(toupper(kmers[i])), 1), 3)
    }), nrow=k, ncol=n) # bases in rows, instances in columns
    for(i in 1:k){ # bitshift base rows left as needed, 2 bits at a time
        m[i,] <- bitwShiftL(m[i,], 2*(k-i))
    }    
    as.integer(colSums(m)) # sum the bit shifted values to the final integer encoding of each kmer
}
 
# return all possible values for bitwXor(x, y) as applied to 
# 2-bit encoded k-mers where Hamming distance == 1
# identical bases Xor to 0==A, base mismatches are non-zero
getHamming1Values <- function(k){
    kmers <- character()
    bases <- rep("A", k)
    for(i in k:1){
        for(base in c("C","G","T")){
            kmer <- bases
            kmer[i] = base
            kmers <- c(kmers, paste(kmer, collapse=""))
        }
    }
    encode2Bit(kmers)
}

# load data
d           <- read.table(file("stdin"), header=FALSE, sep="\t", stringsAsFactors=FALSE,
                          colClasses=c("character","character","integer","character","integer","character"))
grpCol      <- c('tp', 'tgt', 'pos', 'str')
datCol      <- c('nDup', 'umi')
colnames(d) <- c(grpCol, datCol)
keys        <- apply(d[,grpCol], 1, paste, collapse="\t") # collapse within a sample+timepoint+target+position
d           <- d[,datCol]

# convert UMIs to 2bit encoding for tree construction
umi2Bit     <- encode2Bit(d[,"umi"]) # vector of 2-bit encoded UMIs for a sample+timepoint
hamming1    <- getHamming1Values(nchar(d[1,"umi"])) # vector of possible Xor results for single-base mismatches

# add descending edge nodes connected to the index node to the UMI network
# uses 'directional' model of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/ (see below)
nis   <- 0
is    <- integer() # indexes of the UMIs for a target+position
ntwkJ <- logical() # which of those indexes that are flagged as part of a UMI network
add_edge_nodes <- function(j){ # js are index of is; j is the index node
    candJs <- (j+1):nis # only check nodes of lower rank (fewer total reads) than index node
    edgJ <- j + which(
        !ntwkJ[candJs] & # don't re-process nodes already known to be in the network
        d[is[j],'nDup'] > 2 * d[is[candJs],'nDup'] - 1 & # directional count criteria for node descent
        bitwXor(umi2Bit[is[j]], umi2Bit[is[candJs]]) %in% hamming1 # next nodes must be exactly 1-base mismatches
    )
    if(length(edgJ) > 0){
        ntwkJ[edgJ] <<- TRUE
        for(j_ in edgJ[edgJ < nis]){
            add_edge_nodes(j_) # recurse through UMI network
        }        
    }
} 
 
# combine similar UMIs into single summed output line
collapse_UMIs <- function(is_){ # i = index of d, always a sequential black since d is pre-sorted by key
    posKey <- keys[is_[1]]
    nis    <<- length(is_)    
    message(paste(posKey, paste('(', nis, ' UMIs)', sep=''), sep="\t"))
    if(nis <= maxNUmis | # only process control allele, or positions with manageable no. of UMIs
       grepl("\tctl\t", posKey, fixed=TRUE)){
        is <<- is_ # WORKING (not initial) is
        while(nis > 0){
            ntwkJ    <<- rep(FALSE, nis)
            ntwkJ[1] <<- TRUE # add the UMI with the greatest count to the current network
            if(nis > 1) add_edge_nodes(1) # recurse to find an entire network related to this UMI
            write(
                c(posKey, sum(d[is[ntwkJ],'nDup']), paste(d[is[ntwkJ],'umi'], collapse=",")),
                '',
                ncolumns=3,
                append=TRUE,
                sep="\t"
            )
            is  <<- is[!ntwkJ] # collect any further networks; speed will progresively improve
            nis <<- length(is)
        }
    } else {
        message(paste("\tskipping, more than", maxNUmis, "UMIs!", sep=" "))   
    }
}

# process one timepoint+pos at a time
# i.e. the group of potentially overcounted UMIs
options(expressions=500000) # need more recursion levels at heavily-hit bins??
discard <- aggregate(1:nrow(d), by=list(keys=keys), collapse_UMIs)

#-----------------------------------------------------------------------
# sample input
#-----------------------------------------------------------------------
#T0      ctl     529     -       34      TTTGGTTTTTTG
#T0      ctl     529     -       34      TTTGTGTTTTGG
#T0      ctl     529     -       34      TTTTTTGTTTTT
#T0      ctl     529     -       33      TGTGGGACTTTT
#T0      ctl     529     -       32      TCTAGTTTATGT
#T0      ctl     529     -       32      TGTTTGTGTGTT
#T0      ctl     529     -       32      TTGTTGTGGTTT
#T0      ctl     529     -       32      TTGTTGTTGATT
#T0      ctl     529     -       31      ATTTTGATTGTT
#T0      ctl     529     -       31      GGTTTTGTGATT
#T0      ctl     529     -       31      GTTGGTTGTGGT
#T0      ctl     529     -       31      TTAATTTGTTTT
#T0      ctl     529     -       31      TTGTTAGCGGTT
#T0      ctl     529     -       31      TTGTTATTGTTG
#T0      ctl     529     -       31      TTTGATGTTTTT
#T0      ctl     529     -       31      TTTTGTTTTTTG
#T0      ctl     529     -       31      TTTTTTATTGTT

#-----------------------------------------------------------------------
# sample output
#-----------------------------------------------------------------------
#T0      ctl     529     -       289     TTTGGTTTTTTG,TTTGGTTTTTTA,GTTGGTTTTTTG,TTAGGTTTTTTG,TTTGCTTTTTTG,TTTGGTTTGTTG,TTTGGATTTTTG,TTTGTTTT
#TTTG,TTTGGTTTTTTT,TTTCGTTTTTTG,GTAGGTTTTTTG,TCTGGTTTTTTA,TTAGGTTATTTG,TTTGGTTGTTTG,TTTACTTTTTTG,TTTGGTATGTTG,TTTGGTTTTCTA,TCTGGTTTTTTG,TTTC
#CTTTTTTG,TTTGGATTCTTG,TAAGGTTTTTTG,TTTGGTTATTTG,CTTGGTTGTTTG,GTTGGTTTCTTG,TTGGGTTGTTTG,TTTCGTGTTTTG,TTTGGATTTTTC,TCTGGTTATTTA,TTAGGTATGTTG,
#TTAGGTTCTTTG,TTTGCTTTTTTA,CTTGCTTGTTTG,TAAGCTTTTTTG,TAAGGTTTTATG,TGTGCTTTTTTG,CCTGGTTTTTTA,TCTCCTTTTTTG
#T0      ctl     529     -       256     TTTGTGTTTTGG,TTTGTGTTTAGG,TTTGTGTTTTGA,TTGGTGTTTTGG,TGTGTGTTTTGG,TTTGTGATTTGG,TTCGTGTTTTGG,TTTGTCTT
#TTGG,TTTGAGTTTAGG,TTTGTGTGTTGG,TTTGAGTTTTGG,TTTGTATTTAGG,TTTGTGTCTTGG,TTTGTGTTGAGG,TTTGTGTTTTAA,TTTGTGTTTTCG,TCGGTGTTTTGG,GTCGTGTTTTGG,TGTG
#TGTTATGG,TGTGTGTTTGGG,TTGGCGTTTTGG,TTGGTGCTTTGG,TTGGTGTTTAGG,TTAGAGTTTAGG,TTTGTGGTTAGG,TCGGTGGTTTGG,TTTGAGTGTTGG,TTTGTATCTAGG,ATTGTATTTAGG,
#CCGGTGTTTTGG,GTCGTGTTTTAG,TCTGTGTTTTGA,TGCGTGTTATGG,TTAGAGTTTAGA,TTAGATTTTAGG,TTGGTGTTTAGC,TTTGTATTTAAG,CCGGTGTTTTTG,TTGGTGTTTAAC
#T0      ctl     529     -       269     TTTTTTGTTTTT,TTTTTCGTTTTT,ATTTTTGTTTTT,CTTTTTGTTTTT,TTTTTTCTTTTT,TTTTTTGCTTTT,TTTTTTGTCTTT,TATTTCGT
#TTTT,TTTTGCGTTTTT,CTTTTTGTTTAT,TTTTCCGTTTTT,TTTTTTCTCTTT,TTTTTTGCCTTT,ATTTTTGTTTTC,CATTTTGTTTTT,CCTTTTGTTTTT,CTTTTTGTTTCT,TTATTTCTTTTT,TTCT
#TCGTTTTT,TTTTGTGCTTTT,TTTTTTCTTTGT,TTTTTTGCTTAT,TTTTTTATCTTT,ATTTTTGATTTT,TGTTTCGTTTTT,TTTTTCGTTCTT,CATATTGTTTTT,CCTTTGGTTTTT,CCTTTTGGTTTT,
#CTTTCTGTTTCT,CTTTTAGTTTTT,CTTTTTATTTCT,TGTTTTGCTTAT,TTTTCTGCCTTT,CATTTTGTTTAT,CCTTTTGTTTCT,TTCTTTCTCTTT,TTTTTCATCTTT,ATCTTTCTCTTT
#T0      ctl     529     -       144     TGTGGGACTTTT,TGTGGGATTTTT,TGTGGGGCTTTT,TGTTGGACTTTT,TGGGGGACTTTT,TGTGGGTCTTTT,TCTGGGATTTTT,TAGGGGAC
#TTTT,TGGCGGACTTTT,TGTCGGACTTTT,TGTTGAACTTTT,TGTTGGCCTTTT,TGTGGGGCTATT,TCTGCGATTTTT,TGTTGGACTTGT,TTTGGGACTTTT,TACGGGACTTTT,TGGCGGACTTGT,TGGG
#GGACTTTC,AGTTGGACTTGT,TGGGGGACCTTT
#T0      ctl     529     -       72      TCTAGTTTATGT,TATAGTTTATGT,TCTAGTTTATGC,GCTAGTTTATGT,TCTAGTTTGTGT,TCTAGTTTTTGT,TCTATTTTATGT,ACTAGTTT
#ATGT,TCTAGGTTATGC,TCTAGTTTCTGT,ACTTGTTTATGT,GATAGTTTATGT

#-----------------------------------------------------------------------
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/
#-----------------------------------------------------------------------
#we developed a final method, directional. We generated networks from the UMIs at a single locus, in which directional edges connect nodes
#   a single edit distance apart
#   when na > 2nb - 1,
#where na and nb are the counts of node a and node b.
#
#The entire directional network is then considered to have originated from the node with the highest counts. The ratio between the final counts for the true UMI and the erroneous UMI generated from a PCR error is dependent upon which PCR cycle the error occurs and the relative amplification biases for the two UMIs, but should rarely be less than twofold. The ?1 component was included to account for strings of UMIs with low counts, each separated by a single edit distance for which the 2n threshold alone is too conservative.
#
#This method allows UMIs separated by edit distances greater than one to be merged so long as the intermediate UMI is also observed, and with each sequential base change from the most abundant UMI, the count decreases.
#
#For this method, the number of directional networks formed is equivalent to the estimated number of unique molecules.

