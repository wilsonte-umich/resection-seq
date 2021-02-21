#!/bin/bash

#q    require $SAMPLE $TIMEPOINTS $MAX_N_UMIS
#q    require $TMP_DIR_LRG $DATA_DIR $OUT_PREFIX

#$    -N  collapse_$SAMPLE
#$    -wd $DATA_SCRIPTS_DIR
#$    -l  vf=6G
#$    -t 1-$N_SIDS

# report work to be done
getTaskObject TIMEPOINT $TIMEPOINTS
echo "processing $SAMPLE $TIMEPOINT"
echo

# set files
IN_FILE=$OUT_PREFIX.aligned.txt.gz
OUT_FILE=$OUT_PREFIX.collapsed.$TIMEPOINT.txt.gz

# collapse near-neighbor UMIs to avoid overcounting
gunzip -c $IN_FILE |
awk 'BEGIN{OFS="\t"}$2=="'$TIMEPOINT'"{print $2, $5,$6,$7, $4, $3}' |
sort -S 4G -T $TMP_DIR_LRG -k1,1 -k2,2 -k3,3n -k4,4 -k5,5nr |
Rscript $ACTION_SCRIPTS_DIR/collapse/collapse_UMIs.R |
gzip -c |
slurp -s 100M -o $OUT_FILE
checkPipe

echo "done"

#AAAAAATATCAAAGAAAAAGAGTCATCTCAAACATA    T0      AAAAAGTTGCTA    1       dsb     614     -       0       36
#GAATTCATCGAGCGATATTCTATCCTGAAATACATA    T0      AAAAGATTATTA    1       ctl     529     -       0       36

#grpCol      <- c('tp', 'tgt', 'pos', 'str')
#colnames(d) <- c(grpCol, 'nDup', 'umi')

