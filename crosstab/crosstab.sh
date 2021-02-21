#!/bin/bash

#q    require $SAMPLE $DATA_DIR $CTL_LEN $DSB_LEN
#q    require $OUT_PREFIX

#$    -N  crosstab_$SAMPLE
#$    -wd $DATA_SCRIPTS_DIR
#$    -l  vf=4G

# report work to be done
echo "tabulating SAMPLE"

# crosstab by target, collapse UMIs
function crosstab_collapsed {
    TARGET=$1
    MIN=$2
    MAX=$3
    
    # set files
    IN_FILE=$OUT_PREFIX.collapsed.*.txt.gz
    OUT_FILE=$OUT_PREFIX.crosstab.$TARGET.txt

    # crosstab just targets, after UMI collapse
    echo "crosstab $TARGET, collapsed UMIs"
    gunzip -c $IN_FILE |
    perl $ACTION_SCRIPTS_DIR/crosstab/crosstab_collapsed.pl $TARGET $MIN $MAX |
    slurp -o $OUT_FILE
    checkPipe
}

crosstab_collapsed ctl 1 $CTL_LEN
crosstab_collapsed dsb 1 $DSB_LEN

echo
echo "done"

