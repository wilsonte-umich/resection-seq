
# set the project paths
$SAMPLE_DIR     $OUTPUT_DIR/$SAMPLE$SPECIAL_EXTENSION
$DATA_DIR       $SAMPLE_DIR/data
$PLOT_DIR       $SAMPLE_DIR/plots
invoke file/create.q $DIR $SAMPLE_DIR $DATA_DIR $PLOT_DIR

# check the source FASTQ file
# ******** adjust $FASTQ_FILES as needed based on input file naming schema ********
$FASTQ_FILES    $RUNS_DIR/$RUN/$PROJECT/Sample_$SAMPLE_ID/$SAMPLE_ID\_*_R1_001.fastq.gz
#$FASTQ_FILES    $RUNS_DIR/$RUN/$PROJECT/$RUN/$PROJECT/Sample_$SAMPLE_ID/$SAMPLE_ID\_*_R1_001.fastq.gz
invoke file/require.q $FILE $FASTQ_FILES

# collect this sample's timepoint information
$TIMEPOINTS run perl $ACTION_SCRIPTS_DIR/getMetadata.pl TIMEPOINTS
$SIDS       run perl $ACTION_SCRIPTS_DIR/getMetadata.pl SIDS
$N_TIMEPOINTS   run echo $TIMEPOINTS | wc -w
$N_SIDS         run echo $SIDS | wc -w

# set the output file prefixes
$OUT_PREFIX     $DATA_DIR/$SAMPLE
$PLOT_PREFIX    $PLOT_DIR/$SAMPLE
$REF_ID_FILE    $OUT_PREFIX.reference.txt # should be obsolete given metadata table
$PERCENT_CTL_FILE $OUT_PREFIX.percent_control.txt

# set the reference sequence information
$REF_PREFIX     $REFS_DIR/$REFERENCE$SPECIAL_EXTENSION/$REFERENCE
$REF_FILE       $REF_PREFIX.fa
$REF_INFO       $REF_PREFIX.R
invoke file/require.q $FILE $REF_FILE $REF_INFO
$CTL_LEN        run cat $REF_FILE | grep -A 1 '>ctl' | tail -n1 | wc -c
$DSB_LEN        run cat $REF_FILE | grep -A 1 '>dsb' | tail -n1 | wc -c

preserve all

