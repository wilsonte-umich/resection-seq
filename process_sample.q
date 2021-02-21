
# pretty parsing
echo "================================================="

# set needed file paths and extracted information
invoke $ACTION_SCRIPTS_DIR/file_info/file_paths.q

# record the global reference sequence
# should be obsolete given metadata table
$DISCARD  run echo $REFERENCE > $REF_ID_FILE
$DISCARD  run echo $PERCENT_CONTROL > $PERCENT_CTL_FILE

# recover the leader sequence specific to this sample DSB allele
$LEADER_SEQUENCE    run cat $REF_PREFIX.leader_sequence.txt

# align the reads and collapse to UMIs
qsub $ACTION_SCRIPTS_DIR/align/align_and_group.sh
qsub $ACTION_SCRIPTS_DIR/collapse/collapse_UMIs.sh

# analyze the unique reads
qsub $ACTION_SCRIPTS_DIR/crosstab/crosstab.sh

