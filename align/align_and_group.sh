#!/bin/bash

#q    require $SAMPLE $SIDS $TIMEPOINTS
#q    require $LEADER_SEQUENCE $FIXED_LEN_0 $UMI_LEN $FIXED_BASES_1 $BARCODE_LEN $FIXED_BASES_2 $ALIGNMENT_LENGTH
#q    require $TMP_DIR_LRG $FASTQ_FILES $DATA_DIR $REF_PREFIX $OUT_PREFIX

#$    -N  align_$SAMPLE
#$    -wd $DATA_SCRIPTS_DIR
#$    -l  vf=12G

# report work to be done
echo "aligning $SAMPLE reads"
echo
echo $FASTQ_FILES
which bowtie
which groupBy
echo

# set output file
OUT_FILE=$OUT_PREFIX.aligned.txt.gz

# parse leader portion of sequncing reads
# use timpepoint_UMI as read name(plus temporary counter for uniqueness)
gunzip -c $FASTQ_FILES |
perl $ACTION_SCRIPTS_DIR/report/report_count.pl 'uncollapsed (i.e. input) read sequences' 4 |
perl $ACTION_SCRIPTS_DIR/align/split_reads.pl | 

# map simultaneously to target sequences and phiX
bowtie \
-v 3 \
-k 1 \
-m 1 \
--best \
--un  >(gzip -c | slurp -s 100M -o $OUT_PREFIX.unaligned.fq.gz) \
--max >(gzip -c | slurp -s 100M -o $OUT_PREFIX.multi_map.fq.gz) \
--sam \
--sam-nohead \
$REF_PREFIX - |

# LEGACY: count and save phiX alignments
# NB: no longer expect any phiX on NovaSeq
tee >(          
    grep phiX174 |
    perl $ACTION_SCRIPTS_DIR/report/report_count.pl 'uncollapsed reads mapped to PhiX' | 
    gzip -c |
    slurp -s 100M -o $OUT_PREFIX.phiX.sam.gz
) |
grep -v phiX174 |

# count and discard unmapped reads
# NB: bowtie report multi-mappers with flag 0x4 also
# so this counts all _non-productively_ mapped reads
tee >(          
    awk 'and($2,4)==4' |
    perl $ACTION_SCRIPTS_DIR/report/report_count.pl 'uncollapsed unmapped + multi-mapped reads' |
    gzip -c |
    slurp -s 100M -o $OUT_PREFIX.unmapped.sam.gz
) |
awk 'and($2,4)==0' |

# sort by timepoint+UMI (strip away the temporary uniqueness counter)
perl -ne '$_ =~ s/_\d+\t/\t/; print $_' | 
sort -S 6G -T $TMP_DIR_LRG -k1,1 |

# split UMIs and parse to output table of best alignments
perl $ACTION_SCRIPTS_DIR/align/process_bowtie.pl |
perl $ACTION_SCRIPTS_DIR/report/report_count.pl 'unique best alignments' |
gzip -c |
slurp -s 100M -o $OUT_FILE
checkPipe

# report summary counts to log file
echo
echo -e "timepoint\tcount" # by timepoint
gunzip -c $OUT_FILE |
cut -f 2 |
sort -S 8G -T $TMP_DIR_LRG |
groupBy -g 1 -c 1 -o count
checkPipe

echo
echo -e "allele\tstrand\tcount" # by alignment strand
gunzip -c $OUT_FILE |
cut -f 5,7 |
sort -S 8G -T $TMP_DIR_LRG -k1,1 -k2,2 |
groupBy -g 1,2 -c 2 -o count
checkPipe

echo
echo -e "nMismatches\tcount" # by alignment quality
gunzip -c $OUT_FILE |
cut -f 8 |
sort -S 8G -T $TMP_DIR_LRG -k1,1n |
groupBy -g 1 -c 1 -o count
checkPipe

echo "done"

