use strict;
use warnings;

# this script enforces a proper sequence leader configuration
# and extracts the molecule and timepoint index sequences

# expects:
#   read programatically from ResectionMasterStrainTable.txt
#       $ENV{SIDS}        # sample identification barcodes
#       $ENV{TIMEPOINTS}  # time points
#   read from programatically $REF_PREFIX.leader_sequence.txt
#       $ENV{LEADER_SEQUENCE}
#   # set by user in <data-file.q>
#       $ENV{FIXED_LEN_0}      # the last bases of LEADER_SEQUENCE that are used to locate the read
#       $ENV{UMI_LEN}          # number of bases in the molecule identifier
#       $ENV{FIXED_BASES_1}    # sequence of constant base values after UMI
#       $ENV{BARCODE_LEN}      # number of bases in the sample barcode
#       $ENV{FIXED_BASES_2}    # sequence of constant base values after sample barcode
#       $ENV{ALIGNMENT_LENGTH} # the number of genomic bases past the primer used to located the resection point
# where the following values applied to ILV1-L:
#   $FIXED_LEN_0     6 
#   $UMI_LEN         12 
#   $FIXED_BASES_1   A
#   $BARCODE_LEN     6
#   $FIXED_BASES_2   TA
#   $ALIGNMENT_LENGTH 36 
# and the following values applied to ILV1-R:
#   $FIXED_LEN_0     6 
#   $UMI_LEN         10 
#   $FIXED_BASES_1   A
#   $BARCODE_LEN     6
#   $FIXED_BASES_2   TA
#   $ALIGNMENT_LENGTH 36

# some counters
my ($nLine, $nPassedSeq) = (0) x 10;

# parse the timepoint identifiers
my @sids = split(" ", $ENV{SIDS});
my @timepoints = split(" ", $ENV{TIMEPOINTS});
my %sids = map {
    $sids[$_] =~ tr/ACGT/TGCA/; # reverse complement SIDs
    $sids[$_] = reverse $sids[$_];
    $sids[$_] => $timepoints[$_];
} 0..$#sids;
my ($umi, $timepoint, $alignSeq); # for communicating from seq to qual

# declare the read pattern
my $leaderLen = length($ENV{LEADER_SEQUENCE});
my $fix0Len = $ENV{FIXED_LEN_0};
my $fix0    = substr($ENV{LEADER_SEQUENCE}, -$fix0Len);
my $umiLen  = $ENV{UMI_LEN};
my $fix1    = $ENV{FIXED_BASES_1};
my $tiLen   = $ENV{BARCODE_LEN};
my $fix2    = $ENV{FIXED_BASES_2};
my $fix1Len = length($fix1);
my $fix2Len = length($fix2);
my $maxIndel = 2;
my $minVar = $leaderLen - $fix0Len - $maxIndel; # range of bases that must exist upstream of fix0 search space
my $maxVar = $leaderLen - $fix0Len + $maxIndel;
my $alignLen = $ENV{ALIGNMENT_LENGTH};
my $seqRegex = qr/^(.{$minVar,$maxVar})$fix0(.{$umiLen})$fix1(.{$tiLen})$fix2(.{$alignLen})/;
my $fixLen = $fix0Len + $umiLen + $fix1Len + $tiLen + $fix2Len;
my $lenVar = 0; # for communicating from seq to qual

# run the reads
while(my $line = <STDIN>){
    chomp $line;    
    my $lineN = $nLine % 4;
    if($lineN == 1){ # seq
        $timepoint = "";
        if($line =~ m/$seqRegex/){
            $lenVar = length($1) + $fixLen; # leader length accounting for indels
            $umi = $2;
            $timepoint = $sids{$3};
            $alignSeq = $4;            
        }
    } elsif($lineN == 3){ # qual
        if($timepoint){
            $nPassedSeq++;
            $line =~ m/.{$lenVar}(.{$alignLen})/;
            print join("\n",
                '@'.join("_", $timepoint, $umi, $nPassedSeq),
                $alignSeq,
                '+',
                $1
            ), "\n";
        }
    }
    $nLine++;
}

