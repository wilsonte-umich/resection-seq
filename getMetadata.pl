use strict;
use warnings;

# read sample metadata from ResectionMasterStrainTable.txt

# expects:
#   $ENV{DATA_SCRIPTS_DIR}/ResectionMasterStrainTable.txt
#   $ENV{RUN}
#   $ENV{SAMPLE}

# get parsing instructions
my ($getCol) = @ARGV;

# open the connection
my $METADATA = "$ENV{DATA_SCRIPTS_DIR}/ResectionMasterStrainTable.txt";
open my $inH, "<", $METADATA or die "could not open $METADATA: $!\n";

# parse the header
my $header = <$inH>;
chomp $header;
$header =~ s/\r//g;
my $i = 0;
my @header =  split("\t", $header);
my %header = map { $_ => $i++ } @header;

# return the results
my @results;
if($getCol eq 'TIMEPOINTS' or $getCol eq 'SIDS'){
    getSampleMetadata();
} else {
    getColumnMetadata();
}
print join(" ", @results), "\n";

# clean up
close $inH;

# get a regular headed column, returning all rows matching the current RUN
sub getColumnMetadata {
    while (my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;
        my @f = map { $_ =~ s/^\s+//; $_ =~ s/\s+$//; $_ } split("\t", $line);
        if($f[$header{Run}] and
           $f[$header{Run}] eq $ENV{RUN}){
            push @results, $f[$header{$getCol}];
        } 
    }    
}

# get information on the barcoded timepoint series within a single sample
sub getSampleMetadata {
    my @timepoints = map { $_ =~ m/^T\d+$/ ? $_ : () } @header;
    my @sampleTimepoints;
    my @timepointBarcodes;
    while (my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;
        my @f = map { $_ =~ s/^\s+//; $_ =~ s/\s+$//; $_ } split("\t", $line);
        if($f[$header{Lab_ID}] and
           $f[$header{Lab_ID}] eq $ENV{SAMPLE}){
            foreach my $timepoint(@timepoints){
                if($f[$header{$timepoint}] and
                   $f[$header{$timepoint}] ne 'NA'){
                    push @sampleTimepoints, $timepoint;
                    push @timepointBarcodes, $f[$header{$timepoint}];
                }
            }
        } 
    }
    @results = $getCol eq 'TIMEPOINTS' ? @sampleTimepoints : @timepointBarcodes;
}

