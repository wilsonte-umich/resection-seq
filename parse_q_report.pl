use strict;
use warnings;

# parse_q_report.pl runs q report on <data-file>.q
# and writes a crosstable of many values and counts by sample to STDOUT
# must be run from the data_scripts directory with ResectionMasterStrainTable.txt

# working variables
my $resectionDir = "/home/wilsonte_lab/clubhouse/data/Wilson/NovaSeq_resection";
my ($command, $sample, %sampleData,
    %sampleValues, %valueNames, %samples);
my %alleles = (ctl=>'Control', dsb=>'DSB');

# collect sample reference and other metadata from ResectionMasterStrainTable.txt
our %refPos;
my $refDir = "$resectionDir/references";
my @metadataOutputCols = qw(
    Experiment
    YWNumber
    DSB_End
    SuperSample
    Condition
    SpikeIn
    Run
); # metadata columns to repeat in our output
# SampleN	Experiment	ExpReplicate	Timecourse	YWNumber	SuperSample	Condition
# Dextrose	DSB_End	SpikeIn	Lab_ID	AGC_ID	Run	Project
# i5Index	i7Sequence	T0	T35	T60	T75	T90
sub setSampleMetadata {
    
    # collect the metadata line for the current working sample
    open my $inH, "<", "ResectionMasterStrainTable.txt" or die "could not open strain table: $!\n";
    my $header = <$inH>;
    chomp $header;
    $header =~ s/\r//g;
    my $i = 0;
    my @header =  split("\t", $header);
    my %header = map { $_ => $i++ } @header;
    while (my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;
        my @f = split("\t", $line);
        if($f[$header{Lab_ID}] and
           $f[$header{Lab_ID}] eq $sample){
            %{$sampleData{$sample}} = map {  $_ => $f[$header{$_}]} @header;
            last;
        } 
    } 
    close $inH;
    
    # collect the position limits for the control and dsb alleles for this sample
    my $ref = $sampleData{$sample}{DSB_End};
    require "$refDir/$ref/$ref.pl";
    
    # remap certain values from ResectionMasterStrainTable to our output table
    # due to name ordering, these appear at the top of the output
    my $j = 100;
    foreach my $valueName(@metadataOutputCols){
        my $vn = "$j\_$valueName";
        $sampleValues{$vn}{$sample} = $sampleData{$sample}{$valueName};
        $valueNames{$vn}++;
        $j++;
    }
}

# INLINE CODE: analyze q report at top level; act on individual jobs with information we need
open my $inH, "-|", "q report -j all NovaSeqResection.q" or die "could not open q report: $!\n";
while(my $line = <$inH>){
    if($line =~ m|^q: target script:.+/(\w+)_(.+).sh|){ # the q script name signal job start and contains the info we need
        ($command, $sample) = ($1, $2);
        $samples{$sample}++;
        setSampleMetadata();
        $command eq 'align' and extractAlign();
        $command eq 'collapse' and extractCollapse();
        $command =~ m/^crosstab/ and extractCrosstab();
    }
}
close $inH;

# INLINE CODE: print the final results, after entire q report has been parsed
my @valueNames = sort {$a cmp $b} keys %valueNames; # sorted by the number prefix for output ordering
my @samples = sort {$a cmp $b} keys %samples; # sample for now are alphabetical, so WT sort together, etc.
print join("\t", '', @samples), "\n";
foreach my $valueName(@valueNames){
    print $valueName =~ m/^\d+_(.+)/ ? $1 : $valueName; # strip the number prefix used for row sorting
    foreach my $sample(@samples){
        my $value = defined $sampleValues{$valueName}{$sample} ? commify($sampleValues{$valueName}{$sample}) : "NA";
        print "\t$value";
    }
    print "\n";
}

# collect read counts and alignments and calculate percentages
sub extractAlign {
    my %values;
    while(my $line = <$inH>){
        finishSampleCommand($line, \%values) and return;
        if($line =~ m/(\d+)\s+uncollapsed.+read sequences/){
            $values{'201_InputReads'} = $1;
        } elsif($line =~ m/reads processed: (\d+)/){            
            $values{'202_UniqueReads'} = $1;
            $values{'203_PercentUnique'} = int(( $values{'202_UniqueReads'} / $values{'201_InputReads'} ) * 1000 + 0.5) / 10 ;
        } elsif($line =~ m/reads with at least one reported alignment:\s+(\d+)/){
            $values{'204_AlignedReads'} = $1;
            $values{'205_PercentAligned'} = int(( $values{'204_AlignedReads'} / $values{'202_UniqueReads'} ) * 1000 + 0.5) / 10 ;
        } elsif($line =~ m/reads that failed to align:\s+(\d+)/){
            $values{'206_UnalignedReads'} = $1;
            $values{'207_PercentUnaligned'} = int(( $values{'206_UnalignedReads'} / $values{'202_UniqueReads'} ) * 1000 + 0.5) / 10 ;
        } elsif($line =~ m/reads with alignments suppressed due to -m:\s+(\d+)/){
            $values{'208_MultimappedReads'} = $1;
            $values{'209_PercentMultimapped'} = int(( $values{'208_MultimappedReads'} / $values{'202_UniqueReads'} ) * 1000 + 0.5) / 10 ;
        } elsif($line =~ m/(\d+)\s+unique molecule IDs/){
            $values{'210_UniqueUMIs'} = $1;
            $values{'211_PercentUniqUMI'} = int(( $values{'210_UniqueUMIs'} / $values{'204_AlignedReads'} ) * 1000 + 0.5) / 10 ;
        } elsif($line =~ m/(\d+)\s+unique best alignments/){
            $values{'212_UniqueMolecules'} = $1;
            $values{'213_PercentUniqMol'} = int(( $values{'212_UniqueMolecules'} / $values{'204_AlignedReads'} ) * 1000 + 0.5) / 10 ;
        }
    }   
}

# sum position counts into aggregate counts by sample and timepoint
sub extractCollapse { # read the q report lines
    my %values;
    while(my $line = <$inH>){
        finishSampleCommand($line, \%values, \&finishCollapse) and return;
        if($line =~ m/^(T\d+)\s+(dsb|ctl)\s+(\d+)\s+-\s+\((\d+) UMIs\)/){
            my ($timepoint, $allele, $pos, $count) = ($1, $2, $3, $4);
            if($pos >= $refPos{$sampleData{$sample}{DSB_End}}{countStart} and 
               $pos <= $refPos{$sampleData{$sample}{DSB_End}}{countEnd}{$allele} ){
                my $alleleOut = $alleles{$allele};
                $values{"$timepoint\_$alleleOut"} += $count;
            } 
        }
    }     
}
sub finishCollapse { # calculate new composite values
    my ($values) = @_;
    foreach my $key(keys %$values){
        if($key =~ m/^(T\d+)_Control/){
            my $dsb = $$values{"$1\_DSB"};
            my $ctl = $$values{"$1\_Control"};
            my $ctlCorrected = int($ctl / ($sampleData{$sample}{SpikeIn} / 100) + 0.5);
            $$values{"$1\_ControlCorrected"} = $ctlCorrected;
            $$values{"$1\_PercentResected"} = int(( $dsb / $ctlCorrected ) * 1000 + 0.5) / 10 ;
        }
    }
}

# calculate the control background shear level
# get from the crosstab product file (not from q report directly)
sub extractCrosstab { # read the q report lines
    my %values;
    while(my $line = <$inH>){
        finishSampleCommand($line, \%values, \&finishCrosstab) and return;
    }     
}
sub finishCrosstab { # read the crosstab file generated by successful jobs
    my ($values) = @_;
    my $crosstabFile = "$resectionDir/output/$sample/data/$sample.crosstab.ctl.txt";
    open my $inH, "<", $crosstabFile or die "could not open $crosstabFile: $1\n";
    my $header = <$inH>;
    chomp $header;
    $header =~ s/\r//g;
    my $i = 0;
    my @header = split("\t", $header);
    my %header = map { $_ => $i++ } @header;
    my %counts;
    my $countStart = $refPos{$sampleData{$sample}{DSB_End}}{countStart};
    my $countEnd = $refPos{$sampleData{$sample}{DSB_End}}{countEnd};
    my $bkgEnd = $refPos{$sampleData{$sample}{DSB_End}}{backgroundEnd};
    while(my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;
        my @f = split("\t", $line);
        $f[$header{position}] >= $countStart or next;
        foreach my $i(2..$#header){
            my $n = $f[$i];
            $f[$header{position}] <= $countEnd and $counts{$header[$i]}{ctlN} += $n;
            $f[$header{position}] <= $bkgEnd   and $counts{$header[$i]}{bkgN} += $n;  
        }
    }
    close $inH;
    my $nBkgPos = $bkgEnd - $countStart + 1;
    foreach my $i(2..$#header){ # background level in control is a prime indicator of library quality
        my $bkg = ($counts{$header[$i]}{bkgN} / $counts{$header[$i]}{ctlN}) / $nBkgPos;
        $$values{"$header[$i]\_BackgroundControl"} = int($bkg * 1e6 + 0.5) / 1e6;
    }
}

# if a job succeeded, record the job values for subsequent printing
sub finishSampleCommand {
    my ($line, $values, $finishSub) = @_;
    if($line =~ m|^q: exit_status: (\d+)|){ # use the exit status line (always last) as our clue that we reached the end of a job's log
        if($1 == 0){ # only commit successful jobs
            $finishSub and &$finishSub($values); # call an action-specific sub that calculates additional composite values, if requested
            foreach my $valueName(keys %$values){
                $sampleValues{$valueName}{$sample} = $$values{$valueName}; # most recent successful job will over-write old values
                $valueNames{$valueName}++;
            } 
        }
        return 1; # send a signal that this job should stop reading new input lines
    }  
}

# utility sub for pretty printing
sub commify {
    my ( $sign, $int, $frac ) = ( $_[0] =~ /^([+-]?)(\d*)(.*)/ );
    my $commified = (
        reverse scalar join ',',
        unpack '(A3)*',
        scalar reverse $int
    );
    return $sign . $commified . $frac;
}

