use strict;
use warnings;

# group timepoint + UMI + mapped position to declare unique input molecules

# declare variables
my ($prev, @alns) = ("");
my ($nIds, $nDiscarded) = (0) x 100;

# collect all alignments for each timepoint+UMI
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line, 12);
    if($prev and $prev ne $f[0]){ # input is pre-sorted by timepoint+UMI = moleculeID for sample
        printIt();      
        @alns = ();
    }
    $prev = $f[0]; 
    push @alns, \@f;
}                                       
printIt();

# print output file
# ensure that timepoint+position+UMI reports only one best alignment
# and that the alignment is singular
sub printIt {
    @alns or return;
    $nIds++;
    
    # collect all mapped positions claimed by this timepoint+UMI
    my (%nm);
    foreach my $f(@alns){
        my $pos = $$f[3] + length($$f[9]) - 1;        
        my $str = $$f[1] & 16 ? "-" : "+";
        my $aln = "$$f[2]\t$pos\t$str";        
        $$f[11] =~ m/NM:i:(\S+)/ and my $nm = $1 || 0; 
        $nm{$aln}{$nm} or $nm{$aln}{$nm} = $f; # keeps the first encountered equally good alignment
    }     
  
    # finish processing acceptable molecules, i.e. duplicate groups
    # here, let UMI = timepoint+UMI+mapped position
    foreach my $aln(keys %nm){
        my ($nm) = sort {$a <=> $b} keys %{$nm{$aln}}; # keep best alignment as read with the fewest sequencing errors
        my $f = $nm{$aln}{$nm};
        $$f[11] =~ m/MD:Z:(\S+)/ and my $md = $1 || 'NA';
        print join("\t", 
            $$f[9],    
            split("_", $prev), # timepoint, umi
            scalar(@alns),
            $aln, # target, position, strand
            $nm,
            $md
        ), "\n";         
    }
}

print STDERR "$nIds\tunique molecule IDs (timepoint+UMI)\n";
print STDERR "$nDiscarded\tmolecule IDs discarded with >1 alignment\n"; # legacy

#T0_AAAAAAGTAGCT 16      dsb
#575     255     40M     *       0       0
#AGAAAAAAAATATCAAAGAAAAAGAGTCATCTCAAACATA
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF
#XA:i:1  MD:Z:1A38       NM:i:1  XM:i:2

