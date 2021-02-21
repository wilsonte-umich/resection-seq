use strict;
use warnings;

# provide log file information on read counts

my ($msg, $nPerGrp) = @ARGV;
$nPerGrp or $nPerGrp = 1;

my $n = 0;

while(<STDIN>){
    $n++;
    print $_;
}

print STDERR ($n/$nPerGrp)."\t$msg\n";
