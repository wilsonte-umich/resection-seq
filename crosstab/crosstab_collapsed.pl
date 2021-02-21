use strict;
use warnings;

# set variables
my ($target, $min, $max) = @ARGV;
my @timepoints = split(" ", $ENV{TIMEPOINTS});
my (%data);

# assemble hash
while(my $line = <STDIN>){
    chomp $line;
    my ($tp, $tgt, $pos, $str) = split("\t", $line, 5);
    $tgt eq $target or next;    
    $str eq '-' or next; # for NovaSeq scripts (was + in MISeq scripts)
    $data{$pos}{$tp}++;
}

# print the table
if ($max == 0) {
    print "\n";
    exit;
}
print join("\t", qw(target position), @timepoints), "\n";
foreach my $pos($min..$max){
    print join("\t", $target, $pos);
    foreach my $tp(@timepoints){
        my $count = $data{$pos} ? ($data{$pos}{$tp} || 0) : 0;
        print "\t$count";
    }
    print "\n";
}   


#ACCT    120_Dex ctl     1033    +       132     TACCTGCGGG
#ACCT    120_Dex ctl     1033    +       125     GGCACTATGT
#ACCT    120_Dex ctl     1033    +       113     TCCCCACTAA,TTCCCACTAA
#ACCT    120_Dex ctl     1033    +       111     TCCCACGCCT,TCCCATGCCT,TTCCACGCCT
#ACCT    120_Dex ctl     1033    +       89      GTGAGAGATA
#ACCT    120_Dex ctl     1033    +       83      GTCTGGATCC
#ACCT    120_Dex ctl     1033    +       77      CAAACAGGTA
#ACCT    120_Dex ctl     1033    +       77      CACATATTAC
#ACCT    120_Dex ctl     1033    +       70      CGATCGTAAA
#ACCT    120_Dex ctl     1033    +       65      ATATCTGCCC
#ACCT    120_Dex ctl     1033    +       57      GTGAATGAAT
#ACCT    120_Dex ctl     1033    +       47      TGGTTGTACC
#ACCT    120_Dex ctl     1033    +       41      GGCCGCCTCC
#ACCT    120_Dex ctl     1033    +       10      ACATTTTGAC

