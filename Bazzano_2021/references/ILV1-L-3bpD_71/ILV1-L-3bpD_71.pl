use strict;
use warnings;

# set target sequence landmarks used by parse_q_report.pl
use vars qw(%refPos);
$refPos{'ILV1-L-3bpD_71'} = {
	countStart => 190,
	countEnd => {
    		ctl => 538, # NdeI plus 9 bases at positions with poor NdeI digestion
    		dsb => 524  # all resected position, not including HO and 2 more internal bases
	},
    backgroundEnd => 523
};

1;

