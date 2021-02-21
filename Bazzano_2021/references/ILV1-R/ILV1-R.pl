use strict;
use warnings;

# set target sequence landmarks used by parse_q_report.pl
use vars qw(%refPos);
$refPos{'ILV1-R'} = {
	countStart => 103,
	countEnd => {
    		ctl => 459, # RE plus 9 bases at positions with poor RE digestion 529 + 9
    		dsb => 447  # all resected position, not including HO and 2 more internal bases 530 -3
	},
    backgroundEnd => 447
};

1;

