#!/usr/bin/perl

use strict;
use warnings;

#Will need to install HTML::TableExtract Module
use HTML::TableExtract;


#open config file
open my $fhConfig2, '>>', 'config.prm' or die $!;

my $html_file = <*.htm>;

my $te = HTML::TableExtract->new( headers =>['Window Start','Window End'] );
$te->parse_file($html_file);

my @arrayLCR;

#Store parsed htm data
foreach my $ts ($te->tables) {
  foreach my $row ($ts->rows) {
     push @arrayLCR, "@$row\n";
	}
}

#Generate unique array
my @filteredLCR = uniq(@arrayLCR);


#Print number of unique LCRs to config
my $numLCR = scalar @filteredLCR;

print $fhConfig2 "\nNumber_of_LCRs: $numLCR\n";

my $itter = 1;
#Loop through unique array and push to config
foreach my $arrayLCRcycle (@filteredLCR) {
print $fhConfig2 "LCR_range$itter: $arrayLCRcycle";	
$itter++;
}

close $fhConfig2; 


#Unique hash subroutine
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

