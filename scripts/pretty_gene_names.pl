#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;

my $filename = $ARGV[0];

my $fh = new IO::File;
if($fh->open("<$filename")){
    my $header = <$fh>;
    print $header;
    while(my $line = <$fh>){
        my @line = split /,/, $line;
	my $genename = $line[8];
        $genename =~ s/\s*(Gene)?\s*\[.*$//;
	$genename =~ s/\"//g;
	my $line = join ',', @line;
        print $line;
    }
    
    $fh->close;
}
