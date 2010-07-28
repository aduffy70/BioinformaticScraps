#! /usr/bin/perl
# gettaxanames.pl 
# Extracts the taxon names from a nexus file.
# Aaron M Duffy aduffy70{at}gmail.com
# September 2008
 
# Make sure the input genome file exists and can be read
my ($nexusfile) = @ARGV;
unless (-e $nexusfile && -r $nexusfile){
    die "nexus file cannot be found or accessed. \n";
}

#Read the nexus file and extract the taxon names
open FILE, "$nexusfile";
while ($line = <FILE>){ 
   	if ($line =~/\[(\d+)\]$/){
 		$line = <FILE>;
 		if ($line =~/(\w+)/){
	 		print "$1\n";
	 	}
    }	
}
print "\n";
close FILE;
