#! /usr/bin/perl
# structure2fasta
# Aaron M Duffy aduffy70{at}gmail.com
# Created March 2007

# This script takes the output of garnier (from the emboss package) and converts it to a
# fasta-like output so secondary structures can be aligned using alignment programs.
# Hey- I'm not claiming this is really useful...


# Make sure the input garnier file exists and can be read
my ($garnier) = @ARGV;
unless (-e $garnier && -r $garnier)
{
    die "garnier structure file cannot be found or accessed. \n";
}

#Read the garnier file and remove everything but the sequence names and structure symbols
open FILE, "$garnier";
my $needlf = 1;
while ($line = <FILE>)
{
	if ($line =~/Sequence:\s+(.+)from:/)
    {
		my $taxon = $1;
		for ($taxon)
        {
			s/\s+$//;
		}
		if ($needlf == 0)
        {
			print "\n";
		}
		print ">$taxon\n";
		$needlf = 0;
	}
	if ($line =~/helix\s(.*)/)
    {
		$helixline = $1;
		$line = <FILE>;
		$line =~/sheet\s(.*)/;
		$sheetline = $1;
		$line = <FILE>;
		$line =~/turns\s(.*)/;
		$turnsline = $1;
		$line = <FILE>;
		$line =~/coil\s(.*)/;
		$coilline = $1;
		@helix = split //,$helixline;
		@sheet = split //,$sheetline;
		@turns = split //,$turnsline;
		@coil = split //,$coilline;
		$count = 0;
		while ($helix[$count])
        {
			if ($helix[$count] =~/H/)
            {
				print "H";
		    }
			if ($sheet[$count] =~/E/)
            {
				print "E";
			}
			if ($coil[$count] =~/C/)
            {
				print "C";
			}
			if ($turns[$count] =~/T/)
            {
				print "T";
			}
			$count++;
		}
	}
}
close FILE;
