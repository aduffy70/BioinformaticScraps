#! /usr/bin/perl
# countblastqueries.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created July 2009

# This script finds the number of different blast queries that had hits in tabular blast output.  The list must be sorted by query.
# Usage- I have a list of 7500 blast hits against my adiantum/angiopteris genes database.  I want to know how many of the queries in the original blast input file had at least one match in the database.

use strict;
use warnings;

#Verify that a file is on the command line.
if ($#ARGV != 0)
{
    die "Please supply a tabular blast output file\n";
}

# Read  filename from the command line
my ($file1) = @ARGV;

# Make sure the file exists and can be read
unless (-e $file1 && -r $file1)
{
    die "File cannot be found or accessed.\n";
}

#Open the file
open FILE1, "$file1";


# Declare some scalars
my $count;
my $countadiantum;
my $countangiopteris;
my $previousgene;
my $thisgene;
my $line;

$count=1;
$line = <FILE1>;
if ($line =~ /^(\S+)\s+\S+/)
{
	$previousgene = $1;
	print "$count\t$previousgene\n";
}
while ($line = <FILE1>)
{
    if ($line =~ /^(\S+)\s+\S+/)
	{
		$thisgene = $1;
		if ($thisgene ne $previousgene)
		{
			$count++;
			print "$count\t$thisgene\n";
			$previousgene = $thisgene;
		}
	}
}

# Close the file
close FILE1;
