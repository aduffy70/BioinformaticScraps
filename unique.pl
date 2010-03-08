#! /usr/bin/perl
# unique.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Nov 2007

# This script finds the unique items in 2 files containing a list of genes.  The files must be lists with a single gene name starting each row (the output from my listgenes.pl script for example).  Any other information on the line will be ignored.  This script does not compare the files very efficiently so it may take a while for very long files.  Sorry.  I am a biologist, not a computer scientist...

#Verify that 2 files are on the command line.
if ($#ARGV != 1)
{
    die "Please supply two files containing lists of genes\n";
}

# Read 2 filenames from the command line
my ($file1, $file2) = @ARGV;

# Make sure the files exist and can be read
unless (-e $file1 && -r $file1 && -e $file2 && -r $file2)
{
    die "One or both of your files cannot be found or accessed.\n";
}
my @nomatch1;
my @nomatch2;
#Open the files and read the contents into two hashes
open FILE1, "$file1";
open FILE2, "$file2";
$count=0;
while ($line = <FILE1>)
{
    if ($line =~ /^(\S+)/)
	{
		$array1[$count] = $1;
    	$count++;
    }
}
$count=0;
while ($line = <FILE2>)
{
    if ($line =~ /^(\S+)/)
	{
		$array2[$count] = $1;
    	$count++;
    }
}

# Find and print items in the first file that are not in the second file
print "\n\nGenes in $file1 that are not in $file2\n\n";
foreach $test1 (@array1)
{
	$present=0;
	foreach $test2 (@array2)
	{
		if (($test1=~$test2) && ($test2=~$test1))
		{
			$present=1;
		}
	}
	if ($present<1)
	{
		push @nomatch1, $test1;
	}
}
print "@nomatch1\n\n";

# Find and print items in the second file that are not in the first file
print "Genes in $file2 that are not in $file1\n\n";
foreach $test1 (@array2)
{
	$present=0;
	foreach $test2 (@array1)
	{
		if (($test1=~$test2) && ($test2=~$test1))
		{
			$present=1;
		}
	}
	if ($present<1)
	{
		push @nomatch2, $test1;
	}
}
print "@nomatch2\n\n";

# Close the files
close FILE1;
close FILE2;
