#! /usr/bin/perl
# adjacentpairs2.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Nov 2007

# This script determines which pairs of adjacent genes in the first file are also present in the second file.  The order of the genes in the pairs matters- they must be transcribed in the same order.  Both files must be formatted with two genenames in each row separated by white space.

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

#Open the files and read the contents into four arrays
open FILE1, "$file1";
open FILE2, "$file2";
$count=0;
while ($line = <FILE1>)
{
    if ($line =~ /^(\S+)\s+(\S+)/)
    {
	    $pair1a[$count] = $1;
	    $pair1b[$count] = $2;
        $count++;
    }
}
$count=0;
while ($line = <FILE2>)
{
    if ($line =~ /^(\S+)\s+(\S+)/)
    {
	    $pair2a[$count] = $1;
	    $pair2b[$count] = $2;
        $count++;
    }
}

# Check whether each pair in the first file is in the second file
print "1=adjacent pair is present\n0=adjacent pair is not present\n";
$count=0;
foreach $onea (@pair1a)
{
	$oneb=$pair1b[$count];
	$present=0;
	$count2=0;
	foreach $twoa (@pair2a)
    {
		$twob=$pair2b[$count2];
		if (($onea=~$twoa) && ($oneb=~$twob))
        {
			$present=1;
		}
		$count2++;
	}
	print "$onea $oneb    $present\n";
	$count++;
}

# Close the files
close FILE1;
close FILE2;
