#! /usr/bin/perl
# consaminocount.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Mar 2007

# This script takes a consensus sequence (with X's for ambiguous characters) in fasta format
# and counts the number of ambiguous characters (i.e., nonconserved characters) in
# each block of 50 aminoacids.


# Make sure the input fasta file exists and can be read
my ($fasta) = @ARGV;
unless (-e $fasta && -r $fasta)
{
    die "fasta file cannot be found or accessed. \n";
}

# Read the consensus sequence into a scalar
open FILE, "$fasta";
$headerline = <FILE>;
$consensus = <FILE>;
close FILE;


# Transfer the sequences & alignment markers to arrays
@cons = split //,$consensus;

$leng = length($consensus) - 2;
$truelength = $leng + 1;
print "\nThere are $truelength amino acids\n";
$count = 0;
$blockcount = 0;
$block = 0;
$nonconserved = 0;
while ($count <= $leng)
{
	while (($block < 17) && ($count <= $leng)) #17 = number of amino acids in each block
    {
		if ($cons[$count]=~/[X]/)
        {
			$nonconserved++;
		}
		$block++;
		$count++;
	}
	$blockcount++;
	print "$blockcount $nonconserved \t$block $count \n";
	$block = 0;
	$nonconserved = 0;
}
