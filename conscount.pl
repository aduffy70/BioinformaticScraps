#! /usr/bin/perl
# conscount.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Mar 2007

# This script takes a consensus sequence in fasta format with ambiguous characters
# and counts the number of ambiguous characters (i.e., nonconserved characters) in
# each block of 50 nucleotides.


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
print "\nThere are $truelength nucleotides\n";
$count = 0;
$blockcount = 0;
$block = 0;
$nonconserved = 0;
while ($count <= $leng)
{
	while (($block < 50) && ($count <= $leng)) #50 = number of nucleotides in each block
    {
		if ($cons[$count]=~/[^AGCTagct]/)
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











