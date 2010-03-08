#! /usr/bin/perl
# listgenesfasta.pl
# Aaron M Duffy aduffy70{at}gmail.com
# listgenesfasta.pl
# Created February 2009 *This is a modification of my listgenes.pl script


# This script filters an ncbi genome file and a genome fasta file and returns a fasta file of the genes.  For genes with 2 exons and an intron, it returns the sequence of both exons with the intron, not separate entries for each exon.

# Known Issues:
# It behaves strangely for genes like ndhB that have one exon at the start of the genome sequence and the other at the end.  You can fix those manually, or cheat by modifying the genbank file.
# Genes with more than two exons will also cause problems.  They won't appear in the output.
# Do fasta files allow line breaks within the sequence?  If your fasta file has line breaks within the sequence, only the data up to the first line break will be captured.  You can use my fastaclean.pl script to remove extra line breaks from a fasta file.


# Make sure the input genbank file and  fasta genome sequence file exist and can be read
my ($genome, $fasta) = @ARGV;
unless (-e $genome && -r $genome && -e $fasta && -r $fasta)
{
    die "genbank genome or fasta genome file cannot be found or accessed. \n";
}

# Read the genome fasta file into an array with one nucleotide in each element of the array
open FILE, "$genome";
open FILE2, "$fasta";
$line = <FILE2>; #Discard the line with the gene name...
$line = <FILE2>; # ... and save the line with the sequence data.
$sequence=$line;
@sequencechar = split //, $sequence;

# Get the start and end positions of each gene and print the sequences in fasta format
while ($line = <FILE>)
{
    # First see if it is a gene with 2 exons...
    if ($line =~/^\s+CDS\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ || $line=~/^\s+tRNA\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ || $line=~/^\s+rRNA\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ )
	{
		$start=$1 - 1;  # Note- the nucleotide positions start at one, but perl arrays start at zero
		$end=$2 - 1;
		$start2=$3 - 1;
		$end2=$4 - 1;
		$line =<FILE>;
		if ($line=~/gene="(\S+)"/)
		{
	    	print ">Alsophila_$1\n@sequencechar[$start..$end2]\n";
		}
		else
		{
	    	print ">Alsophila_unknown\n@sequencechar[$start..$end2]\n";
		}
    }
	else
	{
	    # ... then see if it has just a single exon
	    if ($line=~/^\s+CDS\D+(\d+)\.\.(\d+)/ || $line=~/^\s+tRNA\D+(\d+)\.\.(\d+)/ || $line=~/^\s+rRNA\D+(\d+)\.\.(\d+)/)
		{
			$start=$1 - 1;
			$end=$2 - 1;
			$line =<FILE>;
			if ($line=~/gene="(\S+)"/)
			{
		    	print ">Alsophila_$1\n@sequencechar[$start..$end]\n";
			}
			else
			{
		    	print ">Alsophila_unknown\n@sequencechar[$start..$end]\n";
			}
	    }
	}
}

# Clean up
close FILE;
close FILE2;
