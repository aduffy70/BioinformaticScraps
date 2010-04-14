#! /usr/bin/perl
# fastanames
# Aaron M Duffy aduffy70{at}gmail.com
# Created March 2010

# This script takes a fasta file and lists the names of the sequences


# Make sure the input genome file exists and can be read
my ($file) = @ARGV;
unless (-e $file && -r $file)
{
    die "fasta file cannot be found or accessed. \n";
}

#Read the file and print the names of the sequences
open FILE, "$file";
while ($line = <FILE>)
{
    if ($line =~/^>(.*)/)
    {
    	print "$1\n";
    }
}
close FILE;
