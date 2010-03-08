#! /usr/bin/perl
# filterfasta.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Aug 2009

# This script takes a file of fasta sequences and a file listing the sequences wanted and pulls just those sequences from the fasta file.

#Verify that 2 files are on the command line.
if ($#ARGV != 1)
{
    die "Please supply a fasta file and a file listing the sequences wanted\n";
}

# Read 2 filenames from the command line
my ($file1, $file2) = @ARGV;

# Make sure the files exist and can be read
unless (-e $file1 && -r $file1 && -e $file2 && -r $file2)
{
    die "One or both of your files cannot be found or accessed.\n";
}

#Open the files and read the list file into an array
open FASTAFILE, "$file1";
open LISTFILE, "$file2";
$count = 0;
while ($line = <LISTFILE>)
{
	if ($line =~ /(.+)/)
	{
		$sequenceswanted[$count] = $1;
		$count++;
	}
}
close LISTFILE;


#Read through the FASTAFILE and extract any sequences that are in the wanted list
while ($line = <FASTAFILE>)
{
	if ($line =~ /^>(.+)/)
	{
		foreach $wanted (@sequenceswanted)
		{
			if ($1 eq $wanted)
			{
				print ">$1\n";
				$line = <FASTAFILE>;
				print "$line";
			}
		}
	}
}
close FASTAFILE;
