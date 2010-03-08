#! /usr/bin/perl
# pop2struct
# Aaron M Duffy aduffy70{at}gmail.com
# Created September 2007

# This script converts a genepop input file to a structure input file.  It only works on files with a specific number of loci and the same number of digits in each allele record.


# Make sure the input genome file exists and can be read
my ($genepopfile) = @ARGV;
unless (-e $genepopfile && -r $genepopfile)
{
    die "genepop input file cannot be found or accessed. \n";
}

#Read genepop file and convert it to structure input file format
open FILE, "$genepopfile";
while ($line = <FILE>)
{
    if ($line =~/^$/)
    {
    	print"\nblank line\n";
    	#leave out blank lines
    }
    else
    {
    	if ($line =~/(\S+)\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})\s+(\d{3})(\d{3})/)
        {
			print"$1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17\n";
		}
	}
}
close FILE;
