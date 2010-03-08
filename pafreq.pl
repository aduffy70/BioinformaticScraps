#! /usr/bin/perl
# pafreq.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created May 2008


# This script summarizes the dominant allele frequencies for each population in an AFLP dataset, where 1 indicates presence of a band and 2 indicates lack of a band.  Sample names must be in the form taxon-population-individual.  Data must be sorted so each population is together.  Data must be whitespace delimited with a header row containing the names of the loci.  Missing data may cause problems.

#Verify that a filename is on the command line.
if ($#ARGV != 0)
{
    die "Please supply a file containing AFLP data\n";
}

# Read a filename from the command line
my ($file1) = @ARGV;

# Make sure the file exists and can be read
unless (-e $file1 && -r $file1)
{
    die "Data file cannot be found or accessed.\n";
}

#Open the file and read the allele names into a scalar and then into an array
open FILE, "$file1";
$headerline = <FILE>;
#print "$headerline\n";
@locusnames = split(" ",$headerline);
$leng = @locusnames;
print "Population @locusnames\n";


#Read and process the rest of the lines
$firsttime = 0;
while ($line = <FILE>)
{
	if ($line =~/^(\S+\-\S+)\-\S+\s+(.+)/)
	{
		@tempphenotypevalues = split(" ",$2);
		if ($1 =~ $pop)
		{
			$count = 0;
			while ($count < $leng)
			{
				if ($tempphenotypevalues[$count] == 1)
				{
					$phenotypevalues[$count]++;
				}
				$count++;
			}
			$samplecount++;
		}
		else
		{
			if ($firsttime == 0)
			{
				$firsttime = 1;
			}
			else
			{
				$count = 0;
				while ($count < $leng)
				{
					$phenotypefreq[$count] = $phenotypevalues[$count]/$samplecount;
					$allelefreq[$count] = 1-((1-$phenotypefreq[$count])**0.5);
					$count++;
				}
				#print "$pop @phenotypefreq\n";
				print "$pop @allelefreq\n";
			}
			$pop = $1;
			$count = 0;
			while ($count < $leng)
			{
				if ($tempphenotypevalues[$count] == 1)
				{
					$phenotypevalues[$count] = 1;
				}
				else
				{
					$phenotypevalues[$count] = 0;
				}
				$count++;
			}
			$samplecount = 1;
		}
	}
}
$count = 0;
while ($count < $leng)
{
	$phenotypefreq[$count] = $phenotypevalues[$count]/$samplecount;
	$allelefreq[$count] = 1-((1-$phenotypefreq[$count])**0.5);
	$count++;
}
#print "$pop @phenotypefreq\n";
print "$pop @allelefreq\n";
close FILE;
