#! /usr/bin/perl
# fastarename
# Aaron M Duffy aduffy70{at}gmail.com
# Created March 2008

# This script takes fasta files with Genbank names (gi|xxxxxxx|) and coverts the names to 4 letters for the genus and 4 letters for the species epithet.  This is a quick and dirty modification of my fastaclean script.

# Make sure the input genome file exists and can be read
my ($dirtyfile) = @ARGV;
unless (-e $dirtyfile && -r $dirtyfile)
{
    die "dirty fasta file cannot be found or accessed. \n";
}

#Read the dirty file and clean it up
open FILE, "$dirtyfile";
my $needslf = 0;
while ($line = <FILE>)
{
    if ($line =~/^$/)
	{
    	#leave out blank lines
    }
	else
	{
	    if ($line =~/^\s+$/)
		{
            #leave out lines with nothing but white space characters
	    }
		else
		{
	    	if ($line =~/(>gi\|.*\|)\s(.{4})\S*\s(.{4})/)
			{
    			if ($needslf==1)
				{
    				print "\n";
    			}
    			print ">$2_$3\n";
    			$needslf=1;
 			}
			else
			{
				for ($line)
				{
					s/^\s+//;  #remove leading whitespace
					s/\s+$//;  #remove trailing whitespace
					s/\s+//g;  #remove internal whitespace
					s/\d+//;   #remove numbers
                    s/\W+//g;   #remove random characters
				}
 				print "$line";
 			}
 		}
    }
}
print "\n";
close FILE;
