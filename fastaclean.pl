#! /usr/bin/perl
# fastaclean
# Aaron M Duffy aduffy70{at}gmail.com
# Created February 2007

# This script takes messy fasta files (with numbers or whitespace) and
# converts them into clean fasta files


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
	        if ($line =~/(>.*)/)
            {
    		    if ($needslf==1)
                {
    		        print "\n";
    		    }
    		    print "$1\n";
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
