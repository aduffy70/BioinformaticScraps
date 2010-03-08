#! /usr/bin/perl
# listgenes2.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created January 2007

# This script filters ncbi genome files and returns just the gene names, locations, and orientations (1-forward, 0-reverse)

# ***WARNING***
# This script will miss any exons beyond the first two.
# Manually check for genes with 2 or more introns!


# Make sure the input genome file exists and can be read
my ($genome) = @ARGV;
unless (-e $genome && -r $genome)
{
    die "genome file cannot be found or accessed. \n";
}

#Read the genome file and print the gene names and locations
open FILE, "$genome";
while ($line = <FILE>)
{
    if ($line =~/^\s+CDS\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ || $line=~/^\s+tRNA\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ || $line=~/^\s+rRNA\D+(\d+)\.\.(\d+)\D(\d+)\.\.(\d+)/ )
	{ #It has an intron
		$start=$1;
		$end=$2;
		$start2=$3;
		$end2=$4;
		if ($line=~/complement/)
		{ #It is on the reverse strand
			$direction=0;
		}
		else
		{ #It is on the forward strand
			$direction=1;
		}
		$line =<FILE>;
		if ($line=~/gene="(\S+)"/)
		{
	    	$genename=$1;
	    	$line = <FILE>;
	    	while ($line!~/product/)
			{
	    		$line = <FILE>;
	    	}
	    	$line=~/product="(.*)"/;
	    	if ($direction==1)
			{
	    		print "$genename\t$start\t$end\t$direction\t$1\n$genename\t$start2\t$end2\t$direction\t$1\n";
	    	}
	    	else
			{
	    		print "$genename\t$end2\t$start2\t$direction\t$1\n$genename\t$end\t$start\t$direction\t$1\n";
	    	}
		}
		elsif ($line=~/product="(.*)"/)
		{
	    	if ($direction ==1)
			{
	    		print "$1\t$start\t$end\t$direction\t$1\n$1\t$start2\t$end2\t$direction\t$1\n";
			}
			else
			{
				print "$1\t$end2\t$start2\t$direction\t$1\n$1\t$end\t$start\t$direction\t$1\n";
			}
		}
    }
    else
	{  #It does not have an intron
	    if ($line=~/^\s+CDS\D+(\d+)\.\.(\d+)/ || $line=~/^\s+tRNA\D+(\d+)\.\.(\d+)/ || $line=~/^\s+rRNA\D+(\d+)\.\.(\d+)/)
		{
			$start=$1;
			$end=$2;
			if ($line=~/complement/)
			{
				$direction=0;
			}
			else
			{
				$direction=1;
			}
			$line =<FILE>;
			if ($line=~/gene="(\S+)"/)
			{
				$genename=$1;
	    		$line = <FILE>;
	    		while ($line!~/product/)
				{
	    			$line = <FILE>;
	    		}
	    		$line=~/product="(.*)"/;
		    	if ($direction == 1)
				{
		    		print "$genename\t$start\t$end\t$direction\t$1\n";
		    	}
		    	else
				{
		    		print "$genename\t$end\t$start\t$direction\t$1\n";
		    	}
			}
			elsif ($line=~/product="(.*)"/)
			{
		    	if ($direction == 1)
				{
		    		print "$1\t$start\t$end\t$direction\t$1\n";
		    	}
		    	else
				{
		    		print "$1\t$end\t$start\t$direction\t$1\n";
		    	}
			}
	    }
	}
}

close FILE;
