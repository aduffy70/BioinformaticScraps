#! /usr/bin/perl
# findprimer.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created Oct-Dec 2006

# This script:
#  1) Performs a multiple nucleotide alignment using clustalw
#  2) Finds all primers with a length range with min conserved 3' and max# of nonconserved nucleotides
#  3) Filters out primers that don't fall within specified Tm and GC content ranges
#  4) Searches against the Lab Database for equivalent primers
#  5) Outputs to standard output, a fasta format file and a amplifx primer list format file
# Input file requirements:
#  1) fasta file with at least two sequences
#  2) text file called LabPrimerDatabase.txt containing a primer database in the format:
#        sequence  reference   any other additional columns
# This script generates the following files (and clobbers them if they already exist):
#   tempalign.aln  *
#   tempprimer.txt  *
#   temprev.txt  *
#   tempgc.txt  *          *Removed after use
#   goodprimers.txt
#   amplifxfile.txt
# The following programs must be in your PATH:
#    clustalw

#Change these values to match your requirements
my $min3prime = 13;      # minimum conserved length at the 3' end of the primer
my $minprimer = 19;     # minimum primer length (actual value is 1+ this number)
my $maxprimer = 26;     # maximum primer length
my $maxmismatch = 3;    # maximum number of mismatched nucleotides (actual value is this number -1)
my $lowTm = 40;         # lowest acceptable Tm
my $highTm = 75;        # highest acceptable Tm
my $lowgc = 0.30;       # lowest acceptable GC content
my $highgc = 0.70;      # highest acceptable GC content

# Make sure the input fasta file exists and can be read
my ($fasta) = @ARGV;
unless (-e $fasta && -r $fasta)
{
    die "fasta file cannot be found or accessed. \n";
}

# Run clustalw alignment
my $clustalwmessage = qx/clustalw $fasta -OUTFILE=tempalign.aln/;

# Verify that clustalw completed successfully and made a valid alignment file
my $alignmentfile = "tempalign.aln";
unless (-e $alignmentfile && -r $alignmentfile)
{
    die "alignment file cannot be accessed\n";
}

# Read one of the sequences & the conserved nucleotide markers (*) into scalars
open FILE, "$alignmentfile";
$headerline = <FILE>;
$headerline = <FILE>;
$headerline = <FILE>;
LINE:while ($line = <FILE>)
{
    if ($line =~ /^$/)
	{
		next LINE;
    }
	else
	{
	    $line =~ /^(\S+\s+)(\S+)/;
	    $sequence .= $2;
	    $taxonname =$1;
	    $namelength = length($taxonname);
	    while ($line =~ /^\S+\s+\S+/)
		{
			$line = <FILE>;
	    }
	    $line =~ /^\s{$namelength}(.+)/;
	    $alignment .= $1;
	    next LINE;
	}
}
close FILE;

# Make a reverse complement copy of the sequence and a reverse copy of the alignment
$rcsequence = reverse $sequence;
$rcsequence =~ tr/ACGTacgt/TGCAtgca/;
$rcalignment = reverse $alignment;

# Transfer the sequences & alignment markers to arrays
@sequ = split //,$sequence;
@align = split //,$alignment;
$alignlength = @sequ;
@rcsequ = split //,$rcsequence;
@rcalign = split //,$rcalignment;

# A better way to find conserved regions
$primerlength = 0;
$bp = $alignlength;
$threeprime = $alignlength;
$threeprime--;
$primercount = 0;
while ($threeprime > 0)
{
    if ($align[$threeprime]=~/\*/)
	{
		while ($primerlength < $min3prime && $align[$bp]=~/\*/)
		{
	    	$tempprimer[$primerlength] = $sequ[$bp];
	    	$primerlength++;
	    	$bp--;
		}
		if ($primerlength == $min3prime)
		{
	    	$mismatch = 0;
	    	while ($mismatch<$maxmismatch && $primerlength<$maxprimer)
			{
				if ($align[$bp]=~/\*/)
				{
		    		$tempprimer[$primerlength]= $sequ[$bp];
		    		$bp--;
		    		$primerlength++;
		    		if ($primerlength>$minprimer)
					{
						$primer[$primercount] = reverse (join "",@tempprimer);
						$length[$primercount] = $primerlength;
						$direction[$primercount] = "forward";
						$primerlocation[$primercount] = $bp + 1;
						$primercount++;
		    		}
				}
				else
				{
		    		$sequ[$bp]=~tr/AGCT/agct/;
		    		$tempprimer[$primerlength] = $sequ[$bp];
		    		$mismatch++;
		    		$blankcheck = $bp;
		    		$bp--;
		    		$primerlength++;
		    		if ($primerlength>$minprimer && $mismatch<$maxmismatch && $sequ[$blankcheck]!="-")
					{
						$primer[$primercount] = reverse (join "",@tempprimer);
						$length[$primercount] = $primerlength;
						$direction[$primercount] = "forward";
						$primerlocation[$primercount] = $blankcheck;
						$primercount++;
		    		}
				}
	    	}
		}
		else
		{
	    	$threeprime--;
	    	$bp = $threeprime;
	    	$primerlength = 0;
	    	@tempprimer = ();
		}
    }
	else
	{
		$threeprime--;
		$bp = $threeprime;
		$primerlength = 0;
		@tempprimer = ();
    }
}

# And now the reverse primers
$primerlength = 0;
$bp = $alignlength;
$threeprime = $alignlength;
$threeprime--;
while ($threeprime > 0)
{
    if ($rcalign[$threeprime]=~/\*/)
	{
		while ($primerlength < $min3prime && $rcalign[$bp]=~/\*/)
		{
	    	$tempprimer[$primerlength] = $rcsequ[$bp];
	    	$primerlength++;
	    	$bp--;
		}
		if ($primerlength == $min3prime)
		{
	    	$mismatch = 0;
	    	while ($mismatch<$maxmismatch && $primerlength<$maxprimer)
			{
				if ($rcalign[$bp]=~/\*/)
				{
		    		$tempprimer[$primerlength]= $rcsequ[$bp];
		    		$bp--;
		    		$primerlength++;
		    		if ($primerlength>$minprimer)
					{
						$primer[$primercount] = reverse (join "",@tempprimer);
						$length[$primercount] = $primerlength;
						$direction[$primercount] = "reverse";
						$primerlocation[$primercount] = $bp + 1;
						$primercount++;
		    		}
				}
				else
				{
		    		$rcsequ[$bp]=~tr/AGCT/agct/;
		    		$tempprimer[$primerlength] = $rcsequ[$bp];
		    		$mismatch++;
		    		$blankcheck = $bp;
		    		$bp--;
		    		$primerlength++;
		    		if ($primerlength>$minprimer && $mismatch<$maxmismatch && $rcsequ[$blankcheck]!="-")
					{
						$primer[$primercount] = reverse (join "",@tempprimer);
						$length[$primercount] = $primerlength;
						$direction[$primercount] = "reverse";
						$primerlocation[$primercount] = $blankcheck;
						$primercount++;
		    		}
				}
	    	}
		}
		else
		{
	    	$threeprime--;
	    	$bp = $threeprime;
	    	$primerlength = 0;
	    	@tempprimer = ();
		}
    }
	else
	{
		$threeprime--;
		$bp = $threeprime;
		$primerlength = 0;
		@tempprimer = ();
    }
}

# Save the sequences in a fasta format file as input to geecee
my $primerfasta = "tempprimer.txt";
open(PRIMERFILE, ">$primerfasta")
    or die "Could not create temporary primer file.\n";
$printcount = 0;
while ($printcount < $primercount)
{
    print PRIMERFILE ">$printcount\n$primer[$printcount]\n";
    $printcount++;
}
close PRIMERFILE;

# Calculate primer melting temps (Tm) and GC content and save them to an array
my $tmcount=0;
while ($tmcount < $primercount)
{
    my $x=0;
    my $g=0;
    my $c=0;
    my $a=0;
    my $t=0;
    my $sequence=$primer[$tmcount];
    $sequence =~ tr/agct/AGCT/;
    while ($x<$length[$tmcount])
	{
		my $nucleotide = chop($sequence);
		if ($nucleotide=~/G/)
		{
	    	$g++;
	    }
		else
		{
			if ($nucleotide=~/A/)
			{
		    	$a++;
			}
			else
			{
		    	if ($nucleotide=~/T/)
				{
					$t++;
		    	}
				else
				{
					if ($nucleotide=~/C/)
					{
			    		$c++;
					}
		    	}
			}
	    }
		$x++;
    }
    if ($length[$tmcount]<14)
	{
		$tm[$tmcount]= ((($a + $t)*2)+ (($c + $g)*4));
    }
	else
	{
		$tm[$tmcount]= (64.9 + (41 * ($g + $c - 16.4))/ ($a + $t + $g + $c));
    }
    $gc[$tmcount]= (($g + $c) / ($g + $c + $a + $t));
    $tmcount++;
}

# Verify the the Lab Primer Database can be opened
my $dbpresent=0;
my $database = "LabPrimerDatabase.txt";
unless (-e $database && -r $database)
{
    print "Lab Primer Database cannot be accessed\n";
    $dbpresent=1;
}

# Print the info on the primers within an acceptable range of Tm and GC and write sequences to files
my $goodprimer = "goodprimers.fasta";
open(GOODPRIMERFILE, ">$goodprimer")
    or die "Could not create good primer file.\n";
my $amplifx = "amplifxfile.txt";
open(AMPLIFXFILE, ">$amplifx")
    or die "Could not create amplifx primer file.\n";
$printcount = 0;
my $badprimer = 0;
print "\n";
printf "%-30s","Sequence";
printf "%-13s","Direction";
printf "%-9s", "Length";
printf "%-11s", "Location";
printf "%-13s", "GC Content";
printf "%-5s", "Tm";
print "  Equivalent";
print "\n";
printf "%-30s","--------";
printf "%-13s","---------";
printf "%-9s", "------";
printf "%-11s", "--------";
printf "%-13s", "----------";
printf "%-5s", "--";
print "  ----------";
while ($printcount < $primercount)
{
    if ($gc[$printcount]>$lowgc)
	{
		if ($gc[$printcount]<$highgc)
		{
	    	if( $tm[$printcount]>$lowTm)
			{
				if ($tm[$printcount]<$highTm)
				{
		    		print "\n";
		    		printf "%-30s",$primer[$printcount];
		    		printf "%-13s", $direction[$printcount];
		    		printf "%-9s", $length[$printcount];
		    		printf "%-11s", $primerlocation[$printcount];
		    		printf "%.3f        ", $gc[$printcount];
		    		printf "%.1f", $tm[$printcount];
		    		print GOODPRIMERFILE ">$printcount\n$primer[$printcount]\n"; #Print to a fasta file
		    		print AMPLIFXFILE ">$primer[$printcount]\t$printcount\t$direction[$printcount]-$primerlocation[$printcount]\n";    #Print in AmplifX file format
                    #check for equivalent primers in the database if it exists
		    		if ($dbpresent!=0)
					{
						open DBFILE, "$database";
						$headerline = <DBFILE>;
						while (my $dbprimerdata = <DBFILE>)
						{
			    			$dbprimerdata =~/(\S+)\s+(\d+)/;
			    			my $dbsequence = $1;
			    			my $dbreference = $2;
			    			if ($dbsequence=~/$primer[$printcount]$/)
							{
								print "  (ref#$dbreference $dbsequence)";
			    			}
							else
							{
								if ($primer[$printcount]=~/$dbsequence$/)
								{
				    				print "  (ref#$dbreference $dbsequence)";
								}
			    			}
						}
		    			close DBFILE;
		    		}
				}
				else
				{
		    		$badprimer++;
				}
	    	}
			else
			{
				$badprimer++;
	    	}
		}
		else
		{
	    	$badprimer++;
		}
    }
	else
	{
		$badprimer++;
    }
    $printcount++;
}
print "\nLowercase nucleotides are not conserved between taxa.\n";
print "$badprimer oligos were rejected for unacceptable GC content and/or Tm.\n\n";
close GOODPRIMERFILE;
close AMPLIFXFILE;

#Delete all temporary files
qx/rm tempalign.aln/;
qx/rm tempprimer.txt/;
qx/rm *.dnd/;


##################################################################################################
# Future improvements:
#   Read the name of the Labprimerdatabase from the command line
#   Read the various primer limitations in from the command line or use the hardcoded defaults
#   Include other measures of primer quality- selfdimers, polynucleotides, etc
#   Delete only the .dnd created during the script, not all .dnd files
#   Use strict
#   Step through arrays using foreach when possible
#   Do a pattern match for *******'s rather than stepping through the alignment bp by bp
#   Make the database comparisons case insensitive
#   Better, substitute all possible bps for lowercase letters when comparing to the database
##################################################################################################
