#! /usr/bin/perl
#truecoalesce.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created November 2007

# Unlike my previous "coalescent" scripts, this one actually models coalescence as opposed to drift and lineage sorting.  It calculates the time to coalescence and the time to get to the last 2 lineages, so you can see how much of the time was spent waiting for the last two lineages to coalesce.  You can select a sample size separate from the population size.

print "Population size (# of gene copies)? ";
$N = <>;
print "Sample size (# of gene copies)? ";
$sample = <>;
$popsize = $sample;
print "Replicates? ";
$replicates = <>;
print "Display histogram of coalescence times (y/n)? ";
$hist = <>;
$repcount=0;
while ($repcount<$replicates)
{
	my $generation=0;
	$count=1;
	my @parentpop;
	$already2=0;
	while ($popsize>1)
	{
		$count=1;
		while ($count <= $popsize)
		{
			$parent=int(rand($N));
		    undef @is_a_parent;
		    for (@parentpop)
			{
				$is_a_parent[$_] = 1;
			}
			if ($is_a_parent[$parent])
			{
				#print "$parent match @parentpop\n";
			}
			else
			{
				#print "$parent no match @parentpop\n";
				push @parentpop, $parent;
			}
			$count++;
		}
		$popsize = @parentpop+0;
		if (($popsize == 2) && ($already2==0))
		{
			push @last2, $generation;
			$already2=1;
		}
		$generation++;
		splice(@parentpop);
		splice(@is_a_parent);
	}
	$repcount++;
	$popsize=$sample;
	push @summary, $generation;
}
sub numerically { $a <=> $b }
@summary = sort numerically @summary;
#print "@summary\n\n\n\n @last2\n";

#Make a histogram of the generation numbers
if ($hist=~/y/)
{
	$biggest = $replicates - 1;
	$binsize = int($summary[$biggest]/50);
	$binb=1;
	$bint=$binsize;
	while ($binb < $summary[$biggest])
	{
		printf "\n%-6sto%6s ", $binb, $bint;
		foreach $value (@summary)
		{
			if ($value >= $binb)
			{
				if ($value <= $bint)
				{
					print "*";
				}
			}
		}
		$binb = $binb + $binsize;
		$bint = $bint + $binsize;
	}
	print "\n";
}

#Calculate the mean number of generations
$sum=0;
foreach $value (@summary)
{
	$sum = $sum + $value;
}
$mean = $sum / $replicates;
print "\nMean coalescent time= $mean\n";
$sum2=0;
foreach $value2 (@last2)
{
	$sum2 = $sum2 + $value2;
}
$mean2 = $sum2 / $replicates;
print "Mean time to coalesce to last 2 lineages= $mean2\n\n";
