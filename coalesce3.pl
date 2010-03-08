#! /usr/bin/perl
# coalesce3.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created November 2007

# This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  There are no limits on the number of offspring any individual can have but the population size is fixed.  This script runs multiple simulation replicates and tracks the number of generations to fixation in each replicate.  Then it plots a histogram of those values.

print "Population size? ";
$pop = <>;
$pop--;
print "Number of replicates? ";
my $repcount = 0;
my $replicates = <>;
while ($repcount < $replicates)
{
	$popsize = $pop+2;
	$quit=0;
	$generation=0;
	$count=1;
	while ($count < $popsize)
    {
		push @population, $count;
		$count++;
	}
	$popsize = @population + 0;
#	print "\n\n\nStarting generation:\n";
#	print "Gen=$generation     @population\n";
	while ($quit==0)
    {
		foreach $individual (@population)
        {
			$offspring = int(rand ($popsize));
			$count=0;
			push @nextgeneration, $population[$offspring];
		}
		$generation++;
		$popsize = @nextgeneration +0;
		splice(@population);
		sub numerically { $a <=> $b }
		@population = sort numerically @nextgeneration;
#		print "Gen=$generation     @population\n";
		splice(@nextgeneration);
		if ($population[0] == $population[$pop])
        {
			$quit = 1;
		}
	}
	$repcount++;
	splice(@population);
	push @summary, $generation;
}
@summary = sort numerically @summary;
#print"\nNumber of generations for each replicate to sort to a single lineage:\n@summary\n\n\nHistogram:\n";

#Make a histogram of the generation numbers
$biggest = $replicates - 1;
$binsize = int($summary[$biggest]/30);
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
print "\n\n";
$sum=0;

#Calculate the mean number of generations
foreach $value (@summary)
{
	$sum = $sum + $value;
}
$mean = $sum / $replicates;
print "Mean = $mean\n\n";
