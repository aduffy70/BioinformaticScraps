#! /usr/bin/perl
# coalesce.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created October 2007

# This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  Each individual has randomly contributes 0, 1, or 2 offspring to the next generation with no limits on the size of the resulting population.  In each generation you can see which individuals are descended from which individual in the starting population.  Both the size of the population and the proportion descended from each starting individual should respond to drift.


my $quit=0;
my $generation=0;
print "\n\n\nPopulation size?  N=";
my $popsize = <>;
$popsize++;
my $count=1;
while ($count < $popsize)
{
	push @population, $count;
	$count++;
}
$popsize = @population + 0;
print "\n\n\nStarting generation:\n";
print "Gen=$generation     N=$popsize     @population\n\n";
print "After each generation, press Enter to continue or 'q' to quit\n\n";
my $continue = <>;
while ($quit==0)
{
	foreach $individual (@population)
    {
#		print "$individual ";
		$offspring = int(rand 3);
#		print "$offspring\n";
		$count=0;
		while ($count < $offspring)
        {
			push @nextgeneration, $individual;
			$count++;
		}
	}
	$generation++;
	$popsize = @nextgeneration +0;

	splice(@population);
	sub numerically { $a <=> $b }
	@population = sort numerically @nextgeneration;
	print "Gen=$generation     N=$popsize     @population\n";
	splice(@nextgeneration);
	$continue = <>;
	if ($continue=~/q/)
    {
		$quit=1;
	}
}









