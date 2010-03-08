#! /usr/bin/perl
# coalesce2.pl
# Aaron M Duffy aduffy70{at}gmail.com
# Created November 2007

# This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  There are no limits on the number of offspring any individual can have but the population size is fixed.  In each generation you can see which individuals are descended from which individual in the starting population.  The proportion descended from each starting individual should respond to drift.


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
		$offspring = int(rand ($popsize));
#		print "$offspring\n";
		$count=0;
		push @nextgeneration, $population[$offspring];
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









