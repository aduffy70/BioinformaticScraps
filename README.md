# NAME
BioinformaticScraps

# DESCRIPTION
A semi-organized collection of my bits and pieces for processing bioinformatics data


# AFLP SCRIPTS:
pafreq.pl
 * This script summarizes the dominant allele frequencies for each population in an AFLP dataset, where 1 indicates presence of a band and 2 indicates lack of a band.  Sample names must be in the form taxon-population-individual.  Data must be sorted so each population is together.  Data must be whitespace delimited with a header row containing the names of the loci.  Missing data may cause problems.

ppfreq.pl
  * This script summarizes the phenotype frequencies for each population in an AFLP dataset, where 1 indicates presence of a band and 2 indicates lack of a band.  Sample names must be in the form taxon-population-individual.  Data must be sorted so each population is together.  Data must be whitespace delimited with a header row containing the names of the loci.  Missing data may cause problems.

tafreq.pl
  * This script summarizes the dominant allele frequencies for each population and taxon in an AFLP dataset, where 1 indicates presence of a band and 2 indicates lack of a band.  Sample names must be in the form taxon-population-individual.  Data must be sorted so each taxon is together.  Data must be whitespace delimited with a header row containing the names of the loci.  Missing data may cause problems.

tpfreq.pl
  * This script summarizes the phenotype frequencies for each population and taxon in an AFLP dataset, where 1 indicates presence of a band and 2 indicates lack of a band.  Sample names must be in the form taxon-population-individual.  Data must be sorted so each taxon is together.  Data must be whitespace delimited with a header row containing the names of the loci.  Missing data may cause problems.


# COALESCENCE SCRIPTS:
coalesce.pl
  * This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  Each individual randomly contributes 0, 1, or 2 offspring to the next generation with no limits on the size of the resulting population.  In each generation you can see which individuals are descended from which individual in the starting population.  Both the size of the population and the proportion descended from each starting individual should respond to drift.

coalesce2.pl
  * This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  There are no limits on the number of offspring any individual can have but the population size is fixed.  In each generation you can see which individuals are descended from which individual in the starting population.  The proportion descended from each starting individual should respond to drift.

coalesce3.pl
  * This script simulates the effects of drift on a population of asexual organisms with non-overlapping generations and no selection.  There are no limits on the number of offspring any individual can have but the population size is fixed.  This script runs multiple simulation replicates and tracks the number of generations to fixation in each replicate.  Then it plots a histogram of those values.

truecoalesce.pl
  * Unlike my previous "coalescent" scripts, this one actually models coalescence as opposed to drift and lineage sorting.  It calculates the time to coalescence and the time to get to the last 2 lineages, so you can see how much of the time was spent waiting for the last two lineages to coalesce.  You can select a sample size separate from the population size.


# ALIGNMENT/CONSENSUS SEQUENCE SCRIPTS
consaminocount.pl
  * This script takes a consensus sequence (with X's for ambiguous characters) in fasta format and counts the number of ambiguous characters (i.e., nonconserved characters) in each block of 50 aminoacids.

conscount.pl
  * This script takes a consensus sequence in fasta format with ambiguous characters and counts the number of ambiguous characters (i.e., nonconserved characters) in each block of 50 nucleotides.

muscletrans.py
  * This script takes a fasta file of dna sequences and performs a translation alignment using MUSCLE.  Stop codons in the translation are converted to unknown amino acids for the alignment (so MUSCLE doesn't drop them), but the original nucleotides are displayed in the final dna alignment.  The final alignment is saved to a file in interleaved phylip format.

assae.py
  * Aaron's Super Simple Alignment Editor.  My commandline-based multiple sequence alignment viewer/editor.


# BLAST SCRIPTS
countblastqueries.pl
  * This script finds the number of different blast queries that had hits in tabular blast output.  The list must be sorted by query.  Usage- I have a list of 7500 blast hits against my adiantum/angiopteris genes database.  I want to know how many of the queries in the original blast input file had at least one match in the database.

countgenes.pl
  * This script finds the number of different genes (sources) in tabular blast output.  The list must be sorted by gene (source) name.  Usage- I have a list of 7500 blast hits against my adiantum/angiopteris genes database.  I want to know how many different genes I have matches in.

msatblast.py
  * This script takes the blast output produced by conservedmsats.py along with the contig fasta files for two taxa and aligns each pair of sequences that showed blast similarity.  Then the alignments can be manually inspected to see if the microsatellite is actually within the aligned region and whether there is enough conserved sequence to design universal primers.

batchblast.py
  * This script takes a fasta file, blasts the sequences against NCBI's nr database and saves the result as an xml file.

blastglance.py
  * Returns the total number of queries and number of queries with no hits in a blast xml file.

microarray2go.py
  * Takes a csv file of microarray markers as input and returns a file showing the top blastn hit of the marker 70mers against an EST database and the top blastx hit of those ESTs against the nr protein database.  Those protein sequences can then be used in a searches against GO databases.

# FASTA FILE SCRIPTS
fastaclean.pl
  * This script takes messy fasta files (with numbers or whitespace) and converts them into clean fasta files

fastarename.pl
  * This script takes fasta files with Genbank names (gi|xxxxxxx|) and coverts the names to 4 letters for the genus and 4 letters for the species epithet.  This is a quick and dirty modification of my fastaclean script.

fastafilter.py
  * This script takes a file of fasta sequences and a file listing the sequences wanted and pulls just those sequences (or everything other than those sequences) from the fasta file.  A first simple script using biopython.

fastafilterq.py
  * This script takes a file of fasta sequences and a desired sequence id and extracts just that one sequence from the fasta file.

fastadiff.py
  * This script performs a bi-directional comparison of two fasta files and returns the id's of sequences in one but missing from the other.  Surprisingly quick for even large fasta files.

fastanames.pl
  * Just takes a fasta file and returns a list of the sequence id's.

fastaglance.py
  * Returns the total number of records, mean sequence length, longest sequence, and shortest sequence for a fasta file.

fastalengths.py
  * Returns the name and length of each sequence in a fasta file.

fastabreaks.py
  * Removes line breaks so that a fasta file has only 2 lines per record: 1 for the id and one for the sequence.  This avoids issues that extra line breaks can cause for some scripts.

sortbycsv.py
  * Divides a single file of fasta sequences into multiple files based on lines of a csv file of fasta descriptions.  All the items on a single row will go into a file.  The name of each file comes from the first item on the row.  Useful for sorting fasta records by gene.  Fasta descriptions with spaces are ok.

# GENOME FILE SCRIPTS
gb2tbl.py
  * This script converts a genbank flatfile to a features table for importing into Sequin and submitting to genbank.

transexcept.py
  * This script was used when we had a list of annotated RNA edits in gff format that needed to be inserted into the correct CDS entries of a features table as trans_except tags.

listgenes2.pl
  * This script filters ncbi genome files and returns just the gene names, locations, and orientations (1-forward, 0-reverse).  _WARNING_ This script will miss any exons beyond the first two.  Manually check for genes with 2 or more introns!

listgenes.py
  * Similar function to listgenes2.py but uses BioPython rather than raw perl parsing so it handles multiple exons better.

listgenesXXXX.py
  * A series of versions of listgenes.py optimized to handle idiosyncrasies of specific genbank genome records:

  * listgenesPsnu.py Psilotum nudum (NC_003386)-like genomes:
    * tRNAs have anticodons included in the /gene name
    * Look for gene names in /gene, then, /note, then /product to capture orf#'s in the /note field.
    * No /pseudo tags

  * listgenesAdca.py Adiantum capillus-veneris (NC_004766)-like genomes:
    * Some tRNAs have anticodon annotations in /note so we try /note first.  If that fails we can get the name without annotation from /gene but then we can't filter for duplicates.
    * For non-tRNAs look for gene names in /gene, then /note, then /product.
    * No /pseudo tags

  * listgenesCyta.py Cycas taitungensis (NC_009618)-like genomes:
    * tRNAs don't have anticodon names in /gene, but it is in /codon_recognized, so we need to read both and combine them
    * Look for gene names in /gene, then, /note, then /product to capture orf#'s in the /note field.
    * Pseudogenes don't have a CDS, tRNA, or rRNA section so we need to check the gene sections for a /pseudo tags

  * listgenesPhap.py Phalaenopsis aphrodite (NC_007499)-like genomes:
    * Same as listgenesCyta.py except for tRNAs we need to combine /product and /note instead of /gene and /note.

  * listgenesPith.py Pinus thunbergii (NC_001631)-like genomes:
    * Pseudogenes don't have a CDS, tRNA, or rRNA section so we need to check the gene sections for a /pseudo tags
    * tRNAs dont have anticodon annotations in ANY field so we can't filter duplicates tRNAs.
    * Look for gene names in /gene, then product, then /note.

listgenesfasta.pl
  * This script filters an ncbi genome file and a genome fasta file and returns a fasta file of the genes.  For genes with 2 exons and an intron, it returns the sequence of both exons with the intron, not separate entries for each exon.
  * Known Issues:
    * It behaves strangely for genes like ndhB that have one exon at the start of the genome sequence and the other at the end.  You can fix those manually, or cheat by modifying the genbank file.
    * Genes with more than two exons will also cause problems.  They won't appear in the output.
    * Do fasta files allow line breaks within the sequence?  If your fasta file has line breaks within the sequence, only the data up to the first line break will be captured.  You can use my fastaclean.pl script to remove extra line breaks from a fasta file.

structure2fasta.pl
  * This script takes the output of garnier (from the emboss package) and converts it to a fasta-like output so secondary structures can be aligned using alignment programs.  Hey- I'm not claiming this is really useful...

# MISC SCRIPTS
adjacentpairs.pl
  * This script determines which pairs of adjacent genes in the first file are also present in the second file.  The order of the genes in the pairs matters- they must be transcribed in the same order.  Both files must be formatted with two genenames in each row separated by white space.

findprimer.pl
  * This script:
    1. Performs a multiple nucleotide alignment using clustalw.
    2. Finds all primers with a length range with min conserved 3' and max# of nonconserved nucleotides.
    3. Filters out primers that don't fall within specified Tm and GC content ranges.
    4. Searches against the Lab Database for equivalent primers.
    5. Outputs to standard output, a fasta format file and a amplifx primer list format file.
  * Input file requirements:
    1. fasta file with at least two sequences.
    2. text file called LabPrimerDatabase.txt containing a primer database in the format: sequence  reference  any_other  additional_columns
  * This script generates the following files, clobbers them if they already exist, and removes the "temp" ones after use:
    * tempalign.aln  
    * tempprimer.txt  
    * temprev.txt  
    * tempgc.txt   
    * goodprimers.txt
    * amplifxfile.txt
  * clustalw must be in your PATH:

pop2struct.pl
  * This script converts a genepop input file to a structure input file.  It only works on files with a specific number of loci and the same number of digits in each allele record.

gettaxanames.pl
  * This script extracts the taxon names from a nexus file.

unique.pl
  * This script finds the unique items in 2 files containing a list of genes.  The files must be lists with a single gene name starting each row (the output from my listgenes.pl script for example).  Any other information on the line will be ignored.  This script does not compare the files very efficiently so it may take a while for very long files.

conservedmsats.py
  * This script finds microsatellites that may be conserved between two taxa.  It takes msatcommander output files and contig fasta files for two taxa and finds contigs containing the same repeat type that also show blast similarity (in areas other than the repeat, of course).  If we only have a sampling of two genomes this script may help us identify microsatellites that are homologous.

mlcminer.py
  * This script extracts dN/dS, dN and dS values from a directory of PAML output files (.mlc files).  It "mines" data from mlc files.  It is 3 separate scripts in one file.  Uncomment one section at a time to get dNdS values, dN values, and dS values.

mlcminer2.py
  * This script extracts dN and dS values from a directory of PAML output files (.mlc files) and writes them to a separate csv file for each gene.  It "mines" data from mlc files.  Unlike mlcminer.py, this script doesn't average the dN and dS values - it extracts the raw values so we can use them in statistical analyses.  It also ignores the outgroup values.

wilcox.r
  * These [R] commands run two wilcoxon rank tests (1 on the first 2 columns of a csv file and one on the 3rd & 4th columns).

wilcox.py
  * This script is a python wrapper for the wilcox.r script.  That script runs two wilcoxon rank tests (1 on the first 2 columns of a csv file and one on the 3rd & 4th columns).  This wrapper lets us run the R script on every csv file in a folder and save the results to one long text file.

wilcoxparse.py
  * This script extracts just the p-values from the output of wilcox.py.  Output of this script is csv format: GeneName,dN_p-value,dS_p-value

compare2lists.py
  * This script compares two files of genenames (or anything, really) and returns counts and lists of items in just the first file, just the second file, or both files (i.e., the data needed to generate a venn diagram).

comparelists.py
  * This script does the same thing as compare2lists.py but works on any number of files, returning the counts and lists of items in just one file or in every combination of files.

extract_hplc_data.py
  * Extracts data from multiple HPLC files in the format needed for analysis in the Coulombe Lab.

filter_structure_file_by_informative_loci.py
  * Removes non-variable and autapomorphic loci from a structure file
