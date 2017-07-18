#!/sw/bin/perl

### Virag Sharma, October 2013.
### Script to remove "e lines" for those regions which overlap with ExcludeRegions (Exclude Regions are those that are not covered by any chain or overlap with assembly gaps).  
### Please note that ExcludeRegions can only be defined for a pair for example hg19-mm9. The coordinates of these regions are with respect to the reference i.e. hg19 in this case.

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;
use MyKentFunctions;

my $genomePath = $ENV{'genomePath'};
die "ERROR: directory specified in environment variable genomePath (set to $genomePath) does not exist\n" if (! -d $genomePath);

$| = 1;		# == fflush(stdout)

# input options/parameters
our $verbose = 0;
my $MafIn   = ""; 	## the input MAF file.
my $MafOut  = "";   # the output file.
my $reference = "";	# reference genome
my $alignment = "";	# the alignment , e.g hg19_46way, mm10_60way
my $excludeRegionsDir = "";	## This is optional, if I set it to some directory, then these exclude files will be used, NOT the standard

sub usage 
{
	die "usage: $0 -input InputMaf -output OutputMaf -reference ReferenceSpecies -alignment Alignment (hg19_46way, mm10_60way etc) -excludeRegionsDir (Optional: use this parameter if you do not want to use the DEFAULT exclude regions)\n";
}

GetOptions ("v|verbose"  => \$verbose, "input=s" => \$MafIn, "output=s" => \$MafOut, "reference=s" => \$reference, "alignment=s" => \$alignment, "excludeRegionsDir:s" => \$excludeRegionsDir) || usage();

my $excludeDir =  "$genomePath/data/PercentIDPipeline/ExcludeRegions/$alignment/";
$excludeDir = $excludeRegionsDir if ($excludeRegionsDir);	## This is the directory now.

usage() if ($MafIn eq "" || $MafOut eq ""); ### die if either of the two- input or output file is not specified
usage () if ( ($reference eq "") || (! -d $excludeDir) ); ### die if no reference species is specified or the exclude directory does not exist.

print "Exclude Regions are being read from '$excludeDir'\n" if ($verbose);

## Just check if the file has nothing but the header lines and no maf blocks (this is the case with a few chromosomes like: chrUn_gl000226.maf -hg19_46way/100way )

my $nol = `wc -l $MafIn|awk '{ print \$1}'|tr -d "\n"`; ## Get the number of lines in the inputMaf.
if ($nol == 1)
{
	system("cp $MafIn $MafOut"); ## Copy the input file to the output file directory.
	exit; ## And just exit. Nothing else needs to be done.
}

my $Ref2SpeciesList = GetSpeciesList($MafIn);
my @species_array = @$Ref2SpeciesList;
@species_array    = grep {$_ ne $reference} @species_array;		## Remove the ReferenceSpecies from the list of species because there of course is no ReferenceSp.ReferenceSp.EXCLUDE.bed
		
#### Now check if there are Exclude region files for the species/assemblies present in @species_array.		

my @SpExcludeUA = ();		## SpeciesExcludeUnAvailable
my $ct = 0;
print "Now checking if all the species have an exclude file\n" if ($verbose);

foreach my $species (@species_array)									
{	
	my $ExFile = "$excludeDir/$reference.$species.EXCLUDE.bed";			## These are the exclude files.
	
	if (! -e $ExFile)													## Check if the exclude file is present
	{
		$SpExcludeUA[$ct] = $species;
		$ct++;		
	}
}

if  ($ct != 0)
{
	my $ExcludeUA = join(",",@SpExcludeUA);
	die "Exclude files are not present for $ExcludeUA\n" if (! $excludeRegionsDir);			### Will die at this point if the DEFAULT exclude files are missing for any assembly/species.
}																							### If I specify a directory and some of the exclude files are missing, then is it not a problem

print "All exclude files are present. Now ready to remove e lines\n";

#### Now create a bed file which contains the list of coordinates from the reference assembly (for the input maf) since the coordinates in the exclude file are also wrt to the Reference assembly (e.g. hg19).

my $referenceBed = create_temp("f");

print "Reading the $MafIn and extracting the coordinates of each aligning block in the reference genome .. \n" if ($verbose);
open(FI,"$MafIn") || die "Error opening input file which is $MafIn\n";
open(FO,">$referenceBed");

while (my $line=<FI>)
{
	$line=~s/\s+$//;	# remove potential whitespace at the end of a line
	
	# focus only on sequence lines for the reference species
	if ( ($line =~/^s/) && ($line =~/$reference/) )
	{		
		my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line); 
		my $name = "$species\_$chr\_$start\_$size\_$strand";		#### The name is the "unique key" here

		my $stop = $start + $size;
		print FO "$chr\t$start\t$stop\t$name\n";					### The fourth field in the bed file i.e. $name helps us to keep a track of each line later.
	}
}
close FI;
close FO;

my %ExcludeHash = ();		## This is a 2 dimensional hash and stores information about those regions in the reference which overlap with the exclude regions
							## The first dimension is the species, the second dimension is the coordinates of the reference.
							
print "Now running bedtools intersect to see if there is an overlap with the exclude regions\n";
foreach my $species(@species_array)
{	
	my $SpeciesExcludeFile = "$excludeDir/$reference.$species.EXCLUDE.bed";
	my $overlapFile = create_temp("f");

	print "This is the file --> $SpeciesExcludeFile ..";
	if (-e $SpeciesExcludeFile)
	{
		print " And yes, it is there";
		my $status = system("bedtools intersect -a $referenceBed -b $SpeciesExcludeFile -wa > $overlapFile"); 
		## The -wa option does the following -> If the element in referenceBed is chr1:100-200 
		## and an element in the $SpeciesExcludeFile is chr1:150-200, what is reported as an overlapping element is chr1:100-200
		
		die "Error running bedtools intersect for $species" if ($status != 0);
		print "Now Reading overlapFile for $species ..." if ($verbose);	
		my $numberOverlap = 0;

		open(FIO,"$overlapFile") || die "Error opening overlap file which is $overlapFile\n";	
		while (my $line = <FIO>)
		{
			$line =~s/\s+$//;
			my @tmp = split(/\t/,$line);
			my $key = $tmp[3];
		
			$ExcludeHash{$species}{$key} = "T";
			$numberOverlap++;
		}
		close FIO;

		print " --> Number of overlapping elements is $numberOverlap\n" if ($verbose);
		
		system("rm $overlapFile");			## Delete these temporary files, no point in filling the tmp directory.
	}
	print "\n";
}
print "DONE\n" if ($verbose);

system("rm $referenceBed");					## Delete the bed file for the reference as well.

## At this stage, read the entire file again ($MafIn) and remove the elines where the reference genome coordinates overlap with the ExcludeRegions.

print "Now reading $MafIn again and removing e lines which overlap with the exclude regions ...\n";
open(FIM,$MafIn) || die "Error opening chromosome maf file\n";
open(FO,">$MafOut") || die "Error opening output file which is $MafOut \n";

my $id_old = "";

while (my $line=<FIM>)
{	
	$line =~s/\s+$//;	# remove potential whitespace at the end of a line
	
	if ( ($line =~/^s/) and ($line =~/$reference/) )
	{		
		my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line); 
		
		$id_old="$species\_$chr\_$start\_$size\_$strand";		## this is the key that is used in the hash (ExcludeHash) created earlier.
	}
	
	if ($line =~/^e/) 												## This is the check point, checks if the line is an "e" line.
	{
		my ($species, $dummy1, $dummy2, $dummy3, $dummy4, $dummy5, $dummy6) = getELineData($line);	## 	Get the species for this e line. dummy1-6 are variables that we do not care about.
			
		if (! $ExcludeHash{$species}{$id_old}) 						## If this line is not present in the ExcludeHash, print it.
		{	
			print FO "$line\n";
		}
		else	
		{ 
			print " --> $id_old and $species ..  $ExcludeHash{$species}{$id_old} \n" if ($verbose);
		}
	}
	
	else { print FO "$line\n"; }										## Keep printing, no need to do any check for anything but the "e" line.
}

close FIM;
close FO;

print "e lines successfully removed from $MafIn and new maf is written to $MafOut\n";


#### Sub routines follow.

sub GetSpeciesList   ### This function returns a reference to an array which contains all species present in the maf.
{
	my $MafIn = shift;		### The only argument is the MAF block.
	my @species_array = ();
	@species_array = `mafSpeciesList $MafIn stdout`; 
	chomp(@species_array); 
	
	die "ERROR: mafSpeciesList must have failed because the list of species is empty\n" if (scalar @species_array < 1); 
	$" = " ";	### separate array elements by a space in the output
	
	print "List of species in the maf: @species_array\n" if ($verbose);
	return (\@species_array);
}
