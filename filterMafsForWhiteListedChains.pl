#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany
## This script filter maf files and keeps only those "s" and "i" lines which come from the whiteListed chains plus those that DO not come from blackListed chains
## Additionally, it also filters out 1) e-lines which lie in the middle of the blacklisted s-lines
## 2) e-lines which lie after non-whiteListed s-lines but the s-lines after the e-lines block are whiteListed because they come from a different scaffold
## 3) e-lines which lie before non-whiteListed s-lines but the s-lines before the e-lines block are whiteListed because they come from a different scaffold

use strict;
use warnings;

use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;
use MyKentFunctions;

sub usage
{
	die "usage: $0 [-v|verbose] -in Input_file -out Output_file -ref ReferenceSpecies -whiteListedChainsDir WhiteListedChainsDirectory -blackListedChainsDir BlackListedChainsDirectory\n";
}

my $verbose = 0;
my $in = my $out = my $ref = my $wlDir = my $blDir = "";
GetOptions ("v|verbose"  => \$verbose, "in=s" => \$in, "out=s" => \$out, "ref=s" => \$ref, "whiteListedChainsDir=s" => \$wlDir, "blackListedChainsDir=s" => \$blDir) || usage();

usage() if ($in eq "");

## Check if all the input parameters have been specified properly
die "ERROR !! No output file specified\n" if ($out eq "");
die "ERROR !! No reference species specified\n" if ($ref eq "");
die "ERROR !! WhiteListed|BlackListed chains directory is not specified|does not exist\n" if ($wlDir eq "" || ! -d $wlDir || $blDir eq "" || ! -d $blDir);

## Create the output directory in case it does not exist
my $dirName = `dirname $out|tr -d "\n"`;
`mkdir -p $dirName`;

## Check if the file has nothing but the header lines and no maf blocks (this is the case with a few chromosomes like: chrUn_gl000226.maf ->hg19_46way/100way)
my $nol = `wc -l $in|awk '{ print \$1}'|tr -d "\n"`; ## Get the number of lines in the inputMaf.
if ($nol == 1)		## If yes, then copy this file to the output directory
{
	`cp $in $out`;
	exit; ## And just exit.
}

## Get a list of species present in the maf
my @speciesList = ();			### create a species_array which contains the list of all the species present in the maf block.
@speciesList = `mafSpeciesList $in stdout`; 
chomp(@speciesList); 

die "ERROR: mafSpeciesList must have failed because the list of species is empty\n" if (scalar (@speciesList) < 1); 
$" = " ";	# separate array elements by a space in the output
print "List of species in the maf: @speciesList\n" if ($verbose);
@speciesList = grep {$_ ne $ref} @speciesList;	### Remove the reference species from the speciesList

### Step 1: Now check if all whiteListed|blackListed chain file directories are present
foreach my $species(@speciesList)
{
	my $whiteListedChainsDir = "$wlDir/$species";
	my $blackListedChainsDir = "$blDir/$species";
	
	die "ERROR!!! The whiteListed chains directory i.e. $whiteListedChainsDir does not exist\n" if (! -d $whiteListedChainsDir);
	die "ERROR!!! The blackListed chains directory i.e. $blackListedChainsDir does not exist\n" if (! -d $blackListedChainsDir);
}
print "SUCCESS !!! All whiteListed and blackListed chain directories are present\n";

my @tempFilesList = ();  ## This array will contain all the tempFiles

my %sp2File = my %spBedFile = ();
foreach my $species(@speciesList)
{
	my $fh = FileHandle->new();        # need to create a new filehandle.
	$sp2File{$species} = $fh;

	my $bedFileMaf = `mktemp /dev/shm/XXXXXXX.bed`; chomp $bedFileMaf;
	push(@tempFilesList,$bedFileMaf);
	
	$spBedFile{$species} = $bedFileMaf;
	$fh->open(">$bedFileMaf") || safeDie("ERROR: cannot write to $bedFileMaf\n");
}

## Now start reading the maf and print the coordinates to the appropriate bed file.
print "Now beginning to read the mafFile\n";

my $chr = my $refString = "";
open(FI,"$in") || safeDie("Error opening input file '$in'\n");
while (my $line = <FI>)
{
	if ($line =~/^s/) 
	{		
		$line =~s/\s+$//;		## Remove white space at the end
		
		if ($line =~/$ref/)		## If this line comes from the reference species
		{
			my ($speciesR,$chrR,$startR,$sizeR,$strandR,$srcSizeR,$seqR) = getSLineData($line); 		## The "R" at the end of each variable implies that we are talking about the Reference.
			my $stopR = $startR + $sizeR;
			
			$refString = "$chrR\t$startR\t$stopR";
			$chr = $chrR;		## The chromosome is simply the chromosome of the reference species
		}
		else
		{
			my ($speciesQ,$chrQ,$startQ,$sizeQ,$strandQ,$srcSizeQ,$seqQ) = getSLineData($line); 		## The "Q" at the end of each variable implies that we are talking about the Query.
			
			my $fh = $sp2File{$speciesQ};
			my $stopQ = "";
			
			if ($strandQ eq "-")
			{
				$stopQ = $srcSizeQ - $startQ;
				$startQ = $stopQ - $sizeQ;
			}
			else
			{
				$stopQ = $startQ + $sizeQ;	
			}
			
			my $queryString = "$chrQ:$startQ:$stopQ:$strandQ";
			print $fh "$refString\t$queryString\n";
		}
	}
}
close FI;

print "First pass of the input maf $in over\nTime for overlap operations\n";

my %allSpeciesHashBlack = my %allSpeciesHashWhite =  ();	## This is a two dimensional hash. The first dimension is the species, the second is the mafRegion (reference species coordinates)
## If this region overlaps with a region that is covered by a blackListed chain in the query species, then this region is dumped into the hash. The value is "T"

foreach my $species(@speciesList)
{
	my $fh =  $sp2File{$species};		## Close previously opened file handles
	$fh->close;

## Now run overlapSelect to see which loci in the maf overlap with the regions covered with whiteListed chains (we work with reference genome here)
## and which ones are covered by blackListed chains
	my $whiteListedChainFile = "$wlDir/$species/$chr.bed";
	my $bedFile = $spBedFile{$species};
	my $tmpOverlap = `mktemp /dev/shm/XXXXXXX.bed`; chomp $tmpOverlap;
	push(@tempFilesList,$tmpOverlap);
	
	if (-e $whiteListedChainFile)
	{	
		print " --> $bedFile ($species) and whiteListedChains file '$whiteListedChainFile' are being compared\n" if ($verbose);
		
		## Sort the bedFile
		my $tmpBed = `mktemp /dev/shm/XXXXXXX.bed`; chomp $tmpBed;
		my $callSorting = "sort -k1,1V -k2,2n $bedFile > $tmpBed";
		system($callSorting) == 0 || safeDie("Error running sorting call\n");
		`mv $tmpBed $bedFile`;
		
		## Now check for overlap using bedtools intersect
		my $call = "bedtools intersect -a $bedFile -b $whiteListedChainFile -wa -wb -sorted > $tmpOverlap";
		system($call) == 0 || safeDie("ERROR running bedtools intersect for '$species'\n");
		
		my $ref2Hash = processOverlapFile($tmpOverlap,$species);		
		print "Overlapping maf coordinates with whiteListed chain coordinates for '$species' DONE\n";
		my %speciesHash = %$ref2Hash;
		
		foreach my $keys( keys(%speciesHash))
		{
			my $mafLocus = $keys;
			$allSpeciesHashWhite{$species}{$mafLocus} = "T";  ## This hash stores everything that overlaps with the whiteListed chains and should be retained
		}
	}
	else
	{
		print "WARNING !! The whiteListed chains bed file does not exist for this $ref-$species $chr pair. So nothing will overlap with whiteListed chains here \n\n" if ($verbose);
	}
	
	## Now doing the same operation for blackListed chains
	my $blackListedChainFile = "$blDir/$species/$chr.bed";
	
	if (-e $blackListedChainFile)
	{	
		print " --> $bedFile ($species) and blackListedChains file '$whiteListedChainFile' are being compared\n" if ($verbose);
		
		## Sort the bedFile
		my $tmpBed = `mktemp /dev/shm/XXXXXXX.bed`; chomp $tmpBed;
		my $callSorting = "sort -k1,1V -k2,2n $bedFile > $tmpBed";
		system($callSorting) == 0 || safeDie("Error running sorting call\n");
		`mv $tmpBed $bedFile`;
		
		## Now check for overlap using bedtools intersect
		my $call = "bedtools intersect -a $bedFile -b $blackListedChainFile -wa -wb -sorted > $tmpOverlap";
		system($call) == 0 || safeDie("ERROR running bedtools intersect for '$species'\n");
		
		my $ref2Hash = processOverlapFile($tmpOverlap,$species);
		print "Overlapping maf coordinates with blackListed chain coordinates for '$species' DONE\n";
		my %speciesHash = %$ref2Hash;
		
		foreach my $keys( keys(%speciesHash))
		{
			my $mafLocus = $keys;
			$allSpeciesHashBlack{$species}{$mafLocus} = "T";  ## This hash stores everything that overlaps with the whiteListed chains and should be retained
		}
	}
	else
	{
		print "WARNING !! The blackListed chains bed file does not exist for this $ref-$species $chr pair. So nothing will overlap with blackListed chains here\n\n" if ($verbose);
	}
	
	`rm -rf $bedFile $tmpOverlap`;
}

print "All the lines that should be thrown away are identified\nNow filtering out the sLines\n";

my %blackListedSLines = (); ## This is a 2 dimensional hash that stores the 's' lines which were filtered out. The first index is the species, the second is the maf block number.
my %eLines	      = (); ## This is again a 2 dimensional hash that stores 'e' lines. The first index is the species, the second is the maf block number.

## Now read the maf for the final time and print only those regions which DO NOT overlap with the blacklisted chains

my $flagILinePrint = "";
my $tmpInterMediate = create_temp("f");
push(@tempFilesList,$tmpInterMediate);

my $blockNumber = 0;
my %blockNumberQuerySpeciesScaffold = ();

open(FI,"$in") || safeDie("Cannot open the input maf --> $in\n");
open(FO,">$tmpInterMediate") || safeDie("Cannot write to the intermediate file --> $tmpInterMediate\n");
while (my $line = <FI>)
{
	$line =~s/\s+$//;

	if ($line =~/^s/) 
	{		
		$line =~s/\s+$//;
		
		if ($line =~/$ref/)		## If this line comes from the reference species
		{
			my ($speciesR,$chrR,$startR,$sizeR,$strandR,$srcSizeR,$seqR) = getSLineData($line); 
			my $stopR = $startR + $sizeR;
			
			$refString = "$chrR\t$startR\t$stopR";
			print FO "$line\n";
			$blockNumber++;
		}
		else
		{
			my ($speciesQ,$chrQ,$startQ,$sizeQ,$strandQ,$srcSizeQ,$seqQ) = getSLineData($line); 
			
			if (exists $allSpeciesHashWhite{$speciesQ}{$refString})  ## everything that overlaps with whiteListed chains is retained
			{
				print "'$refString' is retained for species '$speciesQ' because it comes from whiteListed chain\n" if ($verbose);
				print FO "$line\n";
				$flagILinePrint = "T";
			}
			elsif (! exists $allSpeciesHashBlack{$speciesQ}{$refString}) ## everything that does not overlap with whiteListed chains but does not overlap with
			{															 ## blackListed chains either is also retained
				print "'$refString' is retained for species '$speciesQ' because it DOES NOT from whiteListed chain but it DOES NOT come from a blackListed chain either\n" if ($verbose);
				print FO "$line\n";
				$flagILinePrint = "T";
			}
			elsif (exists $allSpeciesHashBlack{$speciesQ}{$refString})   ## everything that does not overlap with whiteListed chains and overlaps with blackListed chains
			{															 ## instead is retained.
				print "'$refString' excluded for species '$speciesQ' because it comes from blackListed chain\n" if ($verbose);
				$blackListedSLines{$speciesQ}{$blockNumber} = "T";
				$blockNumberQuerySpeciesScaffold{$speciesQ}{$blockNumber} = $chrQ;
			}
		}
	}
	
	elsif ($line =~/^i/)		## The "i" lines need an approval. If the 's' line is blackListed, then the corresponding 'i' line goes away as well.
	{
		if ($flagILinePrint eq "T")
		{
			print FO "$line\n";
			$flagILinePrint = "";
		}
	}
	
	else		### Print everything else, for e-lines, just get the species name and put it in the eLinesHash
	{
		if ($line =~/^e/)
		{
			my ($species,$dm1,$dm2,$dm3,$dm4,$dm5,$dm6) = getELineData($line);      ## All I care about is the species
			$eLines{$species}{$blockNumber} = "T"
		}
		print FO "$line\n";

	}
}
close FI;
close FO;

my %hashExcludeELines = (); ## This is a 2 dimensional hash, the first dimension is the species, the second is the blockNumber (from the reference) which corresponds to an e-Line (for a query) which should not be printed.
#### Now merge e-line blocks for each species. For instance, if the block number 6,7,8 in a maf have an e-line for a species (such as mm9), we make a composite block 6-7-8

foreach my $species (keys (%eLines))
{
	my @tmpArray = (keys (%{$eLines{$species}}));	## This gives me an array which contains all the e-Line numbers for $species

	my $ref2mergedBlocks = clubELines(\@tmpArray);	## Now I merge these e-Lines, see the function for more details
	my %mergedBlocks = %$ref2mergedBlocks;		## For the above example, the hash is %hash = (6 => 3)

	foreach my $eLineBlockStart(keys (%mergedBlocks))
	{
		my $start = $eLineBlockStart;				## start is 6
		my $stop  = $start + $mergedBlocks{$eLineBlockStart} -1;	## stop is 8, mind the "-1" at the end.

		my $flagExclude = "";
		my $scaffoldBlackListedChain = my $scaffoldChainNotBlackListed = "";
		
		if ( (exists $blackListedSLines{$species}{$start -1}) &&  (exists $blackListedSLines{$species}{$stop +1}) )	## check if the block upstream i.e 5
		{														## and block downstream i.e 9 is blacklisted. If yes, then blackList the elines from 6 till 8
			$flagExclude = "T";
		}
		
		elsif (exists $blackListedSLines{$species}{$start -1})		## This could happen when you have some e-lines at the end of a chain. Only the upstream s-lines are blackListed
		{															## The downstream s-lines are not blacklisted, but they come from a different scaffold	
			$scaffoldBlackListedChain 	= $blockNumberQuerySpeciesScaffold{$species}{$start -1} if (exists $blockNumberQuerySpeciesScaffold{$species}{$start -1});
			$scaffoldChainNotBlackListed    = $blockNumberQuerySpeciesScaffold{$species}{$stop +1}  if (exists $blockNumberQuerySpeciesScaffold{$species}{$stop +1});
			
			$flagExclude = "T" if ($scaffoldBlackListedChain ne $scaffoldChainNotBlackListed);	## i.e the two are different
		}
		
		elsif (exists $blackListedSLines{$species}{$stop +1})		## This could happen when you have some e-lines at the beginning of a chain. Only the downstream s-lines are blacklisted
		{															## The upstream s-lines are not blackListed but they come from a different scaffold
			$scaffoldBlackListedChain 	= $blockNumberQuerySpeciesScaffold{$species}{$stop +1}  if (exists $blockNumberQuerySpeciesScaffold{$species}{$stop +1});
			$scaffoldChainNotBlackListed 	= $blockNumberQuerySpeciesScaffold{$species}{$start -1} if (exists $blockNumberQuerySpeciesScaffold{$species}{$start -1});
			
			$flagExclude = "T" if ($scaffoldBlackListedChain ne $scaffoldChainNotBlackListed);
		}
		
		if ($flagExclude eq "T")
		{		
			while ($start <= $stop)			
			{
				$hashExcludeELines{$species}{$start} = "T";
				$start++;
			}
		}
		
	}
}

### Now the third and the final pass
print "Now writing to the output file and filtering out the blackListed eLines at the same time --> $out\n";

open(FI,$tmpInterMediate) || safeDie("FATAL !!! Cannot open the intermediate file '$tmpInterMediate' for reading\n");
open(FO,">$out") || safeDie("FATAL !!! Cannot open the output file '$out' for writing\n");

$blockNumber = 0;

while (my $line = <FI>)
{
	my $flagPrint = "T";

	$blockNumber++ if ( ($line =~/$ref/) && ($line =~/^s/) );		
	
	if ($line =~/^e/)
	{
		my ($species,$dm1,$dm2,$dm3,$dm4,$dm5,$dm6) = getELineData($line);	## All I care about is the species
		$flagPrint = "F" if (exists $hashExcludeELines{$species}{$blockNumber})			
	}

	print FO "$line" if ($flagPrint eq "T")		;
}

close FI;
close FO;

#############################################
##### End of the main body of the code ######
##### Functions begin		       ######
#############################################

sub processOverlapFile
{
	my ($overlapFile,$species) = @_;
		
## The file looks like the following:
## chr1	44440712	44440779	scaffold_8422:9237:9271:-	chr1	44440655	44440714	308663_scaffold_8422:9269-9328_-
## chr1	44440712	44440779	scaffold_8422:9237:9271:-	chr1	44440719	44440739	308663_scaffold_8422:9249-9269_-

## The first 4 fields come from the inputMaf, the last 4 come from the chain.
	
	my %hashReturn = ();
	open(FIO,$overlapFile) || safeDie("The overlap file cannot be opened\n");
	while (my $line = <FIO>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
		
		my($mafChrRef,$mafStartRef,$mafStopRef,$mafCds,$chainCds) = (split /\t/,$line)[0,1,2,3,7];
		
		my $hashEntry = "$mafChrRef\t$mafStartRef\t$mafStopRef";
		my $lineOri = "$mafChrRef\t$mafStartRef\t$mafStopRef\t$mafCds";
		
		## We are splitting something like: scaffold_8422:9237:9271:-
		my ($mafChr,$mafStart,$mafStop,$mafStrand) = (split /:/,$mafCds)[0,1,2,3];
		
		## We are splitting something like: 308663_scaffold_8422:9269-9328_-
		my($chainChr,$chainCdsAndStrand) = (split /:/,$chainCds)[0,1];
		
		# chainChr is "308663_scaffold_8422"
		# chainCds and Strand is "9269-9328_-"
		
		my @tmpChainChr = split(/_/,$chainChr);
		shift(@tmpChainChr);	## Remove the first element which is the block number
		$chainChr = join("_",@tmpChainChr);	# i.e. scaffold_8422
		
		my @tmpChainCdsAndStrand = split(/_/,$chainCdsAndStrand);	# splitting 9269-9328_-
		my $chainStrand = $tmpChainCdsAndStrand[1];					# -
		
		my @tmpCDS = split(/-/,$tmpChainCdsAndStrand[0]);			# splitting 9269-9328
		my $chainStart = $tmpCDS[0];	# 9269
		my $chainStop  = $tmpCDS[1];	# 9328
		
		my $flagOverlap = "";
		if ($mafChr eq $chainChr && $mafStrand eq  $chainStrand)	## This is the first test --> the chromosomes and the strands must match
		{
			$hashReturn{$hashEntry} = "T" if ( ($chainStop >= $mafStart && $chainStop <= $mafStop) || ## chainStop lies between mafStart and mafStop
			($chainStart >= $mafStart && $chainStart <= $mafStop) || ## chainStart lies between mafStart and mafStop
			($chainStart <= $mafStart && $chainStop >= $mafStop) || ## the chain fully encloses the mafStart and mafStop
			($mafStart <= $chainStart && $mafStop >= $chainStop ) ); ## the maf fully encloses the chainStart and chainStop
		}
	}
	close FIO;
	return \%hashReturn;	
}


sub clubELines	## This function clubs the e-lines. An input is an array (a list of line numbers) and the return is a hash
{		## For example if the array looks like 1,3,4,5,7, the return is %hash =(1 => 1, 3 => 3, 7 => 1)
		## The key is the starting point, the value is the length of the merged element (in above example, only 3,4 and 5 are merged)

	my $ref2List = shift;
	my @tmp_array = sort {$a <=> $b} (@$ref2List);	## Arrange the block numbers in ascending order

	my $start_pt = my $posStart = my $mergedLength = "";
	
	my @tmpDumpArray = my %hashCombined = ();
	my $t = 0;

	for (my $i = 0; $i < @tmp_array; $i++)
	{
		if ($i > 0)		## Start checking from the 2nd element onwards
		{
			if ( $tmp_array[$i] != $start_pt+1 )                      ## Now test if the current value and the $start_pt are adjacent
			{
				$mergedLength = scalar(@tmpDumpArray);              

				$posStart  = $tmpDumpArray[0];                  ## Get the key for the hash ==> $posStart which is the first element of @tmpDumpArray
				$hashCombined{$posStart} = $mergedLength;          ## Now put these values in the hash.

				@tmpDumpArray = ();
				$t = 0;
			}
		}

		$tmpDumpArray[$t] = $tmp_array[$i];                     ## Dump every element into the @tmpDumpArray first.
		$t++;

		$start_pt = $tmp_array[$i];                                       ## Initialize the variable $start_pt with the ith index of the input array.
	}

	$mergedLength = scalar(@tmpDumpArray);                                    ## Do this for the last element because it has not been tested for the if statement.
	$posStart  = $tmpDumpArray[0];
	$hashCombined{$posStart} = $mergedLength;
	
	return \%hashCombined;
}

sub safeDie
{
	my $message = shift;
	my $allTempFiles = join(" ",@tempFilesList);
	`rm -rf $allTempFiles`;
	die $message;
}
