#!/sw/bin/perl
## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;

use BerkeleyDB;
use POSIX;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyKentFunctions;
use MyBDB;
use FileHandle;

my $usage = "perl $0 inFile (my GenePred format) ReferenceSpecies speciesList BigBedIndex outDir -verbose [Optional]\n";
die $usage if (scalar(@ARGV) < 5);

my $verbose = 0;
GetOptions ("v|verbose"  => \$verbose); 

my $file 		= $ARGV[0];
my $ref			= $ARGV[1];
my $speciesList 	= $ARGV[2];
my $index  		= $ARGV[3];
my $outDir		= $ARGV[4];

`mkdir -p $outDir`;
my $fileName = `basename $file`; chomp $fileName;

my %sp2File = ();
my @listFiles = qw(threePrimeTable fivePrimeTable threePrimeTableDel fivePrimeTableDel exclude);
## For each species, I create 5 types of files: 3PrimeTable, 5PrimeTable, 3PrimeTableDelete, 5PrimeTableDelete and exclude
foreach my $f(@listFiles) {
	open(FIS,$speciesList) || die "Error opening speciesList file '$speciesList'\n";
	while (my $species = <FIS>) {
		chomp $species;
		my $fh = FileHandle->new();        # create new file handles for each species-outFile pair
		$sp2File{$species}{$f} = $fh;		
		$fh->open(">$outDir/$fileName.$species.$f") || die "ERROR: cannot write to '$outDir/$fileName.$species.$f'\n";
	}
	close FIS;
}

open(FIS,$speciesList) || die "Error opening speciesList file '$speciesList'\n";
my @speciesListArray = <FIS>;
close FIS;

open(FI,$file) || die "Error opening input file '$file'\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$tID,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,3,8];
	my @exonsList = split(",",$exons);
	
	my %hashInfo = my %hashExonNumbers = my %speciesPresentInMaf = ();
	## Loop over all exons and store information about aligning exons in a hash
	my $exonCount = 1;
	for (my $i = 0; $i < scalar(@exonsList); $i++) {
		
		my $exon = $exonsList[$i];
		$hashExonNumbers{$exon} = $i;
		my ($ref2cdsPrint,$ref2chr,$ref2strand,$ref2seqLength) = processMaf($exon,$index,$chr,$ref,$exonCount);  ## This function processes the maf:
		foreach my $species(@speciesListArray) {
			$species =~s/\s+$//;
			$hashInfo{$species}{$i} = "$ref2cdsPrint->{$species}#$ref2chr->{$species}#$ref2strand->{$species}#$ref2seqLength->{$species}"; ## Store everything useful from the maf in the "hashInfo" hash
			$speciesPresentInMaf{$species} = "T";
		}
		$exonCount++;
	}
	
	foreach my $species(@speciesListArray) {
		$species =~s/\s+$//;
		next if (! exists $speciesPresentInMaf{$species});
		## Get all the file handles and use them accordingly
		my $ftp 	= $sp2File{$species}{threePrimeTable};  ## tp is threePrime, fp is fivePrime
		my $ffp 	= $sp2File{$species}{fivePrimeTable};
		my $ftpDel	= $sp2File{$species}{threePrimeTableDel};  ##tpDel is threePrimeDeletion,  fpDel is fivePrimeDeletion
		my $ffpDel	= $sp2File{$species}{fivePrimeTableDel};
		my $fe		= $sp2File{$species}{exclude};
		
		## Make pairs for exons that should be compared: This is important because we could have an alignment to exon1, e-lines for exon2 and an alignment to exon3. 
		## In this case, we have to check that exon1 and exon3 should not overlap
		my @pairsList = ();
		my $j = "";
		for (my $i = 0; $i < scalar(@exonsList)-1; $i++) {
			my $e1 = $i;
			my $e2 = $i + 1;
			my $exon1 = $exonsList[$e1];
			my $exon2 = $exonsList[$e2];
		
			print "Checking the possibility of exonPair --> '$e1' and '$e2'\n" if ($verbose);
		
			my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$i});
			if ($chr1 eq "NA" || $strand1 eq "NA") {
				print "Chr value is '$chr1' and strand value is '$strand1', so excluding\n\n" if ($verbose);
				next;
			}
		
			my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$i+1});
			my $pair = "";
		
			if ($chr1 ne "NA" && $strand1 ne "NA" && $chr1 eq $chr2 && $strand1 eq $strand2){
				$pair = "$e1-$e2"; 
				print "Possible pair $e1 and $e2 will be evaluated. Chr/strands are '$chr1', '$strand1', '$chr2' and '$strand2'\n\n" if ($verbose);
			} else {
				print "Chr value1 is '$chr1' and strand value is '$strand1'; Chr value2 is '$chr2' and strand value is '$strand2' so excluding and checking for downstream pairs\n\n" if ($verbose);
				my $j = $i + 2;
				my $jStart = $j;
			
				print "Scenario2: Evaluating the possibility of '$j' and downstream exons\n\n" if ($verbose);
				if ($j < scalar(@exonsList)){
					for($j = $jStart; $j < scalar(@exonsList); $j++) {
						my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$j});
						if ($chr1 ne "NA" && $strand1 ne "NA" && $chr1 eq $chr2 && $strand1 eq $strand2){
							$pair = "$e1-$j";
							print "This is the new pair '$pair'\n\n" if ($verbose);
							$i = $j;
							last;
						}
					}
				}
			}
			print "'$pair' is a possible exon pair\n" if ($pair ne "" && ($verbose) );
			push(@pairsList,$pair) if ($pair ne "");
		}

		if (scalar(@pairsList) > 0) {
			my %pairsListHash = map{$_ => 1} @pairsList;
			## If the exons are deleted, then the above code would give me a list like: 0-5 6-12 18-35 36-37 37-38 38-39
			## Ensure that the list becomes continuous i.e, 0-5,5-6,6-12,12-18,18-35 and so on
	
			my @pairsListNew = ();
			my %pairChrStrandInfo = ();
			for(my $k = 0; $k < (scalar(@pairsList)-1); $k++) {
				push(@pairsListNew,$pairsList[$k]);  ## i.e. the original pair
				my ($startOriginal,$stopOriginal) = (split /-/,$pairsList[$k])[0,1];
				my($chrOriginal,$strandOriginal) = (split /#/,$hashInfo{$species}{$startOriginal})[1,2];
				$pairChrStrandInfo{$pairsList[$k]} = "$chrOriginal\t$strandOriginal";

				my $start = (split /-/,$pairsList[$k])[1];
				my $stop  = (split /-/,$pairsList[$k+1])[0];				
				my $pairNew = "$start-$stop";

				## Now check if the "new" exon-pairs that I am adding to the list to ensure continuity, come from the same strand and same chromosome
			
				my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$start});
				my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$stop});

				
				if ($chr1 eq $chr2 && $strand1 eq $strand2){
					if (! exists $pairsListHash{$pairNew} && $start != $stop) {
						push(@pairsListNew,$pairNew);
						$pairChrStrandInfo{$pairNew} = "$chr1\t$strand";
					}
				}
			}
		
			push(@pairsListNew,$pairsList[$#pairsList]);  ## Do this for the last element
			###
			my ($pairStart,$pairEnd)    = (split /-/,$pairsList[$#pairsList])[0,1];
			my ($chrLast,$strandLast)   = (split /#/,$hashInfo{$species}{$pairStart})[1,2];
			$pairChrStrandInfo{$pairsList[$#pairsList]} = "$chrLast\t$strandLast";
			 
			########
			
			@pairsList = ();
			@pairsList = @pairsListNew;
			print "These are the continuous exonPairs '@pairsList'\n" if ($verbose);

			my @list = ();
	
			## Now compare adjacent exons from the pairs that I have created above
			for (my $i = 0; $i < scalar(@pairsList); $i++) {
				my($e1,$e2) = (split /-/,$pairsList[$i]);
				my $chrStrandInfo = $pairChrStrandInfo{$pairsList[$i]};
		
				my $exon1 = $exonsList[$e1];
				my $exon2 = $exonsList[$e2];
		
				my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$e1});
				my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$e2});
		
				my $intronLength = $cds2 - $cds1 - $l1;
				my $keyNew = "$exon1<->$exon2<->$intronLength";
				push(@list,$keyNew);
				$pairChrStrandInfo{$keyNew} = $chrStrandInfo;
			}
			
			@list = reverse(@list) if ($strand eq "-");

			if ($verbose) {
				print "####\n\nBelow are the contents and the information associated with exonPairs\n";
				foreach my $element(@list) {
					print "$element\n";
				}
				print "####\n\n";
			}
		 
			### Formatting the output
			foreach my $element(@list) {
				my $e1 = my $e2 = my $length = "";
				my $chrStrandInfo = $pairChrStrandInfo{$element};
			
				if ($strand eq "+") {
					($e1,$e2,$length) = (split /<->/,$element)[0,1,2];
				} else {
					($e2,$e1,$length) = (split /<->/,$element)[0,1,2];
				}
				print "Strand is $strand .. exon1 is '$e1', exon2 is '$e2' while the length is '$length'\n\n" if ($verbose);
			
				if ($length > 29) {
					$length = $length - 30;
					my $lengthIntron = $length/2;
					my $l1 = my $l2 = "";
					
					if ($length%2 == 0) {
						$l1 = $lengthIntron + 23;
						$l2 = $lengthIntron + 7;
					} else {
						$l1 = ceil($lengthIntron) + 23;
						$l2 = floor($lengthIntron) + 7;
					}
					
					print $ftp "$gene\t$tID\t$chr\t$e1\t$species\t$l1\n";
					print $ffp "$gene\t$tID\t$chr\t$e2\t$species\t$l2\n";
					
				} elsif ($length < 29 && $length >= 0){
					my $lengthIntron = $length/2;
					my $l1 = $lengthIntron;
					my $l2 = $lengthIntron;
				
					if ($length%2 != 0){
						$l1 = ceil($lengthIntron);
						$l2 = floor($lengthIntron);
					}
			
					print "Intron deleted!! The lengths are '$l1' and '$l2'\n\n" if ($verbose);
					print $ftp "$gene\t$tID\t$chr\t$e1\t$species\t$l1\n";
					print $ffp "$gene\t$tID\t$chr\t$e2\t$species\t$l2\n";

					print $ftpDel "$gene\t$tID\t$chr\t$e1\t$species\t$l1\t$chrStrandInfo\n";
					print $ffpDel "$gene\t$tID\t$chr\t$e2\t$species\t$l2\t$chrStrandInfo\n";
				} elsif ($length < 0) {
					print $fe "$gene\t$tID\t$chr\t$e1\t$species\t$length\n";
				        print $fe "$gene\t$tID\t$chr\t$e2\t$species\t$length\n";					
					print "Negative intron length !! The length is '$length'\n\n" if ($verbose);
				}
			}
		}
	}
}
close FI;

## close all fileHandles
foreach my $f(@listFiles) {
	open(FIS,$speciesList) || die "Error opening speciesList file '$speciesList'\n";
	while (my $species = <FIS>) {
		chomp $species;
		my $fh = $sp2File{$species}{$f};
		$fh->close;
	}
	close FIS;
}

sub processMaf {
	my($exon,$index,$chr,$ref,$exonCount) = @_;
	
	## Get the maf for the locus of interest using mafExtract
	my $call   = "mafExtract $index -region=$chr:$exon stdout";
	my $tmpMaf = `$call`;
	die "ERROR: '$call' failed while running mafExtract \n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	my @mafExtractResult = split(/\n/,$tmpMaf);

	print "'$call' was executed, the chr is '$chr' and the exon coordinates are '$exon', the exonNumber is '$exonCount'.\n" if ($verbose);
	print "The following maf is extracted\n$tmpMaf\n" if ($verbose);

	my %speciesInList = ();  ## Create a hash that contains only those species that are present in the input list. Useful if you want to process, say 73 mammals from the 145-way alignment
	my %chrHash = my %strHash = my %cdsPrint = my %seqLength = my %ct = ();
	my %chrFinalHash = my %strFinalHash = ();
	
	open(FIS,$speciesList) || die "Error opening speciesList file\n";
	while (my $species = <FIS>) {
		$species =~s/\s+$//;
		$speciesInList{$species} = "T";
		$seqLength{$species} = 0;
		$ct{$species} = 0; ## A species specific counter
		
		## Initialize stuff to NA
		$chrFinalHash{$species} = "NA";
		$strFinalHash{$species} = "NA";
		$cdsPrint{$species} = "NA";
	}
	close FIS;

	my $i = 0;
	foreach my $line(@mafExtractResult) {
		$line =~s/\s+$//;
		if ($line =~/^s/ && $line !~/$ref/){
			my($sp,$chr,$cdsStart,$lengthAlnSeq,$strand,$srcSize,$seq) = getSLineData($line);
			next if (! exists $speciesInList{$sp});
			
			$chrHash{$sp}{$chr} = "T";
			$strHash{$sp}{$strand} = "T";
			
			$cdsPrint{$sp} = $cdsStart if ($ct{$sp} == 0); ## Get the beginning of the CDS. i.e. the first block from the maf
			$seqLength{$sp} = $seqLength{$sp} + $lengthAlnSeq; ## Get the length of the aligning sequence through all blocks in the maf
			$ct{$sp}++;
		}
	}
	
	## Now populate chrFinalHash and strFinalHash.
	## chrFinalHash: The key is the species, the value is the chromosome. If more than one, then the value is "NA"
	foreach my $species(keys (%chrHash)) {
		my @chrList = keys %{$chrHash{$species}};
		my %tmpHash = map{$_ => 1}@chrList;
		
		if (scalar(keys(%tmpHash)) == 1){
			$chrFinalHash{$species} = $chrList[0];
			print "Species '$species' has only one chromosome in the aligning block '$chrList[0]'\n" if ($verbose);
		} else {
			$chrFinalHash{$species} = "NA";
			print "Species '$species' has more than one chromosome in the aligning block '@chrList', so will be excluded\n" if ($verbose);
		}
	}
	
	## strFinalHash: The key is the species, the value is the strand. If more than one, then the value is "NA"
	foreach my $species(keys (%strHash)) {
		my @strList = keys %{$strHash{$species}};
		my %tmpHash = map{$_ => 1}@strList;
		if (scalar(keys(%tmpHash)) == 1){
			$strFinalHash{$species} = $strList[0];
			print "Species '$species' has only one strand in the aligning block '$strList[0]'\n" if ($verbose);
		}else {
			$strFinalHash{$species} = "NA";
			print "Species '$species' has more than one strand in the aligning block '@strList', so will be excluded\n" if ($verbose);
		}
	}
	
	if ($verbose) {
		foreach my $species (keys(%speciesInList)) {
			print "The coordinates are $cdsPrint{$species}, the exonNumber is '$exonCount'\n";
		}
		print "#############\n\n";
	}
	
	return (\%cdsPrint,\%chrFinalHash,\%strFinalHash,\%seqLength);	
}
