#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my $usage 	= "perl allTranscriptsFile 5PrimeTable 3PrimeTable outFile exonsUnionFile species excludeExons -verbose [Optional]\n";
my $purpose	= "This script create exonGroups file i.e. exons which should be grouped together as a single unit and then aligned with CESAR. This script reads in 4 files - the all transcripts file, the file containing information about intron which is 5' of an exon, the file containing information about the intron which is 3' of the exon and the outFile\n";

die $usage.$purpose if (scalar(@ARGV) < 7);
my $verbose = 0;
GetOptions ("v|verbose"  => \$verbose); 

my $allTranscripts  = $ARGV[0];
my $fivePrimeDel    = $ARGV[1];
my $threePrimeDel   = $ARGV[2];
my $outFile         = $ARGV[3];
my $exonUnionFile   = $ARGV[4];
my $species         = $ARGV[5];
my $excludeExons    = $ARGV[6];

my ($ref2Hash3Prime,$ref2Hash3PrimeLocus) = fileHash($threePrimeDel);  ## convert the threePrimeDel file to hash
my ($ref2Hash5Prime,$ref2Hash5PrimeLocus) = fileHash($fivePrimeDel); 

my %hash3Prime = %$ref2Hash3Prime;
my %hash5Prime = %$ref2Hash5Prime;
my %hash5PrimeLocus = %$ref2Hash5PrimeLocus;  ## Get a queryLocusHash that tells us the locus(scaffold plus the strand) for every exon in the query genome.
my %hash3PrimeLocus = %$ref2Hash3PrimeLocus;

open(FI,$allTranscripts) || die "Error opening allTranscripts file '$allTranscripts'\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$transcript,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,3,8];
	my @exonsList = split(",",$exons);
	@exonsList = reverse(@exonsList) if ($strand eq "-");
	
	my %queryLocusAll = ();  ## Store all queryLoci in this hash
	## Let's store everything in queryLocus now:
	my $ct = 0;
	foreach my $ex(@exonsList) {
		my $key = "$transcript#$chr#$ex";
		$queryLocusAll{$ct} = $hash5PrimeLocus{$key} if (exists $hash5PrimeLocus{$key});
		$queryLocusAll{$ct} = $hash3PrimeLocus{$key} if (exists $hash3PrimeLocus{$key});
		$ct++;
	}
	
	my @exonPairs = ();
	my @deletedExons = ();  ## This list contains exons that are deleted
	for(my $j=0; $j < (scalar(@exonsList)-1); $j++) {
		my $key = "$transcript#$chr#$exonsList[$j]";
		print "Checking for key '$key', exonNumber '$j' if it is associated with an intron deletion\n\n" if ($verbose);
		my $queryLocus = $hash3PrimeLocus{$key};
		
		if (exists $hash3Prime{$key}) {  ## If a given exon exists in 3PrimeDel hash i.e. the intron downstream of an exon is deleted
			print "Exon '$key', exonNumber '$j' is associated with a 3Prime intron deletion, now checking it the next exon has a 5Prime intronDeletion\n\n" if ($verbose);
			my $k = $j+1;
			
			my $key3PrimeExon = "$transcript#$chr#$exonsList[$k]"; ## check the key for its downstream exon, it should exist in the 5PrimeDel hash
			if (exists $hash5Prime{$key3PrimeExon}) {	
				push(@exonPairs,$j,$k);
				print "Yes!! The downsteam exon '$key3PrimeExon' is associated with a 5Prime intron deletion\n\n" if ($verbose);
			} else { ## the intron deletion is followed by exon(s) deletion as well.
			
				push(@deletedExons,$k); ## otherwise, there is an exon deletion and possibly the next exon exists in the 3PrimeDel hash
				## What is the queryLocus for a deleted exon? The queryLocus is the same as the upstream exon has
				$queryLocusAll{$k} = $queryLocus;
				
				print "Pushing exon '$k' in deletedExons list in firstPass\n" if ($verbose);
				for (my $l = $k+1; $l < (scalar(@exonsList) -1); $l++) {
					my $key3PrimeExon = "$transcript#$chr#$exonsList[$l]";
					if (exists $hash5Prime{$key3PrimeExon}) {  ## i.e only an intron between two adjacent exons is deleted
						
						print "Exon-pairs list has the following '@exonPairs'\n" if ($verbose);
						
						print "Here I find the exon whose 5Prime intron is deleted. The number is '$l'\n" if ($verbose);
						push(@exonPairs,$j);  ## Push the exon that we started with to the exonPairs list
						push(@exonPairs,@deletedExons); ## Then push all the intermediate exons (that are deleted)
						push(@exonPairs,$l); ## And lastly push the exon whose 5Prime intron is deleted. i.e. exon which closes the pair
						
						$j = $l;
						print "The deletedExonsList has the following exons '@deletedExons'\n" if ($verbose);
						last;
					} else {
						push(@deletedExons,$l);   ## keep pushing these to deletedExons list until you find the exon whose 5' intron is deleted.
						print "Pushing exon '$l' in deletedExons list in subsequent passes\n" if ($verbose);
						$queryLocusAll{$l} = $queryLocus;
					}
				}
			}
		}
	}
	my %deletedExonsHash = ();
	
	
	print "###\nThese exons make it to the exonPairs group '@exonPairs'\n###\n" if ($verbose);
	
	## An example: the human gene GGT1; In rat, the intron downstream of exon2 is deleted, exon3 is deleted and intron downstream of exon4 is deleted.

	if (scalar(@exonPairs) > 0) {
		## Filter exon pairs and keep only those that are on the same strand and scaffold
		
		my @allExonsLocusList = ();                              ## This array partitions exonPairs and groups them bases on the queryLocus. For example if exons1,2 and 3 come from scaffold:1 on plus strand
		my @queryLocusValues = values(%queryLocusAll);           ## and exons 4 and 5 come from scaffold:2 on plus strand, then the exons will be grouped
		my %queryLocusUniqHash = map{$_ => 1} @queryLocusValues; ## Element1 of this array will contain reference to an array which has exons 1,2 and 3 while element2 will contain a reference with exons 4 and 5.
		my @queryLocusValuesUniq = keys(%queryLocusUniqHash);
		
		foreach my $queryLocus(@queryLocusValuesUniq)			 {								
			print "This is the queryLocus --> '$queryLocus'\n" if ($verbose);
			my @allExonsForALocus = ();
			foreach my $keys(keys(%queryLocusAll)) {
				push(@allExonsForALocus,$keys) if ($queryLocusAll{$keys} eq $queryLocus);
			}
			
			push(@allExonsLocusList,\@allExonsForALocus);
		}
		
		foreach my $element(@allExonsLocusList) {
			my $ref2MergedList = mergeAdjacentElements($element);  ## For GGT1, the mergedList will contain 2-3-4
			print "This is the mergedList '@$ref2MergedList'\n" if ($verbose);
			printFunction($ref2MergedList,$outFile,\@exonsList,$gene,$transcript);
		}
	}
}
close FI;

## Now format the output file so that it is ready for  RealignmentPipeline_Pairs.pl
formatFile($outFile,$exonUnionFile,$species,$excludeExons);

### Sub-routines
#############################
sub fileHash {  ## Converts an inFile to hash
	my $file = shift;
	my %hash = my %hashLocus = ();
	open(FI,$file) || die "Error opening inFile '$file' in fileHash function\n";
	while (my $line = <FI>) {
		$line =~s/\s+$//;
		my($transcript,$chr,$exonCds,$distance,$chrQuery,$strQuery) = (split /\t/,$line)[1,2,3,5,6,7];
		my $queryInfo = "$chrQuery\t$strQuery";
		my $key = "$transcript#$chr#$exonCds";
		
		$hash{$key} = "$distance";
		$hashLocus{$key} = $queryInfo;
		
	}
	close FI;
	return (\%hash,\%hashLocus);
}

sub mergeAdjacentElements {
	my $Ref2List = shift;                           
	my @tmp_array = @$Ref2List;
	
	## In case, there are more than 2 exons that need to be merged, the input list will look like the following:
	## @tmp_array = qw(1 2 3 3 4 4 5), it should become @tmp_array = qw(1 2 3 4 5);
	my %tmpHash = map{$_ => 1} @tmp_array;
	@tmp_array = ();
	@tmp_array = sort {$a <=> $b} (keys(%tmpHash));
	
	my @tmp_dump_array = ();
	my $t = 0;

	my $start_pt = my $pos_start = "";
	my @finalList = ();
	for (my $i = 0; $i < @tmp_array; $i++) {
		if ($i > 0) {
			if ($tmp_array[$i] != $start_pt+1) {                      ## Now test if the current value and the $start_pt are adjacent
				my $string = join("-",@tmp_dump_array);
				push(@finalList,$string);
				
				@tmp_dump_array = ();
				$t = 0;
			}
		}
		$tmp_dump_array[$t] = $tmp_array[$i];                     ## Dump every element into the @tmp_dump_array first.
		$t++;

		$start_pt = $tmp_array[$i];                                       ## Initialize the variable $start_pt with the ith index of the input array.
	}

	my $string = join("-",@tmp_dump_array);  ## For the last element of the array
	push(@finalList,$string);
	
	return \@finalList;
}

sub formatFile {
	my($exonGroupsFile,$exonsUnionFile,$species,$excludeExons) = @_;
	
	my %hashAllExons = ();
	open(FI,$exonsUnionFile) || die "Error opening exonUnion file '$exonsUnionFile'\n";
	while (my $line = <FI>) {
		$line =~s/\s+$//;
		my($gene,$accList,$chr,$strand,$phase,$tag,$cds) = split(/\t/,$line);
		my @accListAll = split(",",$accList);

		foreach my $acc(@accListAll) {
			$hashAllExons{"$acc\_$cds"} = "$chr\t$strand\t$phase\t$tag";
		}
	}
	close FI;

	my $tmpFile = `mktemp`; chomp $tmpFile;
	
	## Read in all exclude exons
	my %excludeExons = ();
	open(FIE,"$excludeExons") || die "Error opening exclude exons file\n";
	while (my $line = <FIE>) {
		chomp $line;
		my($chr,$exon) = (split /\t/,$line)[2,3];
		$excludeExons{$chr}{$exon} = "T";
	}
	close FIE;

	open(FI,$exonGroupsFile) || die "Error opening temporary exonGroup file '$exonGroupsFile'\n";
	open(FO,">$tmpFile") || die "Error cannot write to temporary output file\n";
	open(FOE,">>$excludeExons") || die "Error appending to excludeExons file\n";
	
	while (my $line = <FI>) {
		$line =~s/\s+$//;
		my($gene,$acc,$exons) = (split /\t/,$line);
		my @exonsList = split(",",$exons);
		
		my @phaseList = my @tagList = ();
		my $chrPrint = my $strandPrint = "";
		foreach my $ex(@exonsList) {
			my $key = "$acc\_$ex";
			my($chr,$strand,$phase,$tag) = (split /\t/,$hashAllExons{$key});
			$chrPrint = $chr;
			$strandPrint = $strand;
		
			push(@phaseList,$phase);
			push(@tagList,$tag);
		}
		my $phaseListPrint = join(",",@phaseList);
		my $tagListPrint   = join(",",@tagList);

		my $sumAllExons = 0;
		my $allOK_exclude = my $allOK = "";
		foreach my $ex(@exonsList) {
			if (exists $excludeExons{$chrPrint}{$ex}) {
				$allOK_exclude = "T";
				last;
			}
		}
		
		for(my $j = 0; $j < scalar(@exonsList); $j++) {
			my($start,$stop) = (split /-/,$exonsList[$j]);
            $sumAllExons = $sumAllExons + ($stop - $start);			

			if ($tagList[$j] eq "F" || $tagList[$j] eq "L") {
				my $l = $stop - $start;
				if ($l < 6) {
					$allOK = "F";
					last;
				}
			}
		}
		
		$allOK = "F" if ($sumAllExons > 2000);

		if ($allOK eq "F" || $allOK_exclude eq "F") {
			foreach my $ex(@exonsList) {
				print FOE "$gene\t$acc\t$chrPrint\t$strandPrint\t$ex\n";
			}
		} else {
			print FO "$gene\t$acc\t$chrPrint\t$strandPrint\t$exons\t$phaseListPrint\t$tagListPrint\t$species\n";
		}
	}
	close FI;
	close FO;
	close FOE;

	`mv $tmpFile $exonGroupsFile`;
}

sub printFunction {  ## prints to the output file
	my ($ref2MergedList,$outFile,$ref2exonsList,$gene,$transcript) = @_;
	my @exonPairsList = @$ref2MergedList;
	my @exonsList = @$ref2exonsList;
		
	my @printList = ();
	foreach my $pair(@exonPairsList) {
		my @exonCdsList = ();
		my @pairsL = split("-",$pair);
		
		foreach my $p(@pairsL) {
			push(@exonCdsList,$exonsList[$p])
		}
		
		my $string = join(",",@exonCdsList);
		push(@printList,$string);
	}
	
	open(FO,">>$outFile") || die "Error appending to outFile '$outFile' in printFunction\n";
	foreach my $p(@printList) {
		print FO "$gene\t$transcript\t$p\n";
	}
	close FO;
}
