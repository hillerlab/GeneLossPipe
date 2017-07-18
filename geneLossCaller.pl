#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany
## geneLossCaller.pl reads the mutations file produced by GeneLossPipeline.pl and classifies a gene as lost or not lost
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;

$| = 1;		# == fflush(stdout)

my $verbose = my $compFS = my $longIndelFlag = 0;
GetOptions ("verbose" => \$verbose, "compFS" => \$compFS, "longIndel" => \$longIndelFlag);

my $usage = "USAGE: $0 genesList(one gene per line, the first field is the gene name, the second is a string of principal isoforms which are separated by comma) GeneLossPipeBDB annotationFile species ReadingFrameThreshold exonGroupsDir u12IntronList outFile 
Optional parameters: -compFS (Count frameshifts as inactivating mutations even if they are compensated) -verbose -longIndel [Treat every indel longer than 50 bp as an inactivating mutation]\n";
my $purpose = "PURPOSE: $0 reports gene losses given a list of genes, the geneLossPipe outDirectory, the annotation file, the species and the RF threshold (for example 0.6)\n";

die $usage.$purpose if (scalar(@ARGV) != 8);

my $list          = $ARGV[0];
my $geneLossBDB   = $ARGV[1];
my $annoFile      = $ARGV[2];
my $species       = $ARGV[3];
my $rfThreshold   = $ARGV[4];
my $exonGroupsDir = $ARGV[5];
my $u12IntronList = $ARGV[6];
my $outFile       = $ARGV[7];

## Check parameters
die "Input geneList does not exist\n" if (! -e $list);
die "GeneLoss pipeline BDB file does not exist\n" if (! -e $geneLossBDB);
die "GeneLoss pipe BDB ($geneLossBDB) is not a bdb file\n" if ($geneLossBDB !~/bdb$/i);
die "geneAnnotation file does not exist\n" if (! -e $annoFile);

### Length business, get the information and store it in relevant hashes. With this information, I can get the relative position for every mutation in a gene
my %cumLengthHash = my %geneLengthHash = my %exonLengthHash = (); 
my %exonNumbersHash = my %exonCtHash = ();
my %exonSizePropHash = ();

open(FI,$annoFile) || die "Error opening annotationFile '$annoFile'\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($acc,$strand,$exons) = (split /\t/,$line)[1,3,8];
	my @exonsList = split(",",$exons);
	@exonsList = reverse(@exonsList) if ($strand eq "-");
	
	my $ref2CumLengthHash = getCumLength(\@exonsList);  ## cumulative length hash
	my ($geneLength,$ref2ExonLengthHash,$ref2ExonPropLengthHash) = getGeneLength(\@exonsList);

	$cumLengthHash{$acc}  	= $ref2CumLengthHash;
	$geneLengthHash{$acc} 	= $geneLength;
	$exonLengthHash{$acc} 	= $ref2ExonLengthHash;
	$exonSizePropHash{$acc}	= $ref2ExonPropLengthHash;
	
	my $ct = 1;
	foreach my $ex(@exonsList) {
		$exonNumbersHash{$acc}{$ex} = $ct;
		$ct++;
	}
	$exonCtHash{$acc} = $ct - 1;
}
close FI;

## Now read the exonGroups file (where introns are deleted and store this information in a hash)
my %intronDeletionHash = ();

my $exonGroupFile = "$exonGroupsDir/$species.groups";
open(FI,$exonGroupFile) || die "Error opening exonGroups file '$exonGroupFile'\n";
while (my $line = <FI>) {
	my($acc,$exons) = (split /\t/,$line)[1,4];
	my @exonsList = split(/,/,$exons);
	
	my @exonNumbers = ();
	foreach my $ex(@exonsList) {
		my $number = $exonNumbersHash{$acc}{$ex};
		push(@exonNumbers,$number);
	}
	
	@exonNumbers = sort{$a <=> $b} @exonNumbers;
	for(my $j = 0; $j < scalar(@exonNumbers); $j++) {
		my $j1 = my $j2 = "";
		my $ex = $exonNumbers[$j];
		
		$j1 = "$ex\_D" if ($j != scalar(@exonNumbers)-1);
		$j2 = "$ex\_A" if ($j != 0);
		
		$intronDeletionHash{$species}{$acc}{$j1} = "T" if ($j1 ne "");
		$intronDeletionHash{$species}{$acc}{$j2} = "T" if ($j2 ne "");
	}
}
close FI;

### Store the information about U12 introns in a hash:
my %u12IntronHash = ();
open(FI,$u12IntronList) || die "Error opening U12 intron list\n";
while (my $line = <FI>) {
	
	$line =~s/\s+$//;
	my $accPlusJ = (split /\t/,$line)[6];
	my @tmp = split(/_/,$accPlusJ);
	
	my $acc = shift(@tmp);
	my $junction = join("_",@tmp);
	
	$u12IntronHash{$acc}{$junction} = "T";
}
close FI;

## For each species, I compute the number of geneLosses according to the following criteria:
## If the gene has at least 2 mutations, spread over atleast 2 exons and the %age of unaffected reading frame is less than 60% (or any other threshold)
## For single exon gene, I just cannot have 2 mutations spread over 2 exons but I still ask if there are 2 mutations

open(FIGL,$list) || die "Error opening input geneList '$list'\n";
open(FOGL,">$outFile") || die "Error writing to the output file '$outFile'\n";

while (my $lineIn = <FIGL>) {
	$lineIn =~s/\s+$//;
	my($geneName,$isoforms) = (split /\t/,$lineIn);
	
	my @isoformsList = split(",",$isoforms);
	my $lost = 0;
	my @mutPositionsIsoform = ();
	print "###\n\nChecking gene loss for gene '$geneName'. The gene has the following principalIsoforms '$isoforms', the species of interest is '$species' ###\n\n" if ($verbose);
	
	foreach my $transcript(@isoformsList) {
		## get mutations for every transcript
		my $mutationsGeneLoss = readBDB($geneLossBDB,$transcript);
		my @mutationsListGeneLoss = split(/\n/,$mutationsGeneLoss);
		
		next if (scalar(@mutationsListGeneLoss) < 1);
		
		print "Checking for transcript '$transcript'\n" if ($verbose);
		
		my $noe = $exonCtHash{$transcript};
		my $nMut = my $nExons = "";  ### I compute my nMut and nExons based on the number of exons in the transcript
		
		if ($noe == 1) {
			$nExons = 1;
			$nMut = 2;
		} elsif ($noe <= 10) {
			$nExons = $nMut = 2;
		} else {
			$nExons = 0.2 * $noe;
			$nExons = int($nExons + 0.5);
			$nMut = $nExons;
		}
		
		my $ref2ExonPropHash = $exonSizePropHash{$transcript};
		my $geneLength	= $geneLengthHash{$transcript};
		my $ref2ExonLengthHash = $exonLengthHash{$transcript};
		my %exonLengthHash = %$ref2ExonLengthHash;
		my $numberOfExons = scalar(keys(%exonLengthHash));
		
		my $localHash = $cumLengthHash{$transcript};
		my %cumLengthHash = %$localHash;
		my @mutationsPosList = my @mutatedExonsList = ();
		my @allInactivatingMutations = ();
		
		## Get a list of compensatory events 
		my %compEventsHash = ();
		if ($compFS == 0) {
			my $ref2CompEvents = getCompensatoryEvents_NoExonDel(\@mutationsListGeneLoss,$species);			
			 %compEventsHash = %$ref2CompEvents;
		}
		## Otherwise we do not store any information in the compensatory events hash
		
		my $exonDeletionCount = 0;
		my %longIndelCountHash = ();
		my $x = 1;
		while ($x <= $noe) {
			$longIndelCountHash{$x} = 0;
			$x++;
		}
		
		foreach my $line(@mutationsListGeneLoss) {
		
			$line =~s/\s+$//;
			my $speciesInLine = (split /\t/,$line)[0];
			if ($speciesInLine eq $species){
				if ($line !~/Compensation/ && $line !~/Entire exon deletion/ && $line !~/Entire exon sequence missing/){  ## Everything else is a mutation. Exon deletions are dealt separately
				
					my($species,$mutID,$event,$exonNumber,$pos) = (split /\t/,$line)[0,2,3,4,5];
					
					if ($event eq "Insertion" || $event eq "Deletion"){
						my $length = (split /\t/,$line)[6];  ## the length of the indel
						my $posGene    = $pos + $cumLengthHash{$exonNumber};
						my $posGeneEnd = 0;
						
						if ($event eq "Deletion"){  ## find positions of deletions by getting a mean of (start of mutation + (start of mutation + length of deletion))
							$posGeneEnd = $posGene + $length;
							$posGene = sprintf("%0.3f",(($posGene+$posGeneEnd)/2));
						}
								
						my $relPosition = sprintf("%0.3f",($posGene/$geneLength));
						my $longIndel = "T" if ($length > 50 && ($longIndelFlag));
						
						if ( ($length%3 != 0 && ! exists $compEventsHash{$mutID}) || $line =~/StopCodon/ || ($longIndel) ){  ## Proceed with this mutation if 1) it is a frameshift that is not compensated
																     ## 2) if it is a frame-preserving insertion but "brings in a stop codon"	
							push(@allInactivatingMutations,$relPosition);  ## Store all inactivating mutations in this list
							
							if ($relPosition > 0.2 && $relPosition < 0.8) { 
								push(@mutationsPosList,$relPosition); ## Label a mutation as inactivating if it is in the middle 60% of the gene
								push(@mutatedExonsList,$exonNumber); ## If the mutation is inactivating, then the exon from where it comes is also a mutated exon	
								print "This is an inactivating mutation: type '$event', mutationID '$mutID', relPosition '$relPosition', exonNumber '$exonNumber', the length is '$length'\n" if ($verbose);
								$longIndelCountHash{$exonNumber}++ if (($longIndel) && $length%3 == 0 && ($line !~/StopCodon/));
							} else {
								print "The event '$mutID' is a frameshift but is excluded because the relPosition is '$relPosition'\n" if ($verbose);
							}
						} else {
							print "The event '$mutID' is excluded because: 1) either it is not a frameshift 2) it is a frameshift but is compensated by another frameshift. The length is '$length'\n" if ($verbose);
						}
					} elsif ($event eq "Splice site mutation"){  ## Print splice site mutations that are not caused by intron deletions and do not come from U12 introns
						my $junction = (split /\t/,$line)[5];
						my $j = substr($junction,0,1);
						my $check = "$exonNumber\_$j";
					
						if (! exists $intronDeletionHash{$species}{$transcript}{$check} && ! exists $u12IntronHash{$transcript}{$check}){
							my $posGene = "";
							if ($j eq "D" ){
								$posGene = $cumLengthHash{$exonNumber+1};
							}elsif ($j eq "A"){
								$posGene = $cumLengthHash{$exonNumber};
							}else{
								die "Unknown junction '$j' for transcript --> '$transcript', species --> '$species'\n";
							}
							
							my $relPosition = sprintf("%0.3f",($posGene/$geneLength));
							push(@allInactivatingMutations,$relPosition);
							
							if ($relPosition > 0.2 && $relPosition < 0.8) { 
								push(@mutationsPosList,$relPosition); ## Label a mutation as inactivating if it is in the middle 60% of the gene
								push(@mutatedExonsList,$exonNumber); ## If the mutation is inactivating, then the exon from where it comes is also a mutated exon	
								print "This is an inactivating mutation: type '$event', mutationID '$mutID', relPosition '$relPosition', exonNumber '$exonNumber'\n" if ($verbose);
							} else {
								print "The splice site mutation '$mutID' is excluded because its relPosition is '$relPosition'\n" if ($verbose);
							}
							
						} else {
							print "The spliceSite mutation '$mutID' is either a U12 intron or an intron deletion\n" if ($verbose);
						}
					}elsif ($event  eq "In-frame stop codon"){  ## Print in-frame stop codons that are NOT in thre first|last 20% of the sequence 
						my $posGene 	= $pos + $cumLengthHash{$exonNumber};
						my $relPosition = sprintf("%0.3f",($posGene/$geneLength));
						
						push(@allInactivatingMutations,$relPosition);
						
						if ($relPosition > 0.2 && $relPosition < 0.8) { 
							push(@mutationsPosList,$relPosition); ## Label a mutation as inactivating if it is in the middle 60% of the gene
							push(@mutatedExonsList,$exonNumber); ## If the mutation is inactivating, then the exon from where it comes is also a mutated exon
							print "This is an inactivating mutation: type '$event', mutationID '$mutID', relPosition '$relPosition', exonNumber '$exonNumber'\n" if ($verbose);
						} else {
							print "The IFSC '$mutID' was excluded because the relPosition is '$relPosition'\n" if ($verbose);
						}
					}
				} elsif ($line =~/Entire exon deletion/){
					$exonDeletionCount++;
					my($mutID,$exonNumber,$exonLength) = (split /\t/,$line)[2,4,5];
						
					my $posGene = $cumLengthHash{$exonNumber};
					my $posGeneEnd = $posGene + $exonLength;
					$posGene = sprintf("%0.3f",(($posGene+$posGeneEnd)/2));
						
					my $relPosition = sprintf("%0.3f",($posGene/$geneLength));
					push(@allInactivatingMutations,$relPosition);
						
					if ($relPosition > 0.2 && $relPosition < 0.8) { 
						push(@mutationsPosList,$relPosition); ## Label a mutation as inactivating if it is in the middle 60% of the gene
						push(@mutatedExonsList,$exonNumber); ## If the mutation is inactivating, then the exon from where it comes is also a mutated exon
						print "This is an inactivating mutation: type 'Entire exon deletion', mutationID '$mutID', relPosition '$relPosition', exonNumber '$exonNumber'\n" if ($verbose);
					} else {
						print "The exonDeletion '$mutID' was excluded because the relPosition is '$relPosition'\n" if ($verbose);
					}
				}
			}

			## check whether the gene is lost using the criteria desrcibed above:
			@mutationsPosList = sort {$a <=> $b} @mutationsPosList;  ## Sort in ascending order
			my $inactMutCount = scalar(@mutationsPosList);
			print "$species\t$geneName\t$transcript\t$inactMutCount\tInactivatingMutationCountString\n" if ($verbose);
			
			## Get a unique list of exons that are mutated, just convert the list to a hash and get back the keys
			my %tmpHash = map{$_ => 1} @mutatedExonsList;
			@mutatedExonsList = keys(%tmpHash);
			
			if ($verbose) {
				print "The list of inactivating mutations is: @mutationsPosList\n";
				print "The list of mutated exons is: @mutatedExonsList\n\n";
			}
			
			my $geneLost = checkForGeneLoss(\@mutationsPosList,\@mutatedExonsList,$noe,$rfThreshold,$nMut,$nExons,$ref2ExonPropHash,\@allInactivatingMutations,\%longIndelCountHash);  ## This function returns a "T" 
			
			if ($verbose) {
				print "Transcript '$transcript' called as lost\n\n\n" if ($geneLost eq "T");
				print "Transcript '$transcript' NOT called as lost\n\n\n" if ($geneLost eq "");
			}
			
			if ($noe == 1 && $exonDeletionCount == 1) { ## This condition is not covered, so test it here.
				print "Transcript '$transcript' called as lost because it has only 1 exon and the exon is deleted\n\n\n" if ($verbose);
				$geneLost = "T";
			}
			
			$lost++ if ($geneLost eq "T");
		}
	}
	print FOGL "$geneName\n" if ($lost == scalar(@isoformsList));  ## i.e. We call a gene as lost if all the isoforms for a gene are lost.
	print "'$lost' principal isoforms were called as lost\n" if ($verbose);
}
close FIGL;


###################
## Sub-routines  ##
###################
sub getCumLength { ## get cumulative length at the end of each exon. Returns a hash where the key is the exon number, the value is the cumulative length

	my $ref2ExonsList = shift;
	my @exonsList = @$ref2ExonsList;
	
	my %cumLengthHash = ();
	my $lengthExon = 0;
	my $i = 1;
	
	my $j = 0;
	while ($j <= scalar(@exonsList)) {
		$cumLengthHash{$j} = 0;
		$j++;
	}
	
	foreach my $ex(@exonsList) {
		$cumLengthHash{$i} = $cumLengthHash{$i - 1} + $lengthExon;
		
		my ($exonStart,$exonEnd) = (split /-/,$ex)[0,1];
		$lengthExon = abs($exonStart - $exonEnd);
		$i++;
	}
	return \%cumLengthHash;
}

sub getGeneLength { ## Also returns different length statistics, such as the length of the gene, a hash which contains exon length, another hash which contains the proportion of each exon relative to the size of the whole transcript

	my $ref2ExonsList = shift;
	my @exonsList = @$ref2ExonsList;
	
	my $length = 0;
	my %exonLengthHash = ();
	my $e = 1;
	
	foreach my $ex(@exonsList) {
		my ($exonStart,$exonEnd) = (split /-/,$ex)[0,1];
		my $lengthExon = abs($exonStart - $exonEnd);
		$exonLengthHash{$e} = $lengthExon;
		$length = $length + $lengthExon;
		$e++;
	}
	
	$e = 1; my %exonLengthPropHash = ();
	foreach my $ex(@exonsList) {
		my ($exonStart,$exonEnd) = (split /-/,$ex)[0,1];
		my $lengthExon = abs($exonStart - $exonEnd);
		$exonLengthPropHash{$e} = sprintf("%0.3f",($lengthExon/$length));
		$e++;
	}
		
	return ($length,\%exonLengthHash,\%exonLengthPropHash);
}

sub checkForGeneLoss { ## This function decides whether we call a mutation profile for a transcript-species pair as a gene loss or not.

	my ($ref2MutList,$ref2MutExonsList,$noe,$maxRFThreshold,$mutCountThreshold,$mutExonsCountThreshold,$ref2ExonPropHash,$ref2AllMutations,$ref2longIndelCountHash) = @_;
	my @mutList = @$ref2MutList;
	my $nMut = scalar(@mutList);
	my %longIndelCountHash = %$ref2longIndelCountHash;
	
	my @mutatedExons = @$ref2MutExonsList;
	my %exonPropHash = %$ref2ExonPropHash;
	my @allInactivatingMutations = @$ref2AllMutations;
	@allInactivatingMutations = sort {$a <=> $b}@allInactivatingMutations;
	
	my $geneLossS1 = my $geneLossS2 = "";	
	my $longExonFlag = "";
	
	## set geneLossS2 to True if the transcript meets the number of mutated exons criteria.
	if ($noe > 1) { ## For genes with more than one exon, set this to flag to true only if the number of exons affected is more than equal to our mutatedExonsCount threshold

		if (scalar(@mutatedExons) == 1) {  ## If only one exon is mutated
			if ($exonPropHash{$mutatedExons[0]} >= 0.4) {  ## Set this to true if only one exon is affected but this exon is more than 40% of the gene sequence
				$geneLossS2 = "T";
				$longExonFlag = "T";
				print "Even though only exon is mutated, this exon is more than 40% of the total size of the protein. So these mutations are considered to be inactivating\n" if ($verbose);
			}
		} else {
			$geneLossS2 = "T" if (scalar(@mutatedExons) >= $mutExonsCountThreshold);
		}
	} else { ## For single exon genes, set this to true anyway
		$geneLossS2 = "T";	
	}
	
	if ($verbose) {
		print "Long exon flag --> '$longExonFlag'\n";
		print "The number of mut is '$nMut' (includes long frame-preserving indels)\n";
	}
	$nMut = $nMut - $longIndelCountHash{$mutatedExons[0]} if ($longExonFlag eq "T");
	print "After --> The number of mut is '$nMut'\n" if ($verbose);
	
	## set geneLossS1 to True if the number of inactivating mutations threshold and %intact reading frame threshold for calling geneLoss is met
	if ($nMut >= $mutCountThreshold) {  ## The gene has to have at least x mutations, where x is out mutationsCountThreshold
		my $maxRF = $allInactivatingMutations[0];  ## Set maxRF as the position of the first inactivating mutation
		push(@allInactivatingMutations,1);

		for(my $j = 0; $j < scalar(@allInactivatingMutations)-1; $j++) {  ## Loop over the mutationsPosList to find the max intact Reading Frame
			my $rfProp = $allInactivatingMutations[$j+1] - $allInactivatingMutations[$j];
			$maxRF = $rfProp if ($rfProp > $maxRF);
		}
		
		$geneLossS1 = "T" if ($maxRF <= $maxRFThreshold);  ## Set this to true if the max un-affected Reading frame is less than our threshold
		print "The percent of intact Reading frame is $maxRF\n" if ($verbose);
	} else {
		print "We require the gene to have atleast $mutCountThreshold mutations, but the gene has only $nMut\n" if ($verbose);
	}
	
	
	if ($verbose) {
		print "The minNumberExonsHit ($mutExonsCountThreshold) criteria is NOT fulfilled\n" if ($geneLossS2 eq "");
		print "The minNumberExonsHit ($mutExonsCountThreshold) criteria is fulfilled\n" if ($geneLossS2 eq "T"); 
	}
	
	my $geneLoss = "";
	$geneLoss = "T" if ($geneLossS1 eq "T" && $geneLossS2 eq "T");
	return $geneLoss;
}

sub getCompensatoryEvents_NoExonDel {  ## Returns a hash which contains compensatory events, the key being the ID of the frameshift mutation, the value is just set to T
				     ## Exon deletions can also compensate other exon deletions or frameshifts, for example a -1 at the end of an exon can be compensated by a downstream -65 
				     ## (assuming that the length of the downstram exon is 65 nt). However, any stretch of compensatory events that includes an exon deletion is NOT considered as compensatory event	
	my($ref2MutationsList,$species) = @_;

	my %exonDeletions = my %indelToExonHash = ();
	## Step1: Get a list of exons that are deleted.
	foreach my $line(@$ref2MutationsList) {
		$line =~s/\s+$//;
		my ($speciesInFile,$eventID,$exonNumber) = (split /\t/,$line)[0,2,4];
		next if ($species ne $speciesInFile);
		
		$exonDeletions{$exonNumber} = "T" if ($line =~/Entire exon deletion/);
		$indelToExonHash{$eventID} = $exonNumber;
	}
	
	### Step2 --> Get a list of compensatory events. Basically make stretches of comp events and not just pairs. For example if a compensatory event is 
	### a combination of more than 2 events, ID1 - ID2 and ID2 - ID3
	### this step will make one "composite" event --> ID1-ID2-ID3. If any of them happens to be an exon deletion, then this stretch is no longer considered as a compensatory event stretch
	my @compEventsList = ();
	
	foreach my $line(@$ref2MutationsList) {
		$line =~s/\s+$//;
		if ($line =~/Compensation/) {
			my $speciesInLine = (split /\t/,$line)[0];
			push(@compEventsList,$line) if ($speciesInLine eq $species);
		}
	}
	
	my %compEventsFinal = ();
	my %interMediateHash = ();      ## This hash stores all the events that are intermediate and should be ignored while printing.
	
	foreach my $line(@compEventsList) {
		$line =~s/\s+$//;
		my($acc,$eventFirst,$event2) = (split /\t/,$line)[1,4,5];
		
		my %allEvents = ();	## This hash stores all the events.
		$allEvents{$eventFirst} = "T";

		if (! exists $interMediateHash{$eventFirst}) {
			
			my $eventCompPrint = $event2;
			$interMediateHash{$eventCompPrint} = "T";	
			$allEvents{$eventCompPrint} = "T";
	
			my $eventComp = getCompEvent($event2,$ref2MutationsList,$species);		## See if this event is the first event for any other comp event. For e.g. if one Comp Event is ID1 ID2 and another is ID2 ID3, ID2 is the second event in 
															## first line, and first event in the second line.
			if ($eventComp ne "") {								## If there is something
				$interMediateHash{$eventComp} = "T";
				$allEvents{$eventComp} = "T";
				
				do {		### Continue to do so until your compEventFirst is not the second event anymore.
				
					$eventCompPrint = $eventComp;
					$eventComp = getCompEvent($eventComp,$ref2MutationsList,$species); 

					$interMediateHash{$eventComp} = "T" if ($eventComp ne "");
					$allEvents{$eventComp} = "T" if ($eventComp ne "");
				}
				while ($eventComp ne "");
			}		
			
			my @compEventsPrint = keys(%allEvents);
			
			my $flag = "F";
			foreach my $compEvent(@compEventsPrint) {		## A compensatory event stretch is no good if one of the events happen to be an exon deletion
														## In such cases, none of the events in the comp event stretch will be considered as a comp event
				my $exon = $indelToExonHash{$compEvent};
				if (exists $exonDeletions{$exon}){
					$flag = "T";
					last;
				}
			}

			if ($flag eq "F") {
				foreach my $compEvent( keys (%allEvents)) {
					$compEventsFinal{$compEvent} = "T";
				}
			}
		}
	}
	return \%compEventsFinal;
}

sub getCompEvent {
	my ($event,$ref2MutationsList,$species) = @_;
	my $eventComp = "";
	
	foreach my $line(@$ref2MutationsList) {
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
		next if ($tmp[0] ne $species);
		
		$eventComp = $tmp[5] if ($tmp[4] eq $event);
	}
	return $eventComp;
}
