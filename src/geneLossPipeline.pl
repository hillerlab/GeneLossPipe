#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use FileHandle;

use lib "$ENV{'GeneLossPipeCode'}";
use MyFunctions_VS;
use SeqContextPlusHyperlink;
use useful_functions_VS;
use SpliceSiteFunctionsVS;	## This module contains all the functions related to the SpliceSite mutation and Alternative SpliceSite detection functionality
use IndelFunctions;			## This module contains all the functions related to Indel and also to InFrameStop codon detection.

my $verbose = 0;
my $excludeFPI = my $edeFlag = 0;

sub usage
{
	die "usage: $0 [-v|verbose] -fmt (optional, if set to 'lst', an input file needs to be specified) -input Accession/AccessionListFile -species Species (species list, comma separated) -out Output_Directory  -alignment Alignment(46way|100way|60way etc) -ref ReferenceSpecies -tabDir DirectoryWithTabularFiles -pMammals -EDE (Exclude Dubious Events, default is report all events) -excludeFPI (exclude FramePreserving events) -all (Run 3 modules - Indel, InframeStopCodon and SpliceSite) -indel (If only Indel) -ifsc (If only InFrameStop Codon) -ssm (If only SpliceSiteMutation module) -new (if using realigned mafs) -old (if using old mafs)\n";
}

my $fmt = my $input = my $speciesList = my $OutDir = my $alignment = my $reference = my $data_path = my $pMammals = "";
my $allEvents = my $indelPrint = my $ifscPrint = my $ssmPrint = 0;
my $old = my $new = 0;

GetOptions ("v|verbose"  => \$verbose, "fmt:s" =>\$fmt, "input=s" => \$input, "species:s" => \$speciesList, "out=s" => \$OutDir, "alignment:s" => \$alignment, "ref=s" => \$reference, "tabDir=s" => \$data_path, "pMammals:s" => \$pMammals, "EDE" => \$edeFlag, "excludeFPI" => \$excludeFPI, "indel" => \$indelPrint, "ifsc" => \$ifscPrint, "ssm" => \$ssmPrint, "all" => \$allEvents, "old" => \$old, "new" => \$new) || usage();
usage() if  ($input eq "");

## Check and die with appropriate error message if all the input parameters are not specified appropriately
die "Reference species not specified\n" if ($reference eq "");
die "Input alignment not specified\n" if ($alignment eq "");
die "Format is list but no input list of accessions provided\n" if ( ($fmt eq "lst") && (($input eq "") || (! -e $input)) );
die "Directory with tabular files is not specified/does not exist\n" if ( ($data_path eq "") || (! -d $data_path) );
die "pMammals and species list cannot be specified together\n" if ( ($pMammals eq "T") && ($speciesList ne "") ) ;
die "No output directory is specified\n" if ($OutDir eq "");
die "What version of the mafs you have used? Please specify either 'new' or 'old'\n" if (! $new && ! $old);

`mkdir -p $OutDir`;

## Now get a list of accessions and a list of species where you want to check for the presence of different inactivating mutations
my($ref2AccList,$ref2SpeciesList) = ParseParameters($fmt,$input,$pMammals,$speciesList,$reference,$alignment);
my @AccList = @$ref2AccList;
my @placental_mammals = @$ref2SpeciesList;

foreach my $acc(@AccList)
{
	$acc =~s/\s+$//;
	my $glpBDB = "$OutDir/geneLoss.$acc.BDB";
	next if (-e $glpBDB);
	
	my $tabFileAll = "$data_path/$acc.txt";
	if (! -e $tabFileAll){
		print "Skipping, this gene has no tabFiles dir --> i.e. $tabFileAll is missing\n";
		next;
	}else{ ### Split the bigTabFile into per species tabFiles
		splitBigTabFile($data_path,$acc,$tabFileAll);
	}
	
	## Do some checks to ensure that the tabFile is complete for each species;
	my @exonTermsList = ();
	my $noe = get_noe($acc,$reference,$data_path);
	my $i = 1;
	while ($i <= $noe)
	{
		my $term = "Coding_exon\_$i";
		push(@exonTermsList,$term);
		$i++;
	}

	foreach my $species(@placental_mammals)
	{
		my $speciesFile = "$data_path/$acc/$species.txt";
		next if (! -e $speciesFile);
		my %exonTermHash = ();
			
		open(FIR,$speciesFile) || die "Error opening species file '$speciesFile'\n";
		while (my $line = <FIR>)
		{
			$line =~s/\s+$//;
			my $term = (split /\t/,$line)[5];
			if ($term =~/Coding_exon/){
				$exonTermHash{$term} = "T";
			}
		}
		close FIR;
		
		foreach my $ex(@exonTermsList)
		{
			die "The exon '$ex' is missing from the species file '$speciesFile'\n" if (! exists $exonTermHash{$ex});
		}
	}		
	## Tests done, now run the different modules of the GeneLossPipe
	
	my $strandRef = "";
	## Get strandRef from the reference file
	my $refFile = "$data_path/$acc/$reference.txt";
	open(FIR,$refFile) || die "Error opening $refFile\n";
	while (my $line = <FIR>)
	{
		if ($line =~/Coding_exon/){
			my $strand = (split /\t/,$line)[4];
			$strandRef = $strand;
			last;
		}
	}
	close FIR;
	
	my $outFile = "$OutDir/$acc.txt";		## This is the output file to which the pipeline output is printed
	open(FOUT,">$outFile") || die "Error!! Cannot write to the outFile '$outFile'\n";

	foreach my $species(@placental_mammals)
	{	
		my $tabFile = "$data_path/$acc/$species.txt";
		if (-e $tabFile){		## Check if the tabular file exists
			print "Running the pipeline for '$acc' .. '$species'\n" if ($verbose);
			my $refFile = "$data_path/$acc/$reference.txt";
	
			## Most of the processing that is done here is related to Indel and InFrameStopCodon module, since for SpliceSiteMutation, we do not need the CDS.
			my @GeneArray = ();		## An array that practically stores the entire tabular file barring the non-coding/splice site regions
			my $noe = get_noe($acc,$reference,$data_path);		## Get the number of exons
		
			my $i = my $h = 0;
			my $exonCt = 1;
			my %MapToRefGenomeHash = my %RefGenomeHash = my %ExonLength = ();
			my %HashToGappedPos = ();
			my $GeneSequence = "";
						
			while ($exonCt <= $noe)
			{
				my $Chunk = get_exonic_chunk($acc,$species,$exonCt,$reference,$data_path);	## Get the exonic chunk
				my $exonLen = 0;
				
				my @ExonChunk = split(/\n/,$Chunk);
				foreach my $line(@ExonChunk)
				{
					my @tmp = split(/\t/,$line);
					my ($baseRef,$baseQuery) = (split /\t/,$line)[2,10];
					$GeneArray[$i] = $line;							## Store all the useful lines in an array @GeneArray.
					
					if ($baseRef ne "-"){
						$exonLen++;
						my $string = "$exonCt\_$exonLen";
						
						$MapToRefGenomeHash{$i} = $string;			## This is a useful hash. For every line that is a base in the Reference species, the key is the index of GeneArray
						$RefGenomeHash{$i} = $tmp[2];				## the value is exonNumber_exonIndex. Useful to indicate the position of insertions/deletions in terms of the exonic positions
					}												## of the Reference species
					
					if ($baseQuery ne "-"){
						$GeneSequence = $GeneSequence.$tmp[13];
						$HashToGappedPos{$h} = $i;
						$h++;
					}
					
					$i++;
				}
				$ExonLength{$exonCt} = $exonLen;
				$exonCt++;
			}
			
			### Once everything has been dumped into @GeneArray, look for insertions, deletions and missing data.
			my @posInsertions = my @posDeletions = my @posMissingData = ();
			my $j = 0;
			my $pI = my $pD = my $pMD = 0;

			foreach my $line(@GeneArray)
			{
				my ($baseRef,$baseQuery) = (split /\t/,$line)[2,10];
				
				if ($baseRef eq "-" && $baseQuery ne "-"){			## There is sequence in the query species but not in the reference, hence an insertion
					$posInsertions[$pI] = $j;
					$pI++;
				}elsif ( $baseRef ne "-" && $baseQuery eq "-" ){		## There is sequence in the reference species but not in query species, hence a deletion
					$posDeletions[$pD] = $j;
					$pD++;
				}elsif  ($baseQuery eq "?"){ 								## Simply missing data	
					$posMissingData[$pMD] = $j;
					$pMD++;
				}
				
				$j++;
			}
			
			### Now club these individual events to make one event. For instance deletions are recorded as a 1 bp deletion at pos6, 1 bp deletion at pos7 and 1 bp deletion at pos8.
			### This operation clubs these events. For the above example, there is now one 3bp deletion at pos6.
			
			my $Ref2AllIns = my $Ref2AllDel = my $Ref2AllMD = "";
			
			$Ref2AllIns = ClubIndels(\@posInsertions)  if ($pI > 0);   ## Returns a Reference to all insertions hash
			$Ref2AllDel = ClubIndels(\@posDeletions)   if ($pD > 0);    ## Returns a Reference to all deletions hash
			$Ref2AllMD  = ClubIndels(\@posMissingData) if ($pMD > 0);   ## Returns a Reference to all missing data positions hash
			
			### Separate events (deletions and missing data) based on their exon positions. For instance a 3 bp deletion could actually be at the end of 1st exon.
			### Which means ==> 3bp deletion = 2bp deletion at 1st exon, 1 bp deletion at 2nd exon
			
			if ($pMD > 0)	### If there are indeed some missing data stretches, print them here. And that's it.
			{
				$Ref2AllMD  = checkDeletionSpan($Ref2AllMD,\@GeneArray);
				my %missingDataHash = %$Ref2AllMD;
				
				my $EventCt = 1;
				
				foreach my $event( sort {$a <=> $b} keys(%missingDataHash))
				{
					my $length = $missingDataHash{$event};
					my ($posRef,$RefGenomeCds,$RefGenomeParam) = GetReferencePos($event,"MissingData",\%MapToRefGenomeHash,\@GeneArray,$length,\%RefGenomeHash,$tabFile);
					
					my ($exonNumber,$exonPos) = (split /_/,$posRef)[0,1];											## Get the exon number and the position from where the deletion/insertion begins
					
					if ($length == $ExonLength{$exonNumber}){
						print FOUT "$species\t$acc\tMD$EventCt\tEntire exon sequence missing\t$exonNumber\t$length\t-\n";
					}else{
						print FOUT "$species\t$acc\tMD$EventCt\tMissing sequence\t$exonNumber\t$exonPos\t$length\t-\n";
					}

					$EventCt++;
				}
			}
			
			my $Nevents = $pI + $pD;								## Nevents = Number of events --> All insertions plus deletions.
			## All insertions business from now on:
			my $Ref2FDIns = "";
			my $FDIns = 0;
			my %InsertionsHash = ();
			my %EventsExcludedInsertions = ();		### Events excluded because the quality is questionable
			
			if ($pI > 0){
				my %tmpHash = %$Ref2AllIns;
				$pI = 0;
				
				foreach my $keys(keys(%tmpHash))
				{
					$InsertionsHash{$keys} = $tmpHash{$keys} if ($keys != 0);		### Ignore insertions which are upstream of the ATG in the reference gene, as an example check NM_003068/mm9
					$pI++;
				}
				
				if ($edeFlag){	### Do this only if the ede flag is set to true .. "ede" --> Exclude Dubious Events
				
					my ($ref2FilteredInsertions,$ref2ExcludedInsertions) = checkFlankingNts(\%InsertionsHash,\@GeneArray,"Insertion");		## See details of checkFlankingNts in the function itself.
					%InsertionsHash = %$ref2FilteredInsertions;
					
					my %tmpHashEE = %$ref2ExcludedInsertions;
					
					foreach my $keys(keys(%tmpHashEE))
					{
						my $length = $tmpHashEE{$keys};
						
						my $start = $keys;
						my $stop = $keys + $length -1;
						
						while ($start <= $stop)
						{
							$EventsExcludedInsertions{$start} = "T";
							$start++;
						}
					}
				}
				($Ref2FDIns,$FDIns) = checkFDEvents(\%InsertionsHash); ## Returns a reference to frame-disrupting insertions as well as count of frame-disrupting insertions.
			}
			
			my %AllIns = my @AllInsList = ();
			%AllIns = %InsertionsHash if ($pI > 0);
			@AllInsList = keys(%AllIns);
		
### At the end of the Insertion business, I have the following -> %InsertionsHash where the key is the position of the start of the insertion, while the value is the length of the insertion
### Also I have an array --> @AllInsList whose values are the positions of the insertions. Another reference to hash $Ref2FDIns --> same as %InsertionsHash except that
### it has only frame-disrupting insertions, a variable $FDIns stores the number of frame-disrupting insertions
			
			## All deletions business --> more or less the same processing as for insertions above
			my $Ref2FDDel = "";
			my $FDDel = 0;
			my %EventsExcludedDeletions = ();		### Events excluded because the quality is questionable
			
			my $Ref2AllDelPrint = "";	## This is the real deal. if the edeFlag is set to true, this hash will contain all the events that remain after filtering.
			
			if ($pD > 0){
				if ($edeFlag){
					$Ref2AllDelPrint = checkDeletionSpan($Ref2AllDel,\@GeneArray); ## Checks which deletions span over more than one exon --> look at the details in the sub-routine "checkDeletionSpan".
					
					## Get deletions that are entire exon deletions, they should not be excluded as dubious events simply because they have a N on either side (note that this N lies somewhere in the upstream or the downstream exon)
					my %deletionDataHash = %$Ref2AllDelPrint;
					my @exonDeletionEvents = ();
					foreach my $event( sort {$a <=> $b} keys(%deletionDataHash))
					{
						my $length = $deletionDataHash{$event};
						my ($posRef,$RefGenomeCds,$RefGenomeParam) = GetReferencePos($event,"Deletion",\%MapToRefGenomeHash,\@GeneArray,$length,\%RefGenomeHash,$tabFile);
						
						my ($exonNumber,$exonPos) = (split /_/,$posRef)[0,1];											## Get the exon number and the position from where the deletion/insertion begins
						push(@exonDeletionEvents,$event) if ($length == $ExonLength{$exonNumber});
					}

					my ($ref2FilteredDeletions,$ref2ExcludedDeletions) = checkFlankingNts($Ref2AllDelPrint,\@GeneArray,"Deletion");
					### Now remove exonDeletion events from excludedDeletions and put them in filtered deletions
					my %filteredDeletions = %$ref2FilteredDeletions;
					my %excludedDeletions = %$ref2ExcludedDeletions;
					
					foreach my $event(@exonDeletionEvents)
					{
						if (exists $excludedDeletions{$event}){  ## If an exonDeletion was filtered out, put it back in filteredDeletions hash and remove from excludedDeletions hash
							my $length = $excludedDeletions{$event};
							$filteredDeletions{$event} = $length;
							delete $excludedDeletions{$event};
						}
					}
					
					$ref2FilteredDeletions = \%filteredDeletions;
					$ref2ExcludedDeletions = \%excludedDeletions;

					$Ref2AllDelPrint = $ref2FilteredDeletions;	### Here you go. Only the filtered deletions are included here.
					($Ref2FDDel,$FDDel) = checkFDEvents($ref2FilteredDeletions);	## Returns a reference to frame-disrupting deletions as well as count of frame-disrupting deletions
								
					my %tmpHashEE = %$ref2ExcludedDeletions;
					foreach my $keys(keys(%tmpHashEE))
					{
						my $length = $tmpHashEE{$keys};
						
						my $start = $keys;
						my $stop = $keys + $length -1;
						
						while ($start <= $stop)
						{
							$EventsExcludedDeletions{$start} = "T";
							$start++;
						}
					}
				}
				
				else{		### Otherwise all the deletions get into this hash.
					$Ref2AllDelPrint = checkDeletionSpan($Ref2AllDel,\@GeneArray); ## Checks which deletions span over more than one exon --> look at the details in the sub-routine "checkDeletionSpan".
					($Ref2FDDel,$FDDel) = checkFDEvents($Ref2AllDelPrint);	## Returns a reference to frame-disrupting deletions as well as count of frame-disrupting deletions
				}
			}
			
			my %AllDel = my @AllDelList = ();
			%AllDel = %$Ref2AllDelPrint if ($pD > 0);
			@AllDelList = keys(%AllDel);
			
### End of deletion business. Now I have a reference to all deletions hash --> $Ref2AllDel. This is further processed by the checkDeletionSpan function.
### And ultimately we have a hash called %AllDel. A list of deletions stored in an array called @AllDelList 
### Also a reference to frame disrupting deletions -->	$Ref2FDDel
### Number of framedisrupting deletions --> $FDDel

### Now I begin to separate the frame-disrupting indels from the frame-preserving ones.
				
			my @eventsAll = (@AllInsList, @AllDelList); ## Basically all events are combined, whether they are frame-disrupting or not.
			my @eventsAllSorted = sort {$a <=> $b} @eventsAll; ## and now sorted
			
			my %FDInsertions = %$Ref2FDIns if ($FDIns > 0);	### For frame disrupting indels
			my %FDDeletions  = %$Ref2FDDel if ($FDDel > 0);

			my @keysIns = keys(%FDInsertions);
			my @keysDel = keys(%FDDeletions);
					
			my @IndelsAllFD = (@keysIns,@keysDel);		### Combine all events which are FrameDisrupting
			my @IndelsAllSortedFD = sort {$a <=> $b} @IndelsAllFD;
			
			### Get a hash where the key is the position of the event (insertion|deletion) and the value is the length. Deletions have a length that is in minus
			my %FDEvents = ();
			foreach my $events(@IndelsAllFD)
			{
				my $length = "";
				
				if ( exists $FDInsertions{$events} ){
					$length = $FDInsertions{$events};
				}elsif ( exists $FDDeletions{$events} ){
					$length = 0 - $FDDeletions{$events};
				}
				
				$FDEvents{$events} = $length;
			}
			
			### Get a list of exons which are completely deleted.
			my %exonsDeleted = ();
			foreach my $event(keys(%AllDel))			## Now loop through the entire list of events
			{
				my $nature = "Deletion";
				my $length = $AllDel{$event};
				
				my ($posRef,$RefGenomeCds,$RefGenomeParam) = GetReferencePos($event,$nature,\%MapToRefGenomeHash,\@GeneArray,$length,\%RefGenomeHash,$tabFile);		## Get the exon number and the position in the exon from this event happens (everything is relative to the Reference)
				
				my ($exonNumber,$exonPos) = (split /_/,$posRef)[0,1];											## Get the exon number and the position from where the deletion/insertion begins
				$exonsDeleted{$exonNumber} = "T"  if ($length == $ExonLength{$exonNumber}); 	## If the event is deletion and its length is equal to the entire exon, mark it as "entire exon deletion"	
			}
			
			## Get a list of insertions (note that these are frame-preserving insertions that bring a stop codon, it is really important to catch them with the RealignedMafs) which bring in a stop codon:
			my %insertionStopCodon = ();
			foreach my $insertions(@AllInsList)
			{
				if (! exists $FDInsertions{$insertions}){	## It should be an insertion that does not disrupt the reading frame: thats it
					my $length = $InsertionsHash{$insertions};
					my $start = $insertions;
					my $stop  = $insertions + $length;
					
					my $rfStart = (split /\t/,$GeneArray[$start])[15];  ## The codon position at the beginning of the insert
					my $rfStop =  (split /\t/,$GeneArray[$stop])[15];   ## The codon position at the end of the insert
					
					## Why am I doing this
					## Ref		ATC------GTC
					## Query	ATCCGTTGAGTC    a 6bp insert in the query that brings in a stop codon, fair enough
					## Another case		Ref AT------CGTC
					## 		   Query    ATCCGTTGAGTC    a 6bp insert where the sequence of insert (CCGTTG) does not contain stop codon, but if you translate the sequence ATC-CGT-TGA-GTC, it contains a stop codon
					## So the code below ensures that I look for stop codons, not only in the sequence of insert but also pad it up with nucleotides up and down so that I am looking at a stretch that starts at C1 and ends up at C3
					
					if ($rfStart ne "C1"){  ## Go upstream from the start position until you find C1
						for(my $k = $start; $k > 0; $k--)
						{
							my $seqLine = $GeneArray[$k];
							my($base,$codonPosition) = (split /\t/,$seqLine)[13,15];
							if ($base ne "-" && $base ne "?" && $codonPosition eq "C1"){
								$start = $k;
								last;	
							}
						}
					}
					
					if ($rfStop ne "C3"){  ## Go downstream from the stop position until you find C3
						for(my $k = $stop; $k < scalar(@GeneArray); $k++)
						{
							my $seqLine = $GeneArray[$k];
							my($base,$codonPosition) = (split /\t/,$seqLine)[13,15];
							if ($base ne "-" && $base ne "?" && $codonPosition eq "C3"){
								$stop = $k;
								last;	
							}
						}
					}
					
					## Now check if the insert has a stop codon
					my $seq = "";
					while ($start < $stop)
					{
						my $line = $GeneArray[$start];
						my $nt = (split /\t/,$line)[13];
						$seq = $seq.$nt;
						$start++;
					}
					
					my $status = find_stop_codon($seq,"first");
					$insertionStopCodon{$insertions} = "T" if ($status ne "");
				}
			}
			
			my %CompensationHash = ();	 ## This has contains all the information about the events that can compensate for each other.
			my @ListExcludeRegions = (); ## This array contains the list of compensatory regions that are useful when looking for In-Frame stop codons. 
			my $ler = 0;

			if ($FDIns >= 1 || $FDDel >= 1){		## Number of insertions|deletions > 1 OR Atleast 1 insertion and 1 deletion
				my $Ref2CompensationHash = GetCompensatedEvents(\@IndelsAllSortedFD,\@GeneArray,\%FDEvents,\%MapToRefGenomeHash,\%RefGenomeHash,$tabFile,\%EventsExcludedInsertions,\%EventsExcludedDeletions);
				%CompensationHash = %$Ref2CompensationHash;
				
				foreach my $keys (sort {$a <=> $b} keys(%CompensationHash))
				{
					my $compIndel = $CompensationHash{$keys};
					
					my @tmpEvents = split(",",$compIndel);
					$ListExcludeRegions[$ler] = "$keys-$tmpEvents[$#tmpEvents]\n";
					$ler++;
				}
			}

####################################################################################
####  Now I begin to print mutations depending on which modules I have called  #####
####################################################################################

			if ( ($indelPrint) || ($allEvents) ){
				my %printedCDS = ();
				print "Printing indels for '$acc' --> '$species'\n" if ($verbose);
		
				if ($Nevents > 0){	
					my %EventCtHash = ();
					my $EventCt = 1;
					my %SpecialHash = my %SpecialHashExists = ();
					
					foreach my $event(@eventsAllSorted)			## Now loop through the entire list of events
					{
						my($nature,$length) = GetNatureLength($event,\%AllIns,\%AllDel);
						my $stopIntroducingIns = "";
						$stopIntroducingIns = "T" if (exists $insertionStopCodon{$event});
						
						my ($posRef,$RefGenomeCds,$RefGenomeParam) = GetReferencePos($event,$nature,\%MapToRefGenomeHash,\@GeneArray,$length,\%RefGenomeHash,$tabFile);		## Get the exon number and the position in the exon from this event happens (everything is relative to the Reference)
						
						my ($exonNumber,$exonPos) = (split /_/,$posRef)[0,1];				## Get the exon number and the position from where the deletion/insertion begins
						
						### Get the information about the reference. Useful for getting the link to the genome browser
						
						my ($chr,$strand) = (split /\t/,$RefGenomeParam)[3,4];
						my $EventType = "$nature\_$length";

						my $link   = Link_SeqContext($tabFile,$RefGenomeCds,$strand,$EventType);
		
						if ($nature eq "Deletion"){	## Print all entire exon deletions, no matter what their length is
							if  ($length == $ExonLength{$exonNumber}){ 	## If the event is deletion and its length is equal to the entire exon, mark it as "entire exon deletion"
								$nature = "Entire exon deletion";
								print FOUT "$species\t$acc\tID$EventCt\t$nature\t$exonNumber\t$length\t$link\t$RefGenomeCds\n";		## No need to put exon position here, it is obviously 1
								$printedCDS{$RefGenomeCds} = "T";	### Dump this cds in the printedCDS hash.
							}
						}
						
						if (! $excludeFPI){		## Print everything if the excludeFramePreservingIndels is set to "F"
							print FOUT "$species\t$acc\tID$EventCt\t$nature\t$exonNumber\t$exonPos\t$length\t$link\t$RefGenomeCds\n" if (! exists $printedCDS{$RefGenomeCds});
							## Also print (frame-preserving) insertions that bring in a stop codon:
							print FOUT "$species\t$acc\tID$EventCt-StopCodon\t$nature\t$exonNumber\t$exonPos\t$length\t$link\t$RefGenomeCds\n" if ($stopIntroducingIns eq "T");
						}else{
							if ($length%3 != 0){  	## Even if the flag is set to true, print everything that is not frame preserving 
								print FOUT "$species\t$acc\tID$EventCt\t$nature\t$exonNumber\t$exonPos\t$length\t$link\t$RefGenomeCds\n" if (! exists $printedCDS{$RefGenomeCds});
							}elsif ($length > 20){	## And print those frame preserving indels that are longer than 20 nt.
								print FOUT "$species\t$acc\tID$EventCt\t$nature\t$exonNumber\t$exonPos\t$length\t$link\t$RefGenomeCds\n" if (! exists $printedCDS{$RefGenomeCds});
							}
							## Also print (frame-preserving) insertions that bring in a stop codon:
							print FOUT "$species\t$acc\tID$EventCt-StopCodon\t$nature\t$exonNumber\t$exonPos\t$length\t$link\t$RefGenomeCds\n" if ($stopIntroducingIns eq "T");
						}
						
						$EventCtHash{$event} = "ID$EventCt";			## Remember the key is the original event pos, the value is the new count.
						$EventCt++;
					}		
					
					### At this point, all the indels have been printed. Now let us print the information about events which can compensate for each other.
					foreach my $event(@eventsAllSorted)
					{
						if (exists $CompensationHash{$event}){   
							my $EventNumber = $EventCtHash{$event};		## get the event number for the key
							my $compEvent = $CompensationHash{$event};	## get the compensatory event(s) from the hash.
							
							my @tmp = split(/,/,$compEvent);
							my @tmpPrint = ();
							
							my $k = 0;
							foreach my $cIndel(@tmp)
							{
								$tmpPrint[$k] = $EventCtHash{$cIndel};
								$k++;
							}
							
							unshift(@tmpPrint,$EventNumber);		## Now put all the events which are compensatory in nature in the array called @tmpPrint
							
							for(my $a = 0; $a < (scalar(@tmpPrint) -1); $a++)	## If the events are ID1 (index 0), ID2 (index 1) and ID3 (index 2), loop from 0 to 1 
							{													## First print ID1 and ID2, and then ID2 and ID3
								print FOUT "$species\t$acc\t-\tCompensation\t$tmpPrint[$a]\t$tmpPrint[$a+1]\n";
							}
						}
					}
				}
			}
			
			if ( ($ifscPrint) || ($allEvents) ){
				print "Printing inframe stop codons for $acc --> $species\n" if ($verbose);
				my %StopList = ("TAA"  => "T", "TAG"  => "T", "TGA"  => "T");
				
				my $GeneSequenceCopy = $GeneSequence;
				$GeneSequenceCopy =~s/-//g;
				$GeneSequenceCopy =~s/\?//g;

				if (length($GeneSequenceCopy) >= 3){					## If the gene sequence is nothing but a combination of "-"s and "?"s, or it has barely 2 nucleotides, do nothing
					
					my $Ref2SCL = find_stop_codon_all($GeneSequence);				## SCL -> StopCodonList
					my @StopCodonList = @$Ref2SCL;									## All stop codons.
					my $eventCt = 1;	

					foreach my $pos(@StopCodonList)
					{
						$pos = $HashToGappedPos{$pos};

						my $FirstLine  = $GeneArray[$pos];
						my $SecondLine = $GeneArray[$pos+1];
						my $ThirdLine  = $GeneArray[$pos+2];

						my @Line1 = split(/\t/,$FirstLine);
						my $refBase1st = $Line1[7];
						
						if ($Line1[6] eq "C1"){								## The first condition. The "found" stop codon should align to the first codon position in the Reference gene.
							my @Line2 = split(/\t/,$SecondLine);
							my $refBase2nd = $Line2[7];
							my @Line3 = split(/\t/,$ThirdLine);
							my $refBase3rd = $Line3[7];
							
							my $flag = "T";
							foreach my $boundary(@ListExcludeRegions)		## Just check if the position of the stop codon is not within in the IndelCompensatoryRegion. Such codons do not qualify as in-frame stop codons
							{
								$boundary =~s/\s+$//;
								my @bd = split(/-/,$boundary);
								
								if ($pos >= $bd[0] && $pos <= $bd[1]){
									$flag = "F";
								}
							}
							
							if ($flag eq "T" && $refBase2nd ne "-"  && $refBase3rd ne "-"){ ## Also the "found" stop codon should align to a codon in the reference gene. 
																								 ## Like TGA aligning to TCA is fine. TGA -> T-CA or TGA -> TC-A would not work.
								my $refCodon   = $refBase1st.$refBase2nd.$refBase3rd;
								my $speciesCodon = $Line1[13].$Line2[13].$Line3[13];
						
								if ($speciesCodon !~/-/ && ! exists $StopList{$refCodon}){
									my ($exonNumber,$exonPos) = (split /_/,$MapToRefGenomeHash{$pos})[0,1];
									my $RefCds = $Line1[2];
									my $strandRef = $Line1[4];	
									my $link = Link_SeqContext($tabFile,$RefCds,$strandRef,"IFSC");								

									print FOUT "$species\t$acc\tIFSC$eventCt\tIn-frame stop codon\t$exonNumber\t$exonPos\t$refCodon->$speciesCodon\t$link\t$RefCds\n";
									$eventCt++;
								}
							}
						}
					}
				}
			}
			
			###### Splice Site module
			
			if ( ($ssmPrint) || ($allEvents) ){
				print "Printing splice site mutations for $acc <--> $species\n" if ($verbose);

				if ($noe > 1){		## The species has to have more than one exon, only then it can have splice sites.
					my @searchTerm_Species = ();		## Get search term list: for instance 1_D,2_A,2_D,3_A
					my $exonStart = 1;     my $exonEnd   = $noe;
					my $firstAcc  = "1_A"; my $lastDonor = "$exonEnd\_D";
					
					while ($exonStart <= $exonEnd)
					{
						if (! exists $exonsDeleted{$exonStart}){
							my $accTerm   = "$exonStart\_A";
							my $donorTerm = "$exonStart\_D";
							push(@searchTerm_Species,$accTerm);
							push(@searchTerm_Species,$donorTerm);
						}
						
						$exonStart++;	
					}
					@searchTerm_Species = grep{$_ ne $firstAcc} @searchTerm_Species;
					@searchTerm_Species = grep{$_ ne $lastDonor}  @searchTerm_Species;
						
					my $Ref2SSlist_Reference  = GetSSList_Ref($acc,\@searchTerm_Species,$reference,$data_path);  ## Get the list of splice sites for the reference gene.
					my %ss_list_Reference = %$Ref2SSlist_Reference;								  ## What you get is a reference to the hash. The hash contains key-value pairs where key is the junction (like 1_D) and the value is the nt/bases (GT).
					
					my $ss_list_ref_cds = GetSSCds($acc,\@searchTerm_Species,$reference,$data_path);	## Now get a list of splice site coordinates for each junction.
					my %ss_list_cds     = %$ss_list_ref_cds;		
					
					my $EventCt = 1;
						
					foreach my $junction(@searchTerm_Species)
					{
						my($splice_site,$refCds,$speciesCds,$strandRef)   = GetSpliceSite($junction,$tabFile);
						
						my ($exon_index,$junc_type) = (split /_/,$junction)[0,1];
						my $spliceSite_copy = $splice_site;
						$spliceSite_copy =~s/-//g;
						
						my %JuncHash =("A" => "Acceptor", "D" => "Donor");
						my $eventType = "SSM\_$JuncHash{$junc_type}";
						
						if ($junc_type eq "D"){ ## For donor splice site junctions
							## Ignore those situations where the splice site in the reference and the query do not match, because that does not imply a splice site mutation
							## Reference Splice Site = GT, Query = GC or vice-versa.
							
							my $ss_flag = "";
							$ss_flag = "T"  if ( ( ($ss_list_Reference{$junction} eq "GT") && ($splice_site eq "GC") ) || ( ($ss_list_Reference{$junction} eq "GC") && ($splice_site eq "GT") ) );
							$ss_flag = "T"  if ($ss_list_Reference{$junction} eq $splice_site); ### This fixes the non-canonical cases
							$ss_flag = "T"  if ($splice_site eq "GT" || $splice_site eq "GC");  ## Even when the reference has a weird splice site, this should not be reported as a mutation
							
							if ( ($ss_flag ne "T") && ($ss_list_Reference{$junction} ne $splice_site) && ($splice_site !~/N/) && ($splice_site !~/\?/) ){ ## Test 1)  if the splice site for this gene is different from the splice site of the human gene 2) if the splice site does not contain a base of poor quality.
								my $HumanLink = Link_SeqContext($tabFile,$refCds,$strandRef,$eventType);
								if ($new){  ## For the realigned mafs, just print this down as a mutation
									print FOUT "$species\t$acc\tSS$EventCt\tSplice site mutation\t$exon_index\t$JuncHash{$junc_type}\t$splice_site\t-\t-\t$HumanLink\t$refCds\n"; ## write to the mutations file.
									$EventCt++;		
								}elsif ($old){ ## For the older version, check if the splice site is G- or -G, if yes, extract one base pair downstream and check if then the splice site becomes GT/GC
								
									#my $nt_1Down = "";				   	
									#$nt_1Down = get_1nt($species,$acc,$junction,$reference,$data_path) if ( (length($spliceSite_copy) == 1) && ($spliceSite_copy eq "G") );  ### Example: Splice site in query is G- or -G;					### Extract 1nt downstream									
									
									#if ( ($nt_1Down ne "T") && ($nt_1Down ne "C") ) ## write to the mutations file. 			
									#{
										print FOUT "$species\t$acc\tSS$EventCt\tSplice site mutation\t$exon_index\t$JuncHash{$junc_type}\t$splice_site\t-\t-\t$HumanLink\t$refCds\n";
										$EventCt++;
									#}
								}
							}
						}
						else{ ### For Acceptor positions. Check if the splice sites do not match, and do not contain any ambigious character
						
							if ( ($ss_list_Reference{$junction} ne $splice_site) && ($splice_site !~/N/) && ($splice_site !~/\?/) && ($splice_site ne "AG") ){							
								my $HumanLink = Link_SeqContext($tabFile,$refCds,$strandRef,$eventType);
								if ($new){				
									print FOUT "$species\t$acc\tSS$EventCt\tSplice site mutation\t$exon_index\t$JuncHash{$junc_type}\t$splice_site\t-\t-\t$HumanLink\t$refCds\n"; ## write to the mutations file.
								
								}elsif ($old){
								
									#my $nt_1up = "";
									#$nt_1up = get_1nt($species,$acc,$junction,$reference,$data_path) if ( (length($spliceSite_copy) == 1) && ($spliceSite_copy eq "G") );  ### Example: Splice site in query is -G
								
									#if ($nt_1up ne "A")
									#{																						   	
										print FOUT "$species\t$acc\tSS$EventCt\tSplice site mutation\t$exon_index\t$JuncHash{$junc_type}\t$splice_site\t-\t-\t$HumanLink\t$refCds\n" ; ## Write to the mutations file.		
										$EventCt++;
									#} 		
								}
							} 
						}
					}#### end of the search term list loop
				}
			}
		}
	}
	
	close FOUT;
	
	`rm -rf $data_path/$acc`;
	my $nMutAll = `wc -l < $outFile`; chomp $nMutAll;
	if ($nMutAll > 0) {
		formatOutputFile($outFile,$strandRef);
		open(FIGL,$outFile) || die "Error opening geneLossPipe out file '$outFile'\n";
		my @allMutations = <FIGL>;
		close FIGL;
		my $allMutationsString = join("",@allMutations);
		my $outBDB = "$OutDir/geneLoss.$acc.BDB";
		secureWriteBDB($outBDB,$acc,$allMutationsString,1,1);
	}
	`rm -rf $outFile`;
}

### SubRoutine that formats the output file so that the events are listed in the 5' -> 3' order along the gene sequence
sub formatOutputFile
{
	my ($outFile,$strandRef) = @_;  ## Note that this file contains mutations for all the species, the sorting and subsequent printing has to be done per-species
	
	## Get a list of species present in the mutations file;
	my %speciesWithMutations = ();
	open(FIM,$outFile) || die "Error opening the mutations file '$outFile' in formatOutputFile function\n";
	while (my $line = <FIM>)
	{
		my $species = (split /\t/,$line)[0];
		$speciesWithMutations{$species} = "T";
	}
	close FIM;
	
	my @speciesList = keys(%speciesWithMutations);	
	my $outFileFormatted = create_temp("f");
	open(FOM,">$outFileFormatted") || die "Error writing to outFile formatted in formatOutputFile function\n";
	## Now loop over all the species
	foreach my $species(@speciesList)
	{
		my @speciesMutationsList = (); ## this array contains mutations for each species in the file
		open(FIM,$outFile);
		while (my $line = <FIM>)
		{
			my $speciesInFile = (split /\t/,$line)[0];
			push(@speciesMutationsList,$line) if ($species eq $speciesInFile);
		}
		close FIM;
			
		my %hashCds = ();		### Make this as a 2 dimensional hash. Sometimes the same reference genome cds could have more than 1 event.
		my %compEvents = ();	### For example, NM_206837, mm9 assembly. 5th exon - Mutation at acceptor position, Insertion at the beginning of 5th exon as well.
		my $c = my $m = 1;
		my $j = 0;
		my %missingSeq = ();

		foreach my $line(@speciesMutationsList)
		{
			$line =~s/\s+$//;
			my @tmp = split(/\t/,$line);
		
			my $refGenomeCds = $tmp[$#tmp];
	
			if ($line =~/Compensation/){
				$compEvents{$c} = $line;
				$c++;
			}elsif ($line =~/missing/i){
				$missingSeq{$m} = $line;
				$m++;
			}else{
				$hashCds{$refGenomeCds}{$j} = $line;
				$j++;
			}
		}

		my @cdsSorted = sort {$a <=> $b} ( keys(%hashCds) );
		@cdsSorted = reverse(@cdsSorted) if ($strandRef eq "-");

		foreach my $cds(@cdsSorted)   ## Print mutations
		{
			foreach my $number( keys %{$hashCds{$cds}} )
			{
				my $line = $hashCds{$cds}{$number};
				print FOM "$line\n";
			}
		}
		
		if ($c > 1){  ## print compensatory events, if any
			foreach my $comp( sort {$a <=> $b} (keys (%compEvents) ) )
			{
				my $line = $compEvents{$comp};
				print FOM "$line\n";
			}
		}

		if ($m > 1){    ## print missing sequence events, if any
			foreach my $miss( sort {$a <=> $b} (keys (%missingSeq) ) )
			{
				my $line = $missingSeq{$miss};
				print FOM "$line\n";
			}
		}
	}
	close FOM;
	`mv $outFileFormatted $outFile`;
}

sub splitBigTabFile
{
	my($data_path,$acc,$tabFileAll) = @_;
		
	`mkdir -p $data_path/$acc`;
	my %speciesListAll = ();
	open(FIT,$tabFileAll);
	while (my $line = <FIT>)
	{
		$line =~s/\s+$//;
		if ($line ne "FileDone")
		{
			my $species = (split /\t/,$line)[9];
			$speciesListAll{$species} = "T";
		}
	}
	close FIT;

	my %sp2File = ();	## Open separate filehandles for each species
	foreach my $species(keys(%speciesListAll))
	{
		my $fh = FileHandle->new();       
		$sp2File{$species}=$fh;
		$fh->open(">$data_path/$acc/$species.txt");
	}
		
	open(FIT,$tabFileAll);  ## write to each species file via individual fileHandles
	while (my $line = <FIT>)
	{
		$line =~s/\s+$//;
		if ($line ne "FileDone")
		{
			my $species = (split /\t/,$line)[9];
			my $fh = $sp2File{$species};
			print $fh "$line\n";
		}
	}
	close FIT;
		
	foreach my $species(keys(%speciesListAll))  ## close fileHandles
	{
		my $fh = $sp2File{$species};
		$fh->close;
	}
}
