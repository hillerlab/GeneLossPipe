#!/sw/bin/perl

### Virag Sharma, 
### This module contains a set of functions associated with the IndelDetection module
### Dependencies -> Uses functions from MH's MyFunctions.pm, MyFunctions_VS and useful_functions_VS

package IndelFunctions;
use strict;
use warnings;
use Exporter;

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;

use lib "$ENV{'GeneLossPipeCode'}";
use MyFunctions_VS;
use useful_functions_VS;

our @ISA = ('Exporter');
our @EXPORT = qw(ClubIndels GetNatureLength GetReferencePos checkDeletionSpan checkFDEvents GetCompensatedEvents checkFlankingNts);

our $genomePath = $ENV{'genomePath'};


########################################################################################################################################################
#### This function clubs the events - whether insertions or deletions together based on their adjacency.											####
#### For e.g. if the input is an array with the following values 5,6,7,8,9,10,14,15,18,21															####
#### the return will be a hash where %hash = (5 => 6, 14 =>2, 18 => 1, 21 => 1)  (key is 5, the value is 6 (length of this stretch from 5-10))		####
########################################################################################################################################################

sub ClubIndels                                 		
{                                                   
	my $Ref2List = shift;                           
	my @tmp_array = @$Ref2List;
	
	my @tmp_dump_array = ();
	my %GapsPosLength = ();
	my $t = 0;

	my $start_pt = my $pos_start = my $gap_length = "";
	
	for (my $i = 0; $i < @tmp_array; $i++)
	{
		if ($i > 0)
		{
			if ($tmp_array[$i] != $start_pt+1)                      ## Now test if the current value and the $start_pt are adjacent
			{
				$gap_length = scalar(@tmp_dump_array);              ## If not, get the length of the indel (basically the number of elements present in the @tmp_dump_array)
				$pos_start  = $tmp_dump_array[0];                   ## Get the key for the hash $pos_start which is the first element of @tmp_dump_array
				$GapsPosLength{$pos_start} = $gap_length;          ## Now put these values in the hash.

				@tmp_dump_array = ();
				$t = 0;
			}
		}

		$tmp_dump_array[$t] = $tmp_array[$i];                     ## Dump every element into the @tmp_dump_array first.
		$t++;

		$start_pt = $tmp_array[$i];                                       ## Initialize the variable $start_pt with the ith index of the input array.
	}

	$gap_length = scalar(@tmp_dump_array);                                    ## Do this for the last element because it has not been tested for the if statement.
	$pos_start  = $tmp_dump_array[0];
	$GapsPosLength{$pos_start} = $gap_length;
	
	return \%GapsPosLength;
}

#################################################################################################################
#### Given an event, find whether it was an insertion or deletion. Also get the length of this indel		 ####
#### ## Just find which hash this unique key belongs to, then get its value from the respective hash		 ####
#################################################################################################################

sub GetNatureLength			
{							
	my($eventPos,$Ref2InsHash,$Ref2DelHash) = @_;
	
	my %InsHash = %$Ref2InsHash;
	my %DelHash = %$Ref2DelHash;
	
	my $nature = my $length = "";
	
	if (exists $InsHash{$eventPos})
	{
		$nature = "Insertion";
		$length = $InsHash{$eventPos};
	}
	else
	{
		$nature = "Deletion";
		$length = $DelHash{$eventPos};
	}
		
	return($nature,$length);
}

#################################################################################################################
#### Get position in the reference for a query species														 ####
#################################################################################################################

sub GetReferencePos
{
	my($eventPos,$nature,$Ref2RefGenomeHash,$Ref2GeneArray,$length,$Ref2RefGenomeCdsHash,$checkFile) = @_;
	
	my %RefGenomeHash = %$Ref2RefGenomeHash;
	my @GeneArray = @$Ref2GeneArray;
	my %RefGenomeCdsHash = %$Ref2RefGenomeCdsHash;
	
	my $refPos = my $refGenomeCds = my $refGenomeParam = "";
	
	my $index = "";
	my $flagModify = "";
	
	if ( ($nature eq "Deletion") || ($nature eq "MissingData") )
	{
		$index = $eventPos;
	}
	
	### Now there is some tricky business for Insertion
	
	if ($nature eq "Insertion")		### Here I need to discriminate if the insertion is at the end of an exon, or at the beginning of the following exon.
	{
		my $l = $GeneArray[$eventPos];
		
		my @tmp = split(/\t/,$l);
		my $cdsInsertQuery = $tmp[10];	### Get the genomic coordinate of the insertion from the query species
		my $chrInsertQuery = $tmp[11];	### Get the chromosome/scaffold from where this genomic coordinate is coming

		#### Now check the position of this genomic coordinate in the query species tabular file.
		
		my $lNumber = 1;
		my $lNumberPrint = "";

		open(FIC,$checkFile);
		while (my $line = <FIC>)
		{
			$line =~s/\s+$//;
			my @tmpLine = split(/\t/,$line);
	
			if ($tmpLine[10] =~/\d/)
			{
				if ( ($tmpLine[10] == $cdsInsertQuery) && ($tmpLine[11] eq $chrInsertQuery) )	## Also check this comes from the same scaffold. This is important because for short scaffolds, you might encounter the same coordinate
				{										## and as a result pick up thr wrong line.
					$lNumberPrint = $lNumber;
					last;	
				}
			}
			$lNumber++;
		}
		close FIC;		## Get the line number where the insertion genomic coordinate occurs in the tabular file
	
		$lNumberPrint--;
		my $lineCheck = `sed -n $lNumberPrint,$lNumberPrint\\p $checkFile`;	## Now extract a line which is 1nt upstream of this insertion.
		my @tmpLC = split(/\t/,$lineCheck);

		if ($tmpLC[5] =~/_A_/)		### If this line is annotated as an Acceptor coordinate, e.g. 5_A_1, them set the flagModify flag to "T"
		{							### The insertion indeed occurs at the beginning of the exon	
			$flagModify = "T";
			$index = $eventPos + $length;
		}
		else
		{
			$index = $eventPos - 1;
		}
	}
	
	$refPos = $RefGenomeHash{$index};
	$refGenomeCds = $RefGenomeCdsHash{$index};
	$refGenomeParam = $GeneArray[$index];
	
	if ($flagModify eq "T")
	{
		my @tmpPos = split(/_/,$refPos);
		$refPos = "$tmpPos[0]\_0";
		
		my @tmpRefParam = split(/\t/,$refGenomeParam);
		my $refStrand = $tmpRefParam[4];
		
		if ($refStrand eq "+")
		{
			$refGenomeCds--;
		}
		else
		{
			$refGenomeCds++;
		}
	}
	
	return ($refPos,$refGenomeCds,$refGenomeParam);
}

################################################################################################################################################################
#### The purpose of this function is as follows: suppose there is a deletion that is stored as following in the DeletionsHash (5 => 10,  30 =>3)			####
#### Turns out that the first deletion is spread over 2 exons, exon 1 is simply 6 nts long. So this deletion should be split into two parts					####
#### something like 5=>2, 7 =>8  (deletion of 2nts at 5th position followed by deletion of 8nts at 7th position)											####
################################################################################################################################################################

sub checkDeletionSpan	
{						
	my ($Ref2DelList,$Ref2GeneArray) = @_;		
	
	my @GeneArray = @$Ref2GeneArray;
	my %DelList   = %$Ref2DelList;
	
	my @AllDel = keys(%DelList);
	my @AllDelSorted = sort {$a <=> $b} @AllDel;
	
	my %DelLengthHash = ();
	my $flag = "F";	
	
	foreach my $event(@AllDelSorted)	## Loop over all deletions, they have to be sorted first
	{
		my $length = $DelList{$event};
	
		my $startLine = $GeneArray[$event];		## Beginning of deletion	
		my $stopLine  = $GeneArray[$event+$length - 1];		## End of deletion
		
		my @Start = split(/\t/,$startLine);			## The exon number at the start of deletion will be obtained from this line
		my @End   = split(/\t/,$stopLine);			## The exon number at the end of deletion will be obtained from this line
		
		my $NumberExons = 0;
		
		if ($Start[5] ne $End[5])					## If the two are not the same, the deletion spends over more than one exon
		{
			$flag = "T";			## Set this flag to true.
			
			my $exonS =  $Start[5];		
			$exonS =~s/Coding_exon_//;		## This is the exon number where the deletion starts
			my $exonE =  $End[5];
			$exonE =~s/Coding_exon_//;		## This is the exon number where the deletion ends
			
			$NumberExons = $exonE - $exonS;	
			
			my $start = $event;
			my $stop  = $event + $length;
		
			while ($exonS <= $exonE)
			{
				my $term = "Coding_exon\_$exonS";
				my $length = 0;
			
				for (my $j = $start;$j < $stop; $j++)		## This "for loop" gets the length
				{
					my @tmp = split(/\t/,$GeneArray[$j]);
				
					if ($tmp[5] eq $term) 
					{
						$length++;
					}	
				}
				
				my $startPos = "";
				for (my $j = $start;$j < $stop; $j++)		## This "for loop" gets the position
				{
					my @tmp = split(/\t/,$GeneArray[$j]);
				
					if ($tmp[5] eq $term) 
					{
						$startPos = $j;
						last;
					}	
				}
				
				$DelLengthHash{$startPos} = $length;		## Once we have both, we dump them in the DelLengthHash hash.
				
				$exonS++;
			}
		}
		
		else			## Otherwise it is much simpler
		{
			$DelLengthHash{$event} = $length;	
		}
	}
	
	return \%DelLengthHash;					
}

################################################################################################
#### Check events which disrupt the reading frame/ whose length is not a multiple of 3.		####
################################################################################################

sub checkFDEvents			
{
	my $Ref2Events = shift;
	my %EventsList = %$Ref2Events;
	
	my %FDEvents = ();
	my $FDcount = 0;
	
	foreach my $keys( keys(%EventsList))
	{
		my $length = $EventsList{$keys};
		
		if ($length%3 != 0)
		{
			$FDEvents{$keys} = $length;
			$FDcount++;
		}
	}
	
	return(\%FDEvents,$FDcount);
}

############################################################################################################
#### This function returns all the events which can compensate for each other. A real pain in the a$$	####
############################################################################################################

sub GetCompensatedEvents
{
	my ($Ref2IndelsAllSorted,$Ref2GeneArray,$Ref2FDEvents,$Ref2GenomeHash,$Ref2RefGenomeCds,$checkFile,$ref2ExcludedInsertions,$ref2ExcludedDeletions) = @_;		### 
	
	my @IndelsAllSortedFD = @$Ref2IndelsAllSorted;	## an array that contains all the indels (their positions), the array is sorted
	my @GeneArray 		  = @$Ref2GeneArray;		## array that contains all the useful lines from the tabular file
	my %FDEvents 		  = %$Ref2FDEvents;			## a hash that contains all the FD events, the key is the position, the value is the length of the indel ( +ve for insertions, -ve for deletions).
	my %excludedInsertions	  = %$ref2ExcludedInsertions;
	my %excludedDeletions	  = %$ref2ExcludedDeletions;

	my %CompensationHash = ();		### This is the hash that contains all the events which can be compensated, the key is the first event, the value is the event (or events) which can compensate the first event.
	my %ExcludeEvents = ();			### Any event which can compensate for another event is dumped into this hash, so that it is not considered again.

	for (my $i = 0;$i < @IndelsAllSortedFD;$i++)
	{
		my $queryEvent = $IndelsAllSortedFD[$i];	## this is the query, now lets find events which can compensate for this query event.
		
		if (! exists $ExcludeEvents{$queryEvent})
		{
			my $length = $FDEvents{$queryEvent};	## get the length of the query event
			my $natureQuery = "Insertion";
			$natureQuery = "Deletion" if ($length < 0);
			
			my $sumLengthsPCE = $length;
			my @potCompEventList = ();
			my $pce = 0;
			
			for (my $j = $i + 1;$j < @IndelsAllSortedFD;$j++)	## Now start looking for compensatory events. Remember if the index of the query event is i, the compensatory events should be search from i+1
			{
				my $potCompEvent = $IndelsAllSortedFD[$j];
				
				if (! exists $ExcludeEvents{$potCompEvent})		## the potential compensatory event should not have been "consumed" in compensating for some other event.
				{
					my $lengthPCE = $FDEvents{$potCompEvent};
				
					$sumLengthsPCE = $sumLengthsPCE + $lengthPCE;
					$potCompEventList[$pce] = $potCompEvent;
					$pce++; 
					
					if ($sumLengthsPCE%3 == 0)		## if the length of query is +2, length of potential compensating event is -5, $sumLengthsPCE = -3, so the events are compensatory.
					{

						my $start = $queryEvent;
						my $stop  = $potCompEvent + abs($lengthPCE) - 1;		## Go until the end of the event. For instance if the compensatory event is an Insertion, go until the end of the insertion and check if there is a stop codon in the species reading frame.
						
						my $naturePCE = "Insertion";
						$naturePCE = "Deletion" if ($lengthPCE < 0);
						
						my $tmpString = join(",",@potCompEventList);
						
						my $startUpstream = $queryEvent - 1;
						
						my $speciesRF = 0;
						
						$speciesRF = GetAncestralReadingFrame($startUpstream,\@GeneArray);	## What the function returns is the ancestral reading frame
								
						## Now I make it the species reading frame here:
						
						if ($speciesRF == 1)		{	$speciesRF = 2;	}
						elsif ($speciesRF == 2) 	{	$speciesRF = 3;	}
						elsif ($speciesRF == 3)		{	$speciesRF = 1;	}
						
						### Ancestral reading frame is zero when we are trying to find the ancestral reading frame for a deletion that occurs right at the beginning of the gene. In this case, startUpstream = -1
							
						my $stopCodonFlag = "";
						$speciesRF = 1 if ($speciesRF == 0);

						$stopCodonFlag = CheckForStopCodons($start,$stop,$Ref2GeneArray,$ref2ExcludedInsertions,$ref2ExcludedDeletions,$speciesRF);
						
						if ($stopCodonFlag eq "T") ## if there is no stop codon
						{
							foreach my $CE(@potCompEventList)	
							{
								$ExcludeEvents{$CE} = "T";	## then put all the compensatory events in the %ExcludeEvents hash.
							}
							
							my $compEventsString = join(",",@potCompEventList);
							$CompensationHash{$queryEvent} = $compEventsString;	## Also populate the CompensationHash accordingly.
						}
						
						last;		## and leave this loop, since we have already found compensatory event(s) for our query events.
					}
				}		## otherwise continue to look for compensatory events
			}			## End of the FOR loop that looks for compensatory events
		}				## IF statement which checks that the queryEvent has not already been "consumed"/ it is a compensatory event for some other upstream event.
	}					## End of the FOR loop that iterates through the list of all Indels.
	
	
	### Perform the "Return To Ancestral Reading Frame Test"
	
	foreach my $events( sort {$a <=> $b} keys(%CompensationHash))
	{	
		## Get eventStop position, or the position of the last of the compensatory events. Remember, compensatory events could be a list of events
		my @tmpCompEvents = split(",",$CompensationHash{$events});
		my $stop = $tmpCompEvents[$#tmpCompEvents];
		
		my $startUpstream = $events - 1;
		my $speciesRF = GetAncestralReadingFrame($startUpstream,\@GeneArray);	## What the function returns is the ancestral reading frame. since the species reading frame
																				## has been ancestralized, speciesRF = ancestralRF;
		my $length = abs($FDEvents{$stop});
		
		## Only deletions at the beginning of the gene sequence are excluded from the "Return To Ancestral Reading Frame Test"
		if ($speciesRF ne "")
		{		
			my $start = $events;		## start position is where the first event happens
			my $stop = $stop + $length -1;	## stop position is where the last of the compensatory events happen + the length of the event - 1
			
			my $readingFrameRF = "";
			
			while ($start <= $stop)
			{
				if (! exists $excludedInsertions{$start})
				{
				
					my $line = $GeneArray[$start];
					my @tmp = split(/\t/,$line);
				
					### Reading frame for the reference species
				
					$readingFrameRF = $tmp[6] if ($tmp[6] ne "-");	## Keep getting the reading frame of the Reference species for every line where there is sequence in the Reference
				
					### Reading frame for the query species
					$speciesRF++  if ( ($tmp[13] ne "-") || (exists $excludedDeletions{$start}) );				## Same as for Reference. Even if there is a ?, we assume that there is sequence here. That is why, the if condition is 
																## << if ($tmp[13] ne "-") >>  NOT  <<< if ( ($tmp[13] ne "-") AND if ($tmp[13] ne "?") ) >>>
				}
				
				$start++;
			}
			
			my $speciesRFPrint = $speciesRF%3;
			$speciesRFPrint = 3 if ( $speciesRFPrint == 0);
			$speciesRFPrint = "C".$speciesRFPrint;
			
			## If the queryEvent and its compensatory event are spread over exons, we believe that this compensation is very likely not real. So we do not report these events.
			
			die "This compensatory event pair ...$events .. $CompensationHash{$events} does not return the species reading frame to the ancestral reading frame\n" if ( ($readingFrameRF ne $speciesRFPrint) );
		}
	}
	
	return \%CompensationHash;
}

#############################################################################################################
#### This function checks for the presence of stop codons in a region that is apparently compensated for ####
#############################################################################################################

sub CheckForStopCodons					
{										
	my ($start,$stop,$Ref2Gene,$ref2ExcludedInsertions,$ref2ExcludedDeletions,$readingFrame) = @_;
	
	my @GeneArray = @$Ref2Gene;
	my %excludedInsertions = %$ref2ExcludedInsertions;
	my %excludedDeletions  = %$ref2ExcludedDeletions;
	
	my $FlagStop = "";
	my $geneSeq = "";
	
	for(my $k = $start;$k < $stop; $k++)
	{
		if (! exists $excludedInsertions{$k})
		{
			my @tmpLine = split(/\t/,$GeneArray[$k]);
		
			$geneSeq = $geneSeq.$tmpLine[13];
			$geneSeq = $geneSeq."N" if (exists $excludedDeletions{$k});
		}
	}
	
	if  ($geneSeq eq "")  ### This happens when two exon deletions can compensate for each other. For example, an exon deletion of 22nt (equivalent to -1) followed by an exon deletion of 23nt nt(equivalent to -2).
					### In such cases there will be no readingFrame and no geneSeq. And we have to believe that the events are indeed compensatory
	{
		$FlagStop = "T";
	}
	
	else
	{				
		$geneSeq =~s/\?/N/g;
		$geneSeq =~s/-//g;
	
		if ($readingFrame == 2)		{	$geneSeq = "N".$geneSeq;	}
		elsif ($readingFrame == 3)	{	$geneSeq = "NN".$geneSeq;	}	
	
		my $stopCodon = find_stop_codon($geneSeq,"first");
		$FlagStop = "T" if ($stopCodon eq "");
	}
	
	return $FlagStop;
}

####################################################################
#### Get the ancestral reading frame for a codon position 	    ####
####################################################################

sub GetAncestralReadingFrame
{
	my ($startUpstream,$ref2GeneArray) = @_;
	
	my @GeneArray = @$ref2GeneArray;
	my $speciesRF = "";
	
	while ($startUpstream >= 0)		## This is the index in the GeneArray from where I ancestralize the reading frame of the species i.e. a positon just upstream of the indel
	{
		my $line = $GeneArray[$startUpstream];
		my @tmp = split(/\t/,$line);
		
		if  ($tmp[7] ne "-")		## Get the line where there is a character in the Reference species. Remember we need the position where there is a character in the Reference, not in the query species
		{
			$speciesRF = $tmp[6];
			$speciesRF =~s/C//;
			
			last;
		}
		$startUpstream--;
	}
	
	$speciesRF = 0 if ($speciesRF eq "");
	
	return $speciesRF;
}

sub checkFlankingNts		## This function is a new addition to this version of IndelDetectionModule
{							## Based on quality issues.
	
	my ($Ref2Hash,$Ref2GeneArray,$type) = @_;
	
	my %eventsHash = %$Ref2Hash;		## A hash that contains the positions of events (Insertions|Deletions) as keys, and their lengths as values
	my @GeneArray  = @$Ref2GeneArray;	## An array that contains all the lines from the tabular file
	my %HashFiltered = ();				## This hash is what is eventually returned, filtered hash, all the "dubious" indels are filtered out
	my %HashExclude = ();
	
	foreach my $keys(sort {$a <=> $b} keys(%eventsHash))	## Loop through these events, in a sorted manner
	{
		my $length = $eventsHash{$keys};

		my $start = $keys - 1;		
		my $stop  = $keys + $length;
		
		my $ntUpStream = my $ntDownStream = "";
		
		while ($start >= 0)		## Get the nucleotide upstream of this event (insertion|deletion)
		{
			my $line = $GeneArray[$start];
			my @tmp = split(/\t/,$line);
			
			if ( ($tmp[13] ne "-") && ($tmp[13] ne "?") )
			{
				$ntUpStream = $tmp[13];
				last;
			}
			$start--;
		}
		
		while ($stop < scalar(@GeneArray)) ## Get the nucleotide downstream of this event (insertion|deletion)
		{
			my $line = $GeneArray[$stop];
			my @tmp = split(/\t/,$line);
			
			if ( ($tmp[13] ne "-") && ($tmp[13] ne "?") )
			{
				$ntDownStream = $tmp[13];
				last;
			}
			$stop++;
		}
	
		my $seqInsert = "";

		if ($type eq "Insertion")		## Get the sequence of insert from the tabular file.
		{
			my $start = $keys;
			my $stop  = $start + $length - 1;

			while ($start <= $stop)
			{
				my $line = $GeneArray[$start];
				my @tmp = split(/\t/,$line);
				$seqInsert = $seqInsert.$tmp[13];
				$start++;
			}
		}

		## Criteria for NOT reporting an event: based on observations and discussions with Michael (17/06/2014)
		### 1) if either of the flanking nucleotides is "N" (for both insertions and deletions)
		### 2) if the insert has an "N", even a single N.	
		### The logic is: we might see a 13nt insertion/deletion (note that it is frame disrupting), flanked by an "N" or has an "N" (in case of insertion).
		### Now if this "N" is not real, then the indel becomes frame preserving. Hence we are better off NOT reporting such cases.
		
		if ( ($ntUpStream eq "N") || ($ntDownStream eq "N") || ($seqInsert =~/N/) )
		{
			$HashExclude{$keys} = $length;
		}
		else
		{
			$HashFiltered{$keys} = $length;
		}
	}
	
	return (\%HashFiltered,\%HashExclude);
}

