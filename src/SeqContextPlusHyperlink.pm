#!/sw/bin/perl

## Virag Sharma, Final version June 2014
## This module contains all the functionalities related to printing the URL to the genome browser, plus the sequence context around an event (Indel|InframeStopCodon|SpliceSiteMutation)
## Since this is a bit harder than I initially thought, so I decided to write an entire module for this functionality

package SeqContextPlusHyperlink;
use strict;
use warnings;
use Exporter;

our @ISA = ('Exporter');
our @EXPORT = qw(Link_SeqContext); 

sub Link_SeqContext		## The main function that accepts the arguments and returns the values
{
	my ($file,$refCds,$strandRef,$type) = @_;
	
	die "Event type - $type not defined for the function Link_SeqContext \n" if ( ($type ne "IFSC") && ($type !~ /Insertion/) && ($type !~/Deletion/)  && ($type !~/SSM_/) );

	my $seqContextStart = my $seqContextStop = "";
	my $length = 0;		## Length is the length of the event that you want to show
	
	if ($type eq "IFSC"){ ## Length of inframe stop codon is always 3 nt
		$length = 3;
	}elsif ($type =~/SSM/){	## Length of a mutated splice site is always 2nt
		my $junction = (split /_/,$type)[1];
		$length = 2;
		
		if ($junction eq "Acceptor"){
			$refCds = $refCds - 1 if ($strandRef eq "+"); 
			$refCds = $refCds + 1 if ($strandRef eq "-");
		}
	}elsif ($type =~/Insertion/ || $type =~/Deletion/)		## This is variable and depends on the particular event
	{
		$length = (split /_/,$type)[1];		## That is why the parameter $type (or more specifically eventType) is passed differently to this function
								## $type[0] has the value Insertion|Deletion ... $tmp[1] is the length of the indel
	}
	
	$refCds-- if ($type =~/Insertion/ && $strandRef eq "-");		## Has been tested, works
	$refCds++ if ($type =~/Insertion/ && $strandRef eq "+");
	
	## What we need is 6nts on both sides of the event.
	## Hence set the values of $seqContextStart and $seqContextStop accordingly. We work in the reference genome coordinates here
	
	if ($type eq "IFSC" || $type =~/SSM/ || $type =~/Deletion/){
		if ($strandRef eq "+"){
			$seqContextStart = $refCds - 6;
			$seqContextStop  = $seqContextStart + 11 + $length;
		}else{
			$seqContextStart = $refCds - 5 - $length;
			$seqContextStop  = $seqContextStart + 11 + $length;
		}
	}else{
		if ($strandRef eq "+"){
			$seqContextStart = $refCds - 6;
			$seqContextStop  = $seqContextStart + 12;
		}else{
			$seqContextStart = $refCds - 5;
			$seqContextStop  = $seqContextStart + 11;
		}
	}
		
	($seqContextStart,$seqContextStop) = fixContextStartStop($file,$seqContextStart,$seqContextStop,$refCds); ### In this case, the seqContextStart and/or seqContextStop are not
	my @list = (3,3,$length,3,3); ## This array contains the length of the elements in which the sequence will be split. Read below to see the details
	
	## Get line numbers in the tabularFile array where seqContextStart and seqContextStop occur.
	my %refGenomeCdsHash = ();
	
	my $l = 1;
	open(FI,$file) || die "Error opening tabulatr file '$file' in Link_SeqContext function\n";
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		my $cds = (split /\t/,$line)[2];
		$refGenomeCdsHash{$cds} = $l - 1 if ($cds ne "-" && $cds ne "?");
		
		$l++;
	}
	
	my $lineStart = $refGenomeCdsHash{$seqContextStart};
	my $lineStop = $refGenomeCdsHash{$seqContextStop};			### line numbers for seqContextStart and seqContextStop done.
	### Do swapping if lineStart is greater than lineStop. Case with the minus strand
	
	if ($lineStart > $lineStop)
	{
		my $tmpSwap = $lineStart;
		$lineStart = $lineStop;
		$lineStop = $tmpSwap;
	}
	
	open(FI,"$file") || die "Error opening file $file in Link_SeqContext function\n";
	my @fileArray = <FI>;
	close FI;
	
	my $seqRef = my $seqQuery = "";
	my @selectedLinesRef = my @selectedLinesQuery = ();
	my $l_ref = my $l_query = 0;

	while ($lineStart <= $lineStop)
	{
		my $line = $fileArray[$lineStart];
		my($ntRef,$ntQuery) = (split /\t/,$line)[7,13];
	
		push(@selectedLinesRef,$line);
		push(@selectedLinesQuery,$line);
	
		$seqRef = $seqRef.$ntRef;
		$seqQuery = $seqQuery.$ntQuery;
	
		$lineStart++;
	}
	
	my ($LinkRef,$speciesRef)   = ProcessLines(\@selectedLinesRef,"Ref");		## The most important function. 
	my ($LinkQuery,$speciesQuery,$flagCheck) = ProcessLines(\@selectedLinesQuery,"Query");
	$LinkQuery = "-" if ($flagCheck eq "F");
	### The hyperlinks have been obtained, now I need to get a formatted sequence that I can print in the box.
	
	my $formattedSeqRef = my $formattedSeqQuery = "";
	
	if ($type eq "IFSC" || $type =~/SSM/){		### Get the formatted sequence for IFSC/SSM events.
		$formattedSeqRef   = FormatSeq($seqRef,\@list);
		$formattedSeqQuery = FormatSeq($seqQuery,\@list);
	}else{											## This is a bit hard, since indels can vary in length. So the formatting of the sequence here
													## is handled differently by a separate function
		my @tmp = split(/_/,$type);
		$type = $tmp[0];
		
		($formattedSeqRef,$formattedSeqQuery) = FormatSeqIndel($seqRef,$seqQuery,\@list,$type);		
	}
	
	### Just put additional spaces at the end of one of the species' name (either Reference or Query) when the two names are not of equal length,so that the two strings become equal in length.
	if (length($speciesRef) ne length($speciesQuery))
	{
		my $l1 = length($speciesRef);
		my $l2 = length($speciesQuery);
		
		my $diff = abs($l1 - $l2);
		my $spaceString = "";
		
		my $i = 0;
		while ($i < $diff)
		{
			$spaceString = $spaceString." ";
			$i++;
		}
		
		if ($l1 > $l2){	
			$speciesQuery = $speciesQuery.$spaceString;
		}else{	
			$speciesRef   = $speciesRef.$spaceString;	
		}
	}
	
	my $returnString = $LinkQuery."<><cc>(".$speciesQuery.")<cc>".$formattedSeqQuery."<>".$LinkRef."<><cc>(".$speciesRef.")<cc>".$formattedSeqRef;
	return $returnString;
}

#### End of the main body of the code. Sub-routines follow.
#### The following are the functions that are called by the Link_SeqContext function

sub fixContextStartStop		## Get the nearest values of startConextCds and stopContextCds which have a sequence/base in the Reference species, not "-"
{
	my($file,$startContext,$stopContext,$refCds) = @_;
	my %hashCds = ();
	
	open(FI,"$file") || die "Error opening input file '$file' in fixContextStartStop in SeqContextPlusHyperLink module \n";
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		my $cds = (split /\t/,$line)[2];
		$hashCds{$cds} = "T" if ($cds ne "-");
	}
	close FI;
	
	my $returnStart = my $returnStop = "";
	while ($startContext <= $refCds)
	{
		if (exists $hashCds{$startContext})
		{
			$returnStart = $startContext;
			last;
		}
		$startContext++;
	}
	
	while ($stopContext >= $refCds)
	{
		if (exists $hashCds{$stopContext})
		{
			$returnStop = $stopContext;
			last;
		}
		$stopContext--;
	}
	
	return($returnStart,$returnStop);
}

sub processFlanks
{
	my ($length,$direction) = @_;
	
	return(3,3) if ($length >= 6);
	
	if ($length > 3)
	{
		my $rem = $length%3;
		if ($direction eq "L"){
			return($rem,3);
		}else{
			return(3,$rem);
		}
	}
	else{
		if ($direction eq "L"){
			return(0,$length);
		}else{
			return($length,0);
		}	
	}
}



sub ProcessLines				### This function performs the most important part of this entire module. The input is the reference to an array which contains the lines starting from the beginning
{						### of the sequence context to the end of it. Another input is the label - Ref for Reference
	my ($ref2List,$label) = @_;
	my @list = @$ref2List;
	
	my $lineFirst = $list[0]; chomp $lineFirst;
	my $lineLast  = $list[$#list]; chomp $lineLast;
	
	####print "Let me print the lines as well\n... $lineFirst\n .. $lineLast\n";	### Test statement, now commented. But nevertheless, keep it.
	my @tmp1 = split(/\t/,$lineFirst);			## Get the first line and the last line.
	my @tmp2 = split(/\t/,$lineLast);
	my $flagEmpty = "";
	
	my $speciesIndex = 1;
	my $cdsIndex = 2;
	my $chrIndex = 3;
	my $strandIndex = 4;
	
	if ($label ne "Ref"){
		$speciesIndex = 9;
		$cdsIndex = 10;
		$chrIndex = 11;
		$strandIndex = 12;
	}
	
	### Some processing to extract the chromosome, the strand and the cdsStart/cdsStop
	my $chr        = $tmp1[$chrIndex];
	my $strand     = $tmp1[$strandIndex];
	my $cds1       = $tmp1[$cdsIndex];
	my $cds2       = $tmp2[$cdsIndex];
	my $species    = $tmp1[$speciesIndex];
	
	if ($cds1 ne "?" && $cds2 ne "?" && $cds1 ne "-" && $cds2 ne "-"){
		if ($cds1 > $cds2){  			### Swap cds1 and cds2 if they come from minus strand
			my $tmp = $cds2;
			$cds2 = $cds1;
			$cds1 = $tmp;
		}
		
		if ($strand eq "+"){
			$cds1++;
			$cds2++;
		}
	}else{
		$flagEmpty = "F";
	}
	
	my $link = "http://genome.pks.mpg.de/cgi-bin/hgTracks?org=$species&db=$species&position=$chr:$cds1-$cds2";
	
	if ($label eq "Ref"){	
		return ($link,$species);
	}else{	
		return($link,$species,$flagEmpty);	
	}		## Return an additional variable in case the Query species is passed to this function
}


sub FormatSeq		### This function formats the sequence for InframeStopCodons/Splice site mutations which is relatively straightforward
{					### The input is the sequence and an array which has the list of different fragment sizes, 
	my($seq,$Ref2List) = @_;	### so the array is either (3 3 2 3 3) or (3 3 3 3 3)
	
	my @list = @$Ref2List;
	
	my $i = my $start = 0;
	my @seqSplit = ();
	
	foreach my $unit(@list)
	{
		$seqSplit[$i] = substr($seq,$start,$unit);
		$i++;
		$start = $start + $unit;
	}	
	$seqSplit[2] = "<cc>".$seqSplit[2]."<cc>";	### Just put these tags around the event, so that it is highlighted in the final SVG
		
	my $SeqReturn = join(" ",@seqSplit);
	return $SeqReturn;
}

sub FormatSeqIndel		### Same as FormatSeq described above. However, here we need to do some trimming depending on the size of the Insertion|Deletion.
{
	my($seqRef,$seqQuery,$Ref2List,$type) = @_;
	
	my @list = @$Ref2List;
	my $start = 0;
	my @seqSplitRef = my @seqSplitQuery = ();
			
	foreach my $unit(@list)
	{
		push(@seqSplitRef,substr($seqRef,$start,$unit));
		push(@seqSplitQuery,substr($seqQuery,$start,$unit));
		$start = $start + $unit;
	}
	
	my $length = $list[2];
	### For obvious reasons, I cannot print very long inserts/deletions. So need to trim them
	
	if ($length > 10)
	{
		if ($type eq "Insertion")
		{
			$seqSplitRef[2]   = "------------------";		## The element that is shown in the Reference is simply a string of dashes
			$seqSplitQuery[2] = " Insert ($length    nts)";	## The element that is shown in the Query is "Insert (10 nts)"
		}
		else
		{
			$seqSplitQuery[2] = "Deletion ($length   nts)";		
			
		        my $l1 = substr($seqSplitRef[2],0,2);			## The element that is shown in the Query is "Deletion (11 nts)"
			my $l2 = substr($seqSplitRef[2],length($seqSplitRef[2]) - 2,2);
			$seqSplitRef[2]  = "$l1..............$l2";		## If the deletion is longer than 10 nts, simply show the context, 3 nts on either side and .... in between
		}
	}
	
	$seqSplitQuery[2] = "<cc>".$seqSplitQuery[2]."<cc>";	## Put tags around the events
	$seqSplitRef[2]   = "<cc>".$seqSplitRef[2]."<cc>";
		
	my $seqReturnRef   = join(" ",@seqSplitRef);
	my $seqReturnQuery = join(" ",@seqSplitQuery);
	
	return($seqReturnRef,$seqReturnQuery);		## Return the formatted sequence for the query and the Reference
}
