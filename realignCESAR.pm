package realignCESAR;
use strict;
use warnings;

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;
use MyKentFunctions;
use MyBDB;

use Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(writeMafToBDB fixCoordinates getSequenceFrom2Bit fixPhase getFlankingSeq getCoordinatesFromMaf getRefSeqForMaf 
createPairWiseMaf createPairWiseMaf_ExonGroups createMultiMaf getExonStartStop createString flipStopCodons flipBack getAlignmentCdsHash 
qualityMaskMaf getChromSize createOutDir replaceNAtEnds checkMafCorrectness writeOnlyReferenceSequence flipBackN safeDie swapValues);

sub writeMafToBDB		## This function writes a maf to a BDB file
{
	my ($maf,$key,$bdb,$ref2TempList) = @_;
	
	open(FIM,$maf) || safeDie("Error !! Cannot open maf file '$maf' in the writeMafToBDB file\n",$ref2TempList);
	my @mafArray = <FIM>;
	close FIM;
	`rm -rf $maf`;
	
	my $alignmentString = join("",@mafArray);
	## the key is "$reference#$transcriptIDs#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag"
	my @tmp = split("#",$key);
	
	my $tList = splice(@tmp,1,1);
	my $ref = shift(@tmp);
	my $remKey = join("#",@tmp);
	
	my @transcriptsList = split(",",$tList);
	
	foreach my $t(@transcriptsList)
	{
		my $keyNew = "$ref#$t#$remKey";
		secureWriteBDB($bdb,$keyNew,$alignmentString,1,1);	
	}
}

sub fixCoordinates		## This function fixes coordinates for non-middle exons, basically mask start/stop codons.
{
	my ($start,$stop,$strand,$tag) = @_;
	
	if ($tag eq "F"){
		if ($strand eq "+"){			## Suppose the coordinates are 10-20. What I would want to extract is 13-20
			$start = $start + 3;
		}else{							## Example: 80-100. Extract: 80-97
			$stop = $stop - 3;
		}
	}elsif ($tag eq "L"){
		if ($strand eq "+"){		## Example: 100-120: Extract: 100-117
			$stop = $stop - 3;
		}else{					## Example: 50-60. Extract 53-60
			$start = $start + 3;
		}
	}elsif ($tag eq "S"){			## Example 100-200: Extract: 103-197 irrespective of the strand
		$start = $start + 3;
		$stop = $stop - 3;
	}
	
	return ($start,$stop);
}

sub getSequenceFrom2Bit  ## This function gets sequence from2  bit file given the twoBit file, the chrom, the start and the stop coordinate
{
	my ($twoBitFile,$chr,$start,$stop,$ref2TempList) = @_;
	my $seq = `twoBitToFa $twoBitFile -seq=$chr -start=$start -end=$stop stdout|grep -v ">"|tr -d "\n"`;
	safeDie("Cannot get the sequence from 2BitToFa",$ref2TempList) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);	
	return $seq;
}

sub fixPhase 	## This function does the following: If the input sequence is GTGCATCGAAT, phase5Prime is 1, phase3Prime is 1, the function returns gTGCATCGAAt
{
	my ($seq,$phase5Prime,$phase3Prime,$program) = @_;  
	$seq = uc($seq);	## Make entire sequence upper-case before you do anything to it
	
	if ($phase5Prime != 0){ ## Fix 5' end
		my $s1 = substr($seq,0,$phase5Prime);
		my $s2 = substr($seq,$phase5Prime);
		
		$seq = lc($s1).$s2;
		#$seq = lc($s1)."|".$s2 if ($program eq "CESAR2");
	}
	
	if ($phase3Prime != 0){ ## Fix 3' end
		my $s1 = substr($seq,0,(length($seq) - $phase3Prime));
		my $s2 = substr($seq,(length($seq) - $phase3Prime));
		
		$seq = $s1.lc($s2);
		#$seq = $s1."|".lc($s2) if ($program eq "CESAR2");
	}
	return $seq;
}

sub getFlankingSeq	### This function gets flanking sequence upstream and downstream for each exon for the reference species:
{
	my ($exonStart,$exonStop,$chr,$strand,$twoBitFileRef,$tag,$refFlank,$ref2TempFilesList) = @_;
	
	my $startPos = $exonStart - $refFlank;
	my $stopPos = $exonStop + $refFlank;
	
	if ($tag eq "F"){
		if ($strand eq "+"){			## Suppose the coordinates are 10-20. What I would want to extract is 13-20
			$startPos = $startPos - 3;
		}else{						## Example: 80-100. Extract: 80-97
			$stopPos = $stopPos + 3;
		}
	}elsif ($tag eq "L"){
		if ($strand eq "+"){		## Example: 100-120: Extract: 100-117
			$stopPos = $stopPos + 3;
		}else{					## Example: 50-60. Extract 53-60
			$startPos = $startPos - 3;
		}
	}elsif ($tag eq "S"){			## Example 100-200: Extract: 103-197 irrespective of the strand
		$startPos = $startPos - 3;
		$stopPos  = $stopPos + 3;
	}
	
	my $seqUpstream = getSequenceFrom2Bit($twoBitFileRef,$chr,$startPos,$exonStart,$ref2TempFilesList);
	my $seqDownstream = getSequenceFrom2Bit($twoBitFileRef,$chr,$exonStop,$stopPos,$ref2TempFilesList);
	
	if ($strand eq "-"){
		$seqUpstream 	= revComp($seqUpstream);
		$seqDownstream	= revComp($seqDownstream);
	}
	
	$seqUpstream 	= lc($seqUpstream);
	$seqDownstream	= lc($seqDownstream);
	
	my $flankingSeqString = "";
	$flankingSeqString = "$seqUpstream#$seqDownstream" if ($strand eq "+");
	$flankingSeqString = "$seqDownstream#$seqUpstream" if ($strand eq "-");
	
	return $flankingSeqString;
}

sub getRefSeqForMaf		## This function gives us the first part of the s-Line for the reference sequence: For example: s hg19.chr17 55822154 496 +  81195210
{
	my ($refSeq,$refCds,$flankingSeqString,$strand,$tag,$ref2TempList) = @_;
	
	## Occassionally, the refSeq could have gaps
	my $refSeqCopy = $refSeq;
	$refSeqCopy =~s/-//g;
	$refSeqCopy =~s/\|//g;
	my $lengthRef = length($refSeqCopy);
	
	my ($acceptorSeq,$donorSeq) = (split /#/,$flankingSeqString)[0,1];
	my ($refSpecies,$refStart,$refStop, $refChr,$refStrand) = (split /#/,$refCds)[0,2,3,4,5];
	my $sizeChrRef = getChromSize($refSpecies,$refChr,$ref2TempList);
	
	## Get the incoming and the outgoing phases for the exon
	my ($phase5Prime,$phase3Prime) = getExonStartStop($refSeq);		## This function does not really give the phase. Instead it gives the position of the first lower case character after the CDS has terminated
	$phase3Prime = length($refSeq) - $phase3Prime if ($phase3Prime != 0);		## So use this to get the 3' phase of the exon
	
	## Now deal with the flanking sequence up and down of the reference that I insert in the maf: they DO NOT come from the alignment
	$acceptorSeq =~s/\s+$//;
	$donorSeq 	 =~s/\s+$//;
	
	my $refPosition = $refStart;
	$refPosition = $refPosition - length($acceptorSeq);
	$refPosition = $refPosition + 3 if ($tag eq "F" || $tag eq "S");
	
	$lengthRef = $lengthRef + length($acceptorSeq) + length($donorSeq);
	
	### Format the s-Line for the Reference depending on the strand
	## For reference, we print the coordinates as such. When we run mafExtract, it extracts from these coordinates. It is my downstream processing step that will flip the coordinates
	## and revComp the sequence if the gene is on minus strand
	
	my $sLineRef_Part1 = "s $refSpecies.$refChr $refPosition $lengthRef + $sizeChrRef"; ## The last but one field i.e. the strand is always "+" for the reference
	return $sLineRef_Part1;	
}

sub createPairWiseMaf  ## This function creates a pairwise maf from the fasta file: Basically the heart of all the processing. 
{
	my ($reference,$refSeqAligned,$refSeq_1,$querySeqAligned,$queryLine,$flankingSeqString,$refStrand,$phase5Prime,$phase3Prime,$tag,$strand,$seqRefOri,$list2CodonsFlipped,$species,$allPairWiseAlignment,$ref2TempList) = @_;

	my $refSeqAlignedOri   = $refSeqAligned;
	my $querySeqAlignedOri = $querySeqAligned;
	$querySeqAligned =~s/#/-/g;   ## Replace all hashes with dashes, here we are talking about deletions in the CDS
	
	## First operation --> If the reference sequence contained any in-frame stop codon and it was flipped to NNN, put the original stop codon back
	my @listOfInFrameStops = @$list2CodonsFlipped;
	
	if (scalar(@listOfInFrameStops) > 0)	## Do this only if the reference sequence does contain inframe stop codons.
	{
		my $ref2Hash_SeqCds2AlignmentCds = getAlignmentCdsHash($list2CodonsFlipped,$refSeqAligned,$phase5Prime);
		my $refSeqAlignedNew = flipBack($refSeqAligned,$seqRefOri,$ref2Hash_SeqCds2AlignmentCds,$phase5Prime,$phase3Prime);
		$refSeqAligned = $refSeqAlignedNew;
	}
	
	my ($acceptorSeq,$donorSeq) = (split /#/,$flankingSeqString)[0,1];
	$acceptorSeq = uc($acceptorSeq);
	$donorSeq    = uc($donorSeq);
	
	$refSeqAligned = $acceptorSeq.$refSeqAligned.$donorSeq;
	$refSeqAligned = revComp($refSeqAligned) if ($refStrand eq "-"); ## revComp the reference sequence if the gene is on minus strand
	
	$refSeqAligned =~s/\s+//g; ## and remove the empty space
	$refSeqAligned = uc($refSeqAligned);
	my $sLineRef = "$refSeq_1 $refSeqAligned";	## Thats all for the reference:
		
	my ($alignmentStart,$alignmentStop) = getExonStartStop($querySeqAligned); ## Tells us the positons where the alignment starts and ends
	my $querySeqAlignedR = reverse($querySeqAligned);
	my ($alignmentStartR,$alignmentStopR) = getExonStartStop($querySeqAlignedR); ## Tells us the positons where the alignment starts and ends
	
	my $lQSA = length($querySeqAligned);
	
	($phase5Prime,$phase3Prime) =  getExonPhases_Query($refSeqAlignedOri,$querySeqAlignedOri);
	my $querySeqAligned_CDS = uc(getQueryAlignedSeq($refSeqAlignedOri,$querySeqAligned));
	
	#print "The queryAlignedCDS is '$querySeqAligned_CDS'\n";
	#print "The phase is $phase5Prime ==> 5P,$phase3Prime ==> 3P\n";
	
	#my $querySeqAligned_CDS = uc(substr($querySeqAligned,$alignmentStart-$phase5Prime,($alignmentStop-$alignmentStart+$phase5Prime+$phase3Prime)));
	
	#($phase5Prime,$phase3Prime) =  getExonPhases_Query($refSeqAlignedOri,$querySeqAlignedOri);
	### This is important: Imagine the sequence is like this: atcgatagaGTGGTGagtccc
	## The alignment start is 9 (first Capital letter), the phase5Prime is 1, the alignmentStop is 15 and the phase3Prime is 1
	## The above gives you the right CDS which is aGTGGTGa.

	## Append 23bp gaps on either side of the queryCDS, because the max length of flanking sequence is 23
	my $stringGaps = createString(23,"-");
	$querySeqAligned = $stringGaps.$querySeqAligned.$stringGaps;
	
	my $seqUpQuery = substr($querySeqAligned,($alignmentStart-$phase5Prime - length($acceptorSeq)+23),length($acceptorSeq));
	$seqUpQuery =~s/x/-/ig;  ## Now replace all x (deletions in intron) with a dash;
	my $seqUpQueryCopy = $seqUpQuery;
	$seqUpQueryCopy =~s/-//g; 
	my $lAccQuery = length($seqUpQueryCopy); ## These are the number of bases that are available to us, the coordinates should be stated with this length in mind
	
	my $seqDownQuery = substr($querySeqAligned,($alignmentStop+$phase3Prime+23),length($donorSeq));
	$seqDownQuery =~s/x/-/ig;
	my $seqDownQueryCopy = $seqDownQuery;
	$seqDownQueryCopy =~s/-//g; ## Now replace all x (deletions in intron) with a dash;
	my $lDonQuery = length($seqDownQueryCopy);  ## These are the number of bases that are available to us, the coordinates should be stated with this length in mind
	
	$querySeqAligned_CDS =~s/x/-/ig;
	$querySeqAligned_CDS = uc($seqUpQuery).$querySeqAligned_CDS.uc($seqDownQuery);
	
	## Now get the length of the query
	my $querySeqAligned_CDSUngapped = $querySeqAligned_CDS;
	$querySeqAligned_CDSUngapped =~s/-//g;
	my $lengthQuery = length($querySeqAligned_CDSUngapped);
	
	my ($querySpecies,$queryStart,$queryStop,$queryChr,$queryStrand) = (split /#/,$queryLine)[0,1,2,3,4];
	my $sizeChrQuery = getChromSize($querySpecies,$queryChr);
	my $sLineQuery = "";
	
	## Now modify the alignment start and alignmentStop here, the starts include the positions of Ns..we need to exclude those because those are not real bases
	my $seqStart = substr($querySeqAlignedOri,0,$alignmentStart);
	$seqStart =~s/x//g;
	my $alignmentStartN = length($seqStart);
	$alignmentStop = $alignmentStop - ($alignmentStart - $alignmentStartN); ## Move it by the same distance;
	$alignmentStart = $alignmentStartN;
	
	my $seqStartR = substr($querySeqAlignedR,0,$alignmentStartR);
	$seqStartR =~s/x//g;
	my $alignmentStartRN = length($seqStartR);
	$alignmentStopR = $alignmentStop - ($alignmentStartR - $alignmentStartRN); ## Move it by the same distance;
	$alignmentStartR = $alignmentStartRN;
	
	### At this stage, I need to fix the coordinates. Note that only the coordinates require flipping: the sequence is fine.
	if ($refStrand ne $queryStrand){
		$alignmentStart = $alignmentStartR;
		$alignmentStop  = $alignmentStopR;

		my $tmp = $phase5Prime;
	 	$phase5Prime = $phase3Prime;
		$phase3Prime = $tmp;	
	}
	
	## Now the final processing
	
	if ($refStrand eq "+" && $queryStrand eq "+"){
			my $queryPosition = $queryStart + ($alignmentStart - $phase5Prime) - $lAccQuery;
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "+" && $queryStrand eq "-"){
			my $queryPosition = $sizeChrQuery - $queryStart - ($alignmentStart - $phase5Prime) - $lengthQuery + $lDonQuery;
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "-" && $queryStrand eq "-"){
			my $queryPosition = $sizeChrQuery - $queryStart - ($alignmentStart - $phase5Prime) - $lengthQuery + $lAccQuery;	
			$querySeqAligned_CDS = revComp($querySeqAligned_CDS);
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "-" && $queryStrand eq "+"){
			my $queryPosition = $queryStart + ($alignmentStart - $phase5Prime) - $lDonQuery;
			$querySeqAligned_CDS = revComp($querySeqAligned_CDS);
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}
	
	open(FOTM,">>$allPairWiseAlignment") || safeDie ("Cannot write to allPairWise file '$allPairWiseAlignment' in createPairWiseMaf function\n",$ref2TempList);
	print FOTM "##maf version=1 scoring=tba.v8\n";
	print FOTM "$sLineRef\n$sLineQuery\n";
	close FOTM;
}

sub getExonStartStop	## Gives me the position of the start of an exon and stop point of an exon
{
	my $seq = shift;
	my $posStart = my $posStop = 0;
	
	my @seqArray = split("",$seq);
	
	for (my $j = 0;$j < scalar(@seqArray); $j++)
	{
		if ($seqArray[$j] =~/[A-Z]/ || $seqArray[$j] eq "-"){	## i.e. in capitals, or a gapped position would indicate the beginning of the alignment 
			$posStart = $j;
			last;
		}
	}
	
	for (my $j = $posStart;$j < scalar(@seqArray); $j++)
	{
		if ($seqArray[$j] =~/[a-z]/){ 	## i.e. the first lower case letter after the end of the exonic sequence
			$posStop = $j;
			last;
		}
	}
	## In case, there is no intronic sequence in the 3' direction of the exon, then $posStop = Length of the sequence
	$posStop = length($seq) if ($posStop == 0);
		
	return ($posStart,$posStop);
}


sub createString  ## This function creates a string given the character that makes this string and the length
{
	my ($length,$char) = @_;
	my $string = "";
	
	my $i = 1;
	while ($i <= $length)
	{
		$string = $string.$char;
		$i++;
	}
	
	return $string;
}

sub flipStopCodons  ## This function flips the in-frame stop codons to a sense-codon in the reference sequence
{
	my ($seqOri,$phase5Prime,$phase3Prime) = @_;
	
	my $s1  = substr($seqOri,0,$phase5Prime);
	my $seq = substr($seqOri,$phase5Prime,length($seqOri)-$phase5Prime-$phase3Prime) ;
	my $s2  = substr($seqOri,length($seqOri)-$phase3Prime);
	
	my @list = ();
	my $i = 0;
	
	while ($i < length($seq))
	{
		my $codon = substr($seq,$i,1).substr($seq,($i+1),1).substr($seq,($i+2),1);
		if ($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){
			
			my $codonNew = $codon;
			$codonNew =~s/T/C/; ## TGA becomes CGA, TAG becomes CAG and TAA becomes CAA
			$seq = substr($seq,0,$i).$codonNew.substr($seq,($i+3));
			push(@list,$i);
		}
		$i = $i + 3;
	}
	$seq = $s1.$seq.$s2;
	return($seq,\@list);
}

sub flipBack	## This function flips back a sequence that contains a sense-codon (originally a stop codon) back to the sequence with the stop codon at the right position.
{
	my ($alignedSeq,$seqOri,$ref2cdsHash,$phase5Prime,$phase3Prime) = @_;
	
	$alignedSeq =~s/\s//g;	## Remove all spaces from alignedSeq
	my %cdsHash = %$ref2cdsHash;
	
	foreach my $keys(sort {$a <=> $b} keys(%cdsHash))
	{
		my $posStopUnAligned = $keys;
		my $posStopAligned   = $cdsHash{$keys};
		
		my $codonOri = substr($seqOri,($posStopUnAligned),3);
		$alignedSeq = substr($alignedSeq,0,$posStopAligned).$codonOri.substr($alignedSeq,($posStopAligned+3));
	}
	return $alignedSeq;
}

sub flipBackN
{
	my($alignedSeq,$ref2cdsHash) = @_;
	my %cdsHash = %$ref2cdsHash;
	
	foreach my $keys(sort {$a <=> $b} keys(%cdsHash))
	{
		my $posStopAligned   = $cdsHash{$keys};
		$alignedSeq = substr($alignedSeq,0,$posStopAligned)."n".substr($alignedSeq,($posStopAligned+1));
	}
	return $alignedSeq;
}


sub qualityMaskMaf  ## This function masks the input maf (a wrapper around mafQualityMasking.pl)
{
	my ($mafIn,$speciesList,$ref2TempList) = @_;
	my $mafOut = `mktemp`; chomp $mafOut;
	
	my @speciesListOrdered = @$speciesList;
	my $qualityMaskedFilesList = `mktemp`; chomp $qualityMaskedFilesList;
	
	open(FOQM,">$qualityMaskedFilesList");
	foreach my $species(@speciesListOrdered)
	{
		my $twoBitFileQueryMasked = "$ENV{'genomePath'}/gbdb-HL/$species/$species.qualityStrict.2bit";
		$twoBitFileQueryMasked    = "$ENV{'genomePath'}/gbdb-HL/$species/$species.quality.2bit" if(! -e $twoBitFileQueryMasked);	## Just use the species.2bit file if no qualityStrict2bit file is available
		print FOQM "$species\t$twoBitFileQueryMasked\n";
	}
	close FOQM;
	
	my $mafQualityMaskingCall = "mafQualityMasking.pl -input $mafIn -output $mafOut -listFiles $qualityMaskedFilesList -q"; 
	my $statusMasking = system($mafQualityMaskingCall);
	safeDie("Error running mafQualityMasking\n",$ref2TempList) if ($statusMasking != 0);
	`mv $mafOut $mafIn`;
	`rm $qualityMaskedFilesList`;
}

sub getChromSize  ## This function returns the chromosome size for a given species/chromosome pair 
{
	my ($species,$chrInterest,$ref2TempList) = @_;
	
	my $chromSizeFile = "$ENV{'genomePath'}/gbdb-HL/$species/chrom.sizes";
	my $chrSize = "";
	
	open(FICS,"$chromSizeFile") || safeDie("Error opening chromosomes size file '$chromSizeFile' in getChromSize function\n",$ref2TempList);
	while(my $line = <FICS>)
	{
		$line =~s/\s+$//;
		my($chr,$size) = (split /\t/,$line);
		if ($chr eq $chrInterest){
			$chrSize = $size;
			last;
		}
	}
	close FICS;
	return $chrSize;
}

sub createOutDir  ## This function creates an output directory for a given output file
{
	my $file = shift;
	my $dirName = `dirname $file`;
	chomp $dirName;	
	`mkdir -p $dirName` if ($dirName ne ".");
}

sub replaceNAtEnds  ## This function replaces N at the extremities of the sequence with "-"
{
	my ($seqIn,$length,$dir) = @_;
	my $j = 0;
	
	my @tmp = split(//,$seqIn);
	@tmp = reverse(@tmp) if ($dir eq "3P");
	for(my $i = 0; $i < scalar(@tmp); $i++)
	{
		if ($tmp[$i] eq "n"){
			$tmp[$i] = "x";
			$j++;
		}
		if ($tmp[$i] eq "N"){
			$tmp[$i] = "#";
			$j++;
		}
		last if ($j == $length);
	}
	@tmp = reverse(@tmp) if ($dir eq "3P");
	my $seqOut = join("",@tmp);
	return $seqOut;
}

sub checkMafCorrectness  ## This function runs the faTo2Bit test on an input maf
{
	my ($outMaf,$ref2TempList) = @_;
	open(FI_MAF,$outMaf) || safeDie("Error opening outMaf in checkMafCorrectness while performing the twoBitToFa test\n",$ref2TempList);
	while (my $mafLine = <FI_MAF>)
	{
		if ($mafLine =~/^s/){
			$mafLine =~s/\s+$//;
			my ($species,$chr,$start,$stop,$strand,$srcSize,$seq) = getSLineData($mafLine);
			
			#print "$start .. $stop .. $srcSize\n";
			my $begin = my $end = "";
			if ($strand eq "+"){
				$begin = $start;
				$end   = $begin + $stop;
			}else{
				$begin = $srcSize - $start - $stop;
				$end   = $begin + $stop;
			}
			
			my $twoBitFile = "$ENV{'genomePath'}/gbdb-HL/$species/$species.2bit";
			#$twoBitFile = "orcOrc1.quality.2bit" if ($species eq "orcOrc1");  ## Delete this later
			my $seqFrom2Bit = uc(getSequenceFrom2Bit($twoBitFile,$chr,$begin,$end));
			$seqFrom2Bit = revComp($seqFrom2Bit) if ($strand eq "-");
			$seq = uc($seq);
			$seq =~s/-//g;
			
			safeDie("This exon from '$species' fails the 2BitToFa sequence versus the mafSequence test\n2Bit sequence is '$seqFrom2Bit'\nMaf sequence is  '$seq'\nBegin is '$begin', end is '$end'\nCheck $outMaf for $mafLine",$ref2TempList) if ($seqFrom2Bit ne $seq && $begin != $end);
		}
	}
	close FI_MAF;
}

sub createMultiMaf  ### This function combines different pairwise mafs and makes a multimaf out of it by running maf-join
{
	my($reference,$allPairWiseMaf,$ref2ELinesHash,$ref2SpeciesList,$outBDB,$key,$masking,$ref2TempList) = @_;
	
	my %eLinesHash  = %$ref2ELinesHash;
	my @tempFilesList = @$ref2TempList;
	
	my @speciesListOrdered = @$ref2SpeciesList;
	
	## Combine the pairwise maf to multimaf:
	my $tmpMaf = `mktemp /dev/shm/XXXXX.file`; chomp $tmpMaf;
	push(@tempFilesList,$tmpMaf);	
	my $mafJoin = "maf-joinMod.py $allPairWiseMaf > $tmpMaf";
	system($mafJoin) == 0 || safeDie("Error running maf-join, the call '$mafJoin' failed\n",\@tempFilesList);
	
	## create a sLines hash, contains all the sLines from the tmpMaf produced in the above step
	## This stuff is only for aesthetic purpose. I do this only to make sure that the order of the species in this maf is the same as in the original maf
	my %sLinesHash = ();
	
	my $refSLine = "";
	open(FM_IN,"$tmpMaf") || safeDie("Error opening the nearly finished '$tmpMaf' in createMultiMaf function\n",\@tempFilesList);
	while (my $line = <FM_IN>)
	{
		$line =~s/\s+$//;
		if ($line ne "a" && $line ne "") {
			my ($species,$chr,$d2,$d3,$d4,$d5,$d6) = getSLineData($line);
			$sLinesHash{$species} = $line;
			$refSLine = $line if ($species eq $reference);
		}
	}
	close FM_IN;
	`rm $tmpMaf`;	## Delete this maf. All the useful lines have been dumped into the hash
	
	## Final printing:
	
	my $outMaf = `mktemp /dev/shm/XXXXX.file`; chomp $outMaf;
	open(FM_OUT,">$outMaf") || die "Cannot write to the final maf file '$outMaf' in createMultiMaf function\n";
    	print FM_OUT "##maf version=1\n";
	print FM_OUT "a score=0.000000\n";    
	print FM_OUT "$refSLine\n";
	
	foreach my $species(@speciesListOrdered)
	{
		if (exists $sLinesHash{$species}) {
			print FM_OUT "$sLinesHash{$species}\n";
		} elsif (exists $eLinesHash{$species})  {
			print FM_OUT "$eLinesHash{$species}\n";
		}
	}
	print FM_OUT "\n";
	close FM_OUT;

	## Now check if everything is OK with the maf, if not die at this point.
	my $mafSpeciesListCall = "mafSpeciesList $outMaf stdout";
	my $speciesList_MafCall = `$mafSpeciesListCall`;
	safeDie ("All is NOT well with the final maf produced for the above exon. mafSpeciesList failed\n",\@tempFilesList) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);

	checkMafCorrectness($outMaf,\@tempFilesList); ## Run twoBitToFa test right here:
	qualityMaskMaf($outMaf,$ref2SpeciesList,\@tempFilesList) if ($masking); ## qualityMask the maf here
	writeMafToBDB($outMaf,$key,$outBDB); ## Write this alignment to BDB file:
}

sub writeOnlyReferenceSequence
{
	my ($ref,$flankingSeqString,$seqRefPrint,$refSeq_1,$strand,$refID,$outBDB,$seqRefOri,$list2CodonsFlipped,$phase5Prime,$phase3Prime,$masking,$ref2ELinesHash,$ref2SpeciesList,$ref2TempList) = @_;
	
	my %eLinesHash = %$ref2ELinesHash;
	my @speciesListOrdered = @$ref2SpeciesList;
	my @tempFilesList = @$ref2TempList;
	
	## First operation --> If the reference sequence contained any in-frame stop codon and it was flipped to NNN, put the original stop codon back
	my @listOfInFrameStops = @$list2CodonsFlipped;
	if (scalar(@listOfInFrameStops) > 0)	## Do this only if the reference sequence does not contain inframe stop codons.
	{
		my $refSeqAligned = $seqRefPrint; ## No alignment in this case so the unaligned sequence is the same as the aligned sequence, pretty much
		my $ref2Hash_SeqCds2AlignmentCds = getAlignmentCdsHash($list2CodonsFlipped,$refSeqAligned,$phase5Prime);
		my $refSeqAlignedNew = flipBack($refSeqAligned,$seqRefOri,$ref2Hash_SeqCds2AlignmentCds,$phase5Prime,$phase3Prime);
		$seqRefPrint = $refSeqAlignedNew;
	}
	
	my ($acceptorSeq,$donorSeq) = (split /#/,$flankingSeqString)[0,1];
	my $refSeqAligned = uc($acceptorSeq).uc($seqRefPrint).uc($donorSeq);
	$refSeqAligned = revComp($refSeqAligned) if ($strand eq "-");
	
	my $outMaf = `mktemp /dev/shm/XXXXX.file`; chomp $outMaf;
	push(@tempFilesList,$outMaf);
	
	open(FM_OUT,">$outMaf") || safeDie("Cannot write to the final maf file in writeOnlyReferenceSequence function\n",$ref2TempList);
	print FM_OUT "##maf version=1\n";
	print FM_OUT "a score=0.000000\n";    
	print FM_OUT "$refSeq_1 $refSeqAligned\n";
	
	foreach my $species(@speciesListOrdered)
	{
		print FM_OUT "$eLinesHash{$species}\n" if (exists $eLinesHash{$species});
	}
	print FM_OUT "\n";
	close FM_OUT;
	
	my @speciesList = ();
	push(@speciesList,$ref);

	qualityMaskMaf($outMaf,\@speciesList,\@tempFilesList) if ($masking);
	writeMafToBDB($outMaf,$refID,$outBDB); ## Write this alignment to BDB file:
}


sub getCoordinatesFromMaf ## This function gives me the coordinates (the start/stop/chr/strand) for the different species in the maf
{
	my ($maf,$reference,$ref2TempFilesList) = @_;

	## Get the list of species present in maf
	my @speciesList = ();			### create a species_array which contains the list of all the species present in the maf block.
	@speciesList = `mafSpeciesList $maf stdout`; 
	chomp(@speciesList); 
	die "ERROR: mafSpeciesList must have failed because the list of species is empty\n" if (scalar @speciesList < 1); 
	$" = " ";	# separate array elements by a space in the output and remove the reference from the list
	@speciesList = grep{$_ ne $reference} @speciesList;

	my %cdsStartHash = my %cdsStopHash = ();
	my %chrQuery = my %strandQuery = ();
	my %chrCountHash = my %strandCountHash = ();

	my %cdsListHash = my %eLineListHash = ();
	foreach my $species(@speciesList)
	{
		my @cdsList = my @eLineList = ();
		$cdsListHash{$species} 		= \@cdsList;
		$eLineListHash{$species}	= \@eLineList;
		$strandQuery{$species}		= "NA";
		$chrQuery{$species}		= "NA";
	}
	
	my %speciesPresent = ();
	open(FIM,"$maf") || safeDie ("Error !! Cannot open mafFile '$maf' in getCoordinatesFromMaf function\n",$ref2TempFilesList); 
	while (my $line = <FIM>)
	{
		$line =~s/\s+$//;
		next if ($line =~/$reference/);
		if ($line =~/^s/){	## For s lines
			my ($species,$chr,$start,$stop,$strand,$srcSize,$seq) = getSLineData($line);
			$speciesPresent{$species} = "T";
			
			if ($species ne $reference){
				$chrQuery{$species} = $chr;
				$strandQuery{$species} = $strand;
				
				my $cdsStartQuery = my $cdsStopQuery = "";
																	## An example is
																	## ##maf version=1 scoring=mafExtract
				if ($strand eq "-"){									## a score=0.000000
																	## s hg19-galGal4.chr19  3261974 109 + 61431566 TTCAGGCAGATGGTGGCTGAGGCAGTAGCGGTGGCTACAGTGCATGCAGAACTGGCCCAGAGTGGTGGTGCTGGCCGAGCACTTGGAGAAGCTACAGGTGTTGTCAGCC	
					$cdsStartQuery = $srcSize - $start - $stop;		## s galGal4.chr19      49577681 109 - 51267647 TTCAGGCAGATGGTGGCTGAGACAGTAGCGGTGGCTACAGTGCATGCAGAACTGCCCTAGGGTGGTGGTACTGGCTGAACACTTGGAGAAGCTACAGGTCTTGTCAGCC
					$cdsStopQuery  = $srcSize - $start;				##	
				}else{													## a score=0.000000
					$cdsStartQuery = $start;												## s hg19-galGal4.chr19  3262083 61 + 61431566 ---TTCACCACAGCAGACACCAGGGCGTCGAAGTCCTCCTCACAGGGCAGAGCCATAACTGGAC
					$cdsStopQuery  = $start + $stop;												## s galGal4.chr19      49577790 64 - 51267647 TTCTTCACCACAGCAGACACCAGGGCGTCAAAATCCTCCTCACAGGGCAGTGCCATAACTGGAC
				}													##
																	## startFirstBlock  = 51267647 - 49577681 - 109		==>	1689857		
																	## stopFirstBlock   = 51267647 - 49577681			==> 	1689966	-- We want this
																	##
																	## startLastBlock  = 51267647 - 49577790 -64		==> 	1689793	-- We want this
																	## stopLastBlock   = 51267647 - 49577790			==> 	1689857
				my $string = "$cdsStartQuery-$cdsStopQuery";
				
				my $ref2CdsList = $cdsListHash{$species};  ## get the relevant cdsList
				my @cdsList = @$ref2CdsList;
				push(@cdsList,$string);  ## push the string to the cdsList
				$cdsListHash{$species} = \@cdsList;  ## and update the value of the cdsList in the cdsListHash
				
				$chrCountHash{$species}{$strand}++;
				$strandCountHash{$species}{$chr}++;
			}
		}elsif ($line =~/^e/){
			my($species,$dm1,$dm2,$dm3,$dm4,$dm5,$dm6) = getELineData($line);
			my $ref2eLineList = $eLineListHash{$species};
			my @eLineList = @$ref2eLineList;
			
			push(@eLineList,$line);
			$eLineListHash{$species} = \@eLineList;
			$speciesPresent{$species} = "T";
		}
	}
	close FIM;
	
## + strand is straight-forward as always. Minus is a pain, as always !!!
##maf version=1 scoring=mafExtract
# a score=0.000000
# s hg19.chr11 118532101 14 + 135006516 TGTAGAAAGGCGGT
# s mm9.chr9    79584924 14 - 124076172 TGAAGGAAGGCGAC
# i mm9.chr9   N 0 C 0

# a score=2010932.000000
# s hg19.chr11 118532115 53 + 135006516 GTCATTGGTGTGAGTCAAGTAGCAATCCATCATGAGGGTCAAGAGTGGGGGCT
# s mm9.chr9    79584938 53 - 124076172 ATCCTTGGTATGAGCTACATATCGATCCATCATGAGAGTCAGGAGTGGGGGCT
# i mm9.chr9   C 0 C 0

# a score=0.000000
# s hg19.chr11 118532168 50 + 135006516 GGCTCCGCTGCAGGTAGTACACGCGCCCACCATTGGGGACATGCCCATAG
# s mm9.chr9    79584991 50 - 124076172 GGCTCCGTTGCAGGTAATATATGCGTCCACCGTTGGGGATATGTCCGTAG
# i mm9.chr9   C 0 N 0
	
## With the code above, the firstElement (see below) is '44491234-44491248' and the last element is  '44491131-44491181'. The cdsStart is 44491234 and stop is 44491181. 
## This needs to be flipped. So I reverse the array for minus strand
	
	my %cdsStartQuery = my %cdsStopQuery = ();
	foreach my $species(@speciesList)
	{
		my $ref2CdsList = $cdsListHash{$species};
		my @cdsList = @$ref2CdsList;
		
		if (scalar(@cdsList) > 0){
			@cdsList = reverse(@cdsList) if ($strandQuery{$species} eq "-");

			my $firstElement = $cdsList[0];
			my $lastElement  = $cdsList[$#cdsList];
		
			$cdsStartQuery{$species} = (split /-/,$firstElement)[0];
			$cdsStopQuery{$species}  = (split /-/,$lastElement)[1];
		}else{
			$cdsStartQuery{$species} = "NA";
			$cdsStopQuery{$species} = "NA";
		}
	}
	
	my %moreThanOneStrand = my %moreThanOneChr = my %eLinePrint = ();
	foreach my $species(@speciesList)
	{
		$moreThanOneStrand{$species} 	 = "F";
		$moreThanOneChr{$species} 	 = "F";
		$eLinePrint{$species}		 = "NA";
	}
	## check which species come from more than one strand, species which more come from than one chromosome/scaffold and species which have nothing but eLines
	foreach my $species(@speciesList)
	{
		my @chrList 	= keys(%{$chrCountHash{$species}});
		my @strandList  = keys(%{$strandCountHash{$species}});
		
		$moreThanOneStrand{$species} = "T"  if (scalar(@chrList) > 1);
		$moreThanOneChr{$species} = "T"     if (scalar(@strandList) > 1);
		
		my $queryStart = $cdsStartQuery{$species};
		my $queryStop  = $cdsStopQuery{$species};
		
		if ($queryStart eq "NA" && $queryStop eq "NA"){
			my $ref2ELineList = $eLineListHash{$species};
			my @eLineList = @$ref2ELineList;
			$eLinePrint{$species} = $eLineList[0];
		}
	}
	
	return(\%cdsStartQuery,\%cdsStopQuery,\%strandQuery,\%chrQuery,\%moreThanOneStrand,\%moreThanOneChr,\%eLinePrint,\%speciesPresent);
}

sub safeDie
{
	my ($message,$ref2TempFilesList) = @_;
	
	my @tempFilesList = @$ref2TempFilesList;
	my $tempFilesString = join(" ",@tempFilesList);
	`rm -rf $tempFilesString`;
	die $message;
}

sub swapValues
{
	my($v1,$v2) = @_;
	my $tmp = $v1;
	$v1 = $v2;
	$v2 = $tmp;
	return($v1,$v2);
}

sub getAlignmentCdsHash			## This function converts sequence coordinates into alignment coordinates
{					## Takes two arguments -the sequence coordinate (ungapped) list and the aligned sequence
	my($ref2CdsList,$alignedSeq,$phase5Prime) = @_;

	$alignedSeq =~s/\s//g;	## Remove all spaces from alignedSeq:
	my @tmp_array = split(//,$alignedSeq);
	my @allCds = @$ref2CdsList;
	
	my %hash = ();
	
	foreach my $element(@allCds)
	{
		my $cdsUnAligned = my $cdsAligned = 0;
		$element  = $element + $phase5Prime;
	
		foreach my $base(@tmp_array)
		{
			if ($base ne "-"){
				if ($cdsUnAligned == $element){
					$hash{$element} = $cdsAligned;
					last;
				}
				$cdsUnAligned++;
			}
			$cdsAligned++;	
		}
	}
	
	return \%hash;
}

sub createPairWiseMaf_ExonGroups  ## This function creates a pairwise maf from the fasta file: Basically the heart of all the processing. 
{
	my ($reference,$refSeqAligned,$refSeq_1,$querySeqAligned,$queryLine,$flankingSeqString,$refStrand,$phase5Prime,$phase3Prime,$tag,$strand,$seqRefOri,$list2CodonsFlipped,$species,$ref2TempList) = @_;
	
	my $refSeqAlignedOri   = $refSeqAligned;
	my $querySeqAlignedOri = $querySeqAligned;
	$querySeqAligned =~s/#/-/g;   ## Replace all hashes with dashes, here we are talking about deletions in the CDS

	## First operation --> If the reference sequence contained any in-frame stop codon and it was flipped to NNN, put the original stop codon back
	my @listOfInFrameStops = @$list2CodonsFlipped;
	
	if (scalar(@listOfInFrameStops) > 0)	## Do this only if the reference sequence does contain inframe stop codons.
	{
		my $ref2Hash_SeqCds2AlignmentCds = getAlignmentCdsHash($list2CodonsFlipped,$refSeqAligned,$phase5Prime);
		my $refSeqAlignedNew = flipBack($refSeqAligned,$seqRefOri,$ref2Hash_SeqCds2AlignmentCds,$phase5Prime,$phase3Prime);
		$refSeqAligned = $refSeqAlignedNew;
	}
	
	my ($acceptorSeq,$donorSeq) = (split /#/,$flankingSeqString)[0,1];
	$acceptorSeq = uc($acceptorSeq);
	$donorSeq    = uc($donorSeq);
	
	$refSeqAligned = $acceptorSeq.$refSeqAligned.$donorSeq;
	$refSeqAligned = revComp($refSeqAligned) if ($refStrand eq "-"); ## revComp the reference sequence if the gene is on minus strand
	
	$refSeqAligned =~s/\s+//g; ## and remove the empty space
	$refSeqAligned = uc($refSeqAligned);
	my $sLineRef = "$refSeq_1 $refSeqAligned";	## Thats all for the reference:
		
	my ($alignmentStart,$alignmentStop) = getExonStartStop($querySeqAligned); ## Tells us the positons where the alignment starts and ends
	my $querySeqAlignedR = reverse($querySeqAligned);
	my ($alignmentStartR,$alignmentStopR) = getExonStartStop($querySeqAlignedR); ## Tells us the positons where the alignment starts and ends

	## Recompute phase5Prime and phase3Prime by looking at the alignment because the phase5Prime in the reference and the query may not be the same all the time
	## Below is the example:
#>referenceExon
#                             gATCCACCTGATGGCTGGCCGAGTACCCCAGGGAGCTGATCGAATAGCAGTCAAGGCTGAGATGGAGGCCGTTTTTCTGGAGAACCTGAGGCATGCAGCTGGGGTTTTGGCTCAGGAGGACCTCGTGGGACTGCTGGAGCCCATCAACACCCGCATCACTGACCCCCAGTACTTCCTGGACACGCCCCAGCAGg        
#>sarHar1#1945223#1945453#chr1_GL834466_random#-
#cgtacgcacgggagctcggctgcccccag-GTGCATCTGATGGCAGGTCGCGTCCCGCAAGGTGCGGAGCGCTCCGCCGTGGCCCGGGACATGGAGACCGTGTTCGTGGAGAACCTGCGGCACGCCGCCGACGTGCTGGCCCGCGAGAACCTCGTGGGCCTCGTGGAACCCATCAACAGCCGGCTCACGGACCCGCGCTACTTCCTGGACACGCCTGAGCAGgctgcctct	

## In this case, the phase5Prime is 1 in ref but 0 in query. phase3Prime = 1 in both reference and query 	
	
	($phase5Prime,$phase3Prime) =  getExonPhases_Query($refSeqAlignedOri,$querySeqAlignedOri);
	#print "Alignment start is '$alignmentStart', stop is '$alignmentStop'\n";
	#print "Phase5P is '$phase5Prime', phase3 is '$phase3Prime'\n";
	#my $querySeqAligned_CDS = uc(substr($querySeqAligned,$alignmentStart-$phase5Prime,($alignmentStop-$alignmentStart+$phase5Prime+$phase3Prime)));
	my $querySeqAligned_CDS = uc(getQueryAlignedSeq($refSeqAlignedOri,$querySeqAligned));
	#print "This is the sequence\n'$querySeqAligned_CDS'\n"; exit;
	
	### This is important: Imagine the sequence is like this: atcgatagaGTGGTGagtccc
	## The alignment start is 9 (first Capital letter), the phase5Prime is 1, the alignmentStop is 15 and the phase3Prime is 1
	## The above gives you the right CDS which is aGTGGTGa.

	## Append 23bp gaps on either side of the queryCDS, because the max length of flanking sequence is 23
	my $stringGaps = createString(23,"-");
	$querySeqAligned = $stringGaps.$querySeqAligned.$stringGaps;
	
	my $seqUpQuery = substr($querySeqAligned,($alignmentStart-$phase5Prime - length($acceptorSeq)+23),length($acceptorSeq));
	$seqUpQuery =~s/x/-/ig;  ## Now replace all x (deletions in intron) with a dash;
	my $seqUpQueryCopy = $seqUpQuery;
	$seqUpQueryCopy =~s/-//g; 
	my $lAccQuery = length($seqUpQueryCopy); ## These are the number of bases that are available to us, the coordinates should be stated with this length in mind
	#print "The length is '$lAccQuery', the sequence is '$seqUpQueryCopy'\n";
	
	my $seqDownQuery = substr($querySeqAligned,($alignmentStop+$phase3Prime+23),length($donorSeq));
	$seqDownQuery =~s/x/-/ig;
	my $seqDownQueryCopy = $seqDownQuery;
	$seqDownQueryCopy =~s/-//g; ## Now replace all x (deletions in intron) with a dash;
	my $lDonQuery = length($seqDownQueryCopy);  ## These are the number of bases that are available to us, the coordinates should be stated with this length in mind
	
	$querySeqAligned_CDS =~s/x/-/ig;
	#print "Comp1\t$seqUpQuery\nComp2\t$querySeqAligned_CDS\nComp3\t$seqDownQuery\n"; #exit;
	$querySeqAligned_CDS = uc($seqUpQuery).$querySeqAligned_CDS.uc($seqDownQuery);
	
	## Now get the length of the query
	my $querySeqAligned_CDSUngapped = $querySeqAligned_CDS;
	$querySeqAligned_CDSUngapped =~s/-//g;
	my $lengthQuery = length($querySeqAligned_CDSUngapped);
	
	my ($querySpecies,$queryStart,$queryStop,$queryChr,$queryStrand) = (split /#/,$queryLine)[0,1,2,3,4];
	my $sizeChrQuery = getChromSize($querySpecies,$queryChr);
	my $sLineQuery = "";
	
	## Now modify the alignment start and alignmentStop here, the starts include the positions of Ns..we need to exclude those because those are not real bases
	my $seqStart = substr($querySeqAlignedOri,0,$alignmentStart);
	$seqStart =~s/x//g;
	my $alignmentStartN = length($seqStart);
	$alignmentStop = $alignmentStop - ($alignmentStart - $alignmentStartN); ## Move it by the same distance;
	$alignmentStart = $alignmentStartN;
	
	my $seqStartR = substr($querySeqAlignedR,0,$alignmentStartR);
	$seqStartR =~s/x//g;
	my $alignmentStartRN = length($seqStartR);
	$alignmentStopR = $alignmentStop - ($alignmentStartR - $alignmentStartRN); ## Move it by the same distance;
	$alignmentStartR = $alignmentStartRN;
	
	### At this stage, I need to fix the coordinates. Note that only the coordinates require flipping: the sequence is fine.
	if ($refStrand ne $queryStrand){
		$alignmentStart = $alignmentStartR;
		$alignmentStop  = $alignmentStopR;

		my $tmp = $phase5Prime;
	 	$phase5Prime = $phase3Prime;
		$phase3Prime = $tmp;	
	}
	## Now the final processing
	
	if ($refStrand eq "+" && $queryStrand eq "+"){
			my $queryPosition = $queryStart + ($alignmentStart - $phase5Prime) - $lAccQuery;
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "+" && $queryStrand eq "-"){
			my $queryPosition = $sizeChrQuery - $queryStart - ($alignmentStart - $phase5Prime) - $lengthQuery + $lDonQuery;
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "-" && $queryStrand eq "-"){
			my $queryPosition = $sizeChrQuery - $queryStart - ($alignmentStart - $phase5Prime) - $lengthQuery + $lAccQuery;	
			$querySeqAligned_CDS = revComp($querySeqAligned_CDS);
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}elsif ($refStrand eq "-" && $queryStrand eq "+"){
			my $queryPosition = $queryStart + ($alignmentStart - $phase5Prime) - $lDonQuery;
			$querySeqAligned_CDS = revComp($querySeqAligned_CDS);
			$sLineQuery = "s $querySpecies.$queryChr $queryPosition $lengthQuery $queryStrand $sizeChrQuery $querySeqAligned_CDS";
	}
	
	
	my $tmpMaf = `mktemp /dev/shm/XXXXX.file`; chomp $tmpMaf;
	open(FOTM,">$tmpMaf") || safeDie ("Cannot write to tmpMaf file '$tmpMaf' in createPairWiseMaf function\n",$ref2TempList);
	print FOTM "##maf version=1 scoring=tba.v8\n";
	print FOTM "$sLineRef\n$sLineQuery\n";
	close FOTM;
	return $tmpMaf;
}

sub getExonPhases_Query
{
	my ($refSeqAligned,$querySeqAligned) = @_;

	my @tmpSplit = split(/[A-Z]/,$refSeqAligned);
	@tmpSplit = grep{$_ ne ""}@tmpSplit;
	my @tmpSplitNew = ();
	foreach my $element(@tmpSplit)
	{
		if ($element =~/-/) {
			$element =~s/-//g;
			push(@tmpSplitNew,$element) if ($element ne "");
		} else {
			push(@tmpSplitNew,$element);
		}
	}
	@tmpSplit = @tmpSplitNew;
	
	my $i = 1;
	my $phase5Prime = my $phase3Prime = 0;
	my $phase5PrimeNt = my $phase3PrimeNt = "";
	my $stringNoSpace = $refSeqAligned;
	$stringNoSpace =~s/\s+//g;
 
	if (scalar(@tmpSplit) == 1) {
		my $element = $tmpSplit[0];
		$element =~s/\s+//g;	
		my $pos = index($stringNoSpace,$element);
		
		if ($pos == 0) {
			$phase5Prime = length($element);
			$phase5PrimeNt = $element;	
			$phase3Prime = 0;
		} else {
			$phase5Prime = 0;
			$phase3Prime = length($element);
			$phase3PrimeNt = $element;
		}
	} elsif (scalar(@tmpSplit) == 2) {
		$phase5PrimeNt = $tmpSplit[0];
		$phase3PrimeNt = $tmpSplit[1];
		$phase5PrimeNt =~s/\s+//g;
		$phase3PrimeNt =~s/\s+//g;
		
		$phase5Prime = length($phase5PrimeNt);
		$phase3Prime = length($phase3PrimeNt);
	}
	
	#print "## This is what I get ##\n$phase5Prime\n$phase3Prime\n";
	
	## Now get the positions of elements, then get these elements in the query sequence
	my $phase5PrimeQuery = my $phase3PrimeQuery = 0;
	my $pos5Prime = 0;
	my $l = 0;
	
	if ($phase5Prime != 0) {
		$pos5Prime = index($refSeqAligned,$phase5PrimeNt);
		$l = length($phase5PrimeNt);
		
		my $seqInQuery = substr($querySeqAligned,$pos5Prime,length($phase5PrimeNt));
		$seqInQuery =~s/-//g;
		$phase5PrimeQuery = length($seqInQuery);

		if ($seqInQuery =~/x/) {
			$seqInQuery =~s/x//g;
			$phase5PrimeQuery = length($seqInQuery);
		}
	}
	
	if ($phase3Prime != 0) {
		my $pos = index($refSeqAligned,$phase3PrimeNt,$pos5Prime+$l);
		my $seqInQuery = substr($querySeqAligned,$pos,length($phase3PrimeNt));
		$seqInQuery =~s/-//g;
		$phase3PrimeQuery = length($seqInQuery);
		
		if ($seqInQuery =~/x/) {
			$seqInQuery =~s/x//g;
			$phase3PrimeQuery = length($seqInQuery);
		}
	}
	
	#print "## This is what I get in the query ##\n$phase5PrimeQuery\n$phase3PrimeQuery\n";
	
	die "phase5Prime query cannot be greater than 2, here you have '$phase5PrimeQuery'\n" if ($phase5PrimeQuery > 2);
	die "phase3Prime query cannot be greater than 2, here you have '$phase3PrimeQuery'\n" if ($phase3PrimeQuery > 2);
	return($phase5PrimeQuery,$phase3PrimeQuery);
}

sub getQueryAlignedSeq
{
	my ($refSeq,$querySeq) = @_;
	#print "'$refSeq'\n"; exit;
	my $start = my $stop = 0;
	my $i = my $j = 0;
	my @seqArray = split("",$refSeq);
	foreach my $base(@seqArray)
	{
		if ($base ne " "){
			$start = $i;
			last;
		}
		$i++;
	}
	$refSeq =~s/\s+$//;
	$stop = length($refSeq);
	my $queryAlignedSeq = substr($querySeq,$start,($stop-$start));
	return $queryAlignedSeq;
}
