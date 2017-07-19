#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

The input file is a list of exons with the coordinates, the strand, the chromosome, the incoming-outgoing exon phases and the exon tag(First|Middle|Last|SingleExon Gene)
Additionally, the tool requires a bigBedIndex from where it extracts the coordinates of the aligning sequence, and an output parameter (a bdb file) to which the realigned maf is written
The script also needs look up tables for every species that you want to align (produced by createIntronTables.pl)\n#####\n\n";

use strict;
use warnings;

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;
use MyKentFunctions;
use MyBDB;

use lib "/home/vsharma/cesarScripts";
use realignCESAR;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my $verbose = 0;
my $inter = 0;
my $in = my $bbIndex = my $reference = my $intronLengthBDBDir = my $outBDB = my $clade = my $program = my $parameters = "";
my $masking = 0;  ## If this is specified, then qualityMask the files:

## --set i3_i1_acc=0.25 i3_i1_do=0.25 cd_acc=0.15 cd_do=0.15
sub usage
{
	die "usage: $0 -in inFile -bbIndex bigBedIndexFile -ref RefereceSpecies -intronLookUpTables DirectoryContaininingIntronLookUpTables -outBDB outFile.bdb -clade human|mouse|drosophila|celegans -program (CESAR|CESAR2) -parameters [Optional: Only if you do not want to use the default transition probabilities. Works only with cCESAR]
	-inter[Optional - useful for bug fixes, keeps temporary files]. -mask (if you want to qualityMask the file)) \n$purpose";
}
GetOptions ("v|verbose"  => \$verbose, "in=s" => \$in, "bbIndex=s" => \$bbIndex, "ref=s" => \$reference, "intronLookUpTables=s" => \$intronLengthBDBDir, "outBDB=s" => \$outBDB, "clade=s" => \$clade, "program=s" => \$program, "parameters:s" => \$parameters, "inter" => \$inter) || usage();
usage() if ($in eq "");

### Check different parameters:
die "BigBedIndex file '$bbIndex' does not exist|is not specified\n" if (! -e $bbIndex || $bbIndex eq "");
die "Reference species not specified\n" if ($reference eq "");
die "outBDB file not specified|does not look like a bdb file\n" if ($outBDB eq "" || $outBDB !~/.bdb$/i);
die "No value specified for Clade\n" if ($clade eq "");
die "intronLookUpTables directory not specified|does not exist\n" if ($intronLengthBDBDir eq "" || ! -d $intronLengthBDBDir);
die "The parameter value for program '$program' is not recognized\n" if ($program ne "CESAR" && $program ne "CESAR2");

## Create output directories if not already present
createOutDir($outBDB);

my $twoBitFileRef = "$ENV{'genomePath'}/gbdb-HL/$reference/$reference.2bit";
die "Two bit file for the reference species '$reference' missing\n" if (! -e $twoBitFileRef);

my @tempFilesList = ();  ## this directory contains all files/dir produced by running mktemp
	
open(FI,$in) || die "Error opening the input file '$in'\n";
while (my $line = <FI>)
{
	$line =~s/\s+$//;	
	my($gene,$isoform,$chr,$strand,$exons,$phases,$tags,$species) = (split /\t/, $line)[0,1,2,3,4,5,6,7];
	my $strandRef = $strand; ## Just create a copy of this variable.
	
	print "--> Aligning the exon '$line'\n";
	my $exonDone = 0; ## set this count to zero.
	
	my @exonsList = split(",",$exons);
	my @phaseList = split(",",$phases);
	my @tagsList  = split(",",$tags);
	
	for(my $e = 0;$e < scalar(@exonsList);$e++)
	{
		my $exon 	= $exonsList[$e];
		my $phase 	= $phaseList[$e];
		my $tag  	= $tagsList[$e];
		
		my ($phase5Prime,$phase3Prime) = (split /-/,$phase)[0,1];
		my ($exonStart,$exonStop) = (split "-",$exon)[0,1];
		my $refID = "$reference#$isoform#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag";             ## This id will serve as the key in the BDB file
		
		if (-e $outBDB) {
			my $mafValue = readBDB($outBDB,$refID);
			$exonDone++ if ($mafValue ne "");
		}
	}
	
	if ($exonDone == scalar(@exonsList)) {
		print "Omitting key present in '$line', already present in the outFile'$outBDB'\n";
	} else {	## Do this only if this key is not present in the bdb file already. If it is, then do not realign the exon
		my @refIDList = my @flankingSeqList = my @sLineFirstPartList = my @lengthExonList = my @seqList = ();  ## These arrays contain values for every exon and are used later when splitting the composite exon alignment into individual exon alignments
		my $phase5PrimeFinal = my $phase3PrimeFinal = "";
		my $acceptorAC = my $donorAT = "";
		my $exonStartComp = my $exonStopComp = "";
		my $accSeqFinal = my $donorSeqFinal = "";
		my $firstExon = my $lastExon = "";
		
		my $lengthExon = 0;
		foreach (my $e = 0; $e < scalar(@exonsList); $e++)
		{
			my ($phase5Prime,$phase3Prime) = (split /-/,$phaseList[$e])[0,1];
			my ($exonStart,$exonStop) = split(/-/,$exonsList[$e]);
			my $tag = $tagsList[$e];
			
			my $refID = "$reference#$isoform#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag";
			push(@refIDList,$refID);
			
			my $lExon = $exonStop - $exonStart;
			$lengthExon = $lengthExon + $lExon;  ## The cumulative length of exons until the current exon
			push(@lengthExonList,$lExon);
			
			($exonStart,$exonStop) = fixCoordinates($exonStart,$exonStop,$strand,$tag) if ($tag ne "M"); ## Fix these coordinates right here for non-middle exons
			my $seqRefPrint = "";
			if ($exonStart != $exonStop) {
				$seqRefPrint = getSequenceFrom2Bit($twoBitFileRef,$chr,$exonStart,$exonStop,\@tempFilesList);
			} else {
				$seqRefPrint = "---";
			}
			
			$seqRefPrint = revComp($seqRefPrint) if ($strand eq "-");	## this is the sequence that is used for alignment, but its revComp is written to maf if the gene comes from minus strand
			print "### Reference sequence for exon '$exonsList[$e]' is '$seqRefPrint'\n###\n" if ($verbose);
			push(@seqList,$seqRefPrint);
			
			my $flankingSeqString = getFlankingSeq($exonStart,$exonStop,$chr,$strand,$twoBitFileRef,$tag,20,\@tempFilesList); ## Seq upstream and downstream of the exon	
			my ($accSeq,$donorSeq) = (split "#",$flankingSeqString)[0,1];
			
			my $refSeq_1 = getRefSeqForMaf($seqRefPrint,$refID,$flankingSeqString,$strand,$tag,\@tempFilesList);		## This gives me the first part of the s-Line for the Reference species
			$flankingSeqString = "$donorSeq#$accSeq" if ($strandRef eq "-");
			push(@flankingSeqList,$flankingSeqString);
			
			print "### Reference sequence information for $exonsList[$e] is '$refSeq_1' and the sequence is  '$seqRefPrint'\n###\n" if ($verbose);
			push(@sLineFirstPartList,$refSeq_1);
			
			## Now look at the features that change when you align more than one exon together, 
			## for this I need to check if the exon is the first or the last exon of the block.
			if ($e == 0) {
				$phase5PrimeFinal = $phase5Prime;
				$acceptorAC = "T" if ($accSeq =~/ac$/);
				$exonStartComp = $exonStart;
				$accSeqFinal = $accSeq;
				$firstExon = $exonsList[$e];
			}
			
			if ($e == scalar(@exonsList) -1) {
				$phase3PrimeFinal = $phase3Prime;
				$donorAT = "T" if ($donorSeq =~/^at/);
				$exonStopComp = $exonStop;
				$donorSeqFinal = $donorSeq;
				$lastExon = $exonsList[$e];
			}
		}
		
		if ($strandRef eq "-") {
			$exonStopComp  = (split /-/,$exonsList[0])[1];
			$exonStartComp = (split /-/,$exonsList[$#exonsList])[0];
		}
		
		my $refSeq = "";
		foreach my $element(@seqList)
		{
			$refSeq = $refSeq.$element;
		}
		
		my $flankingSeqFinal = "$accSeqFinal#$donorSeqFinal";
		$refSeq = fixPhase($refSeq,$phase5PrimeFinal,$phase3PrimeFinal,$program);  ## Fix the phases at the end
		
		## Temporary bug fix for CESAR
		my $refSeqWithN = $refSeq;
		$refSeq =~s/n/c/g;
		## Locate positions which had the 'n' character
		my @smallNPositions = ();
		my $nPos = 0;
		while ($nPos < length($refSeqWithN))
		{
			my $char = substr($refSeqWithN,$nPos,1);
			push(@smallNPositions,$nPos) if ($char eq "n");
			$nPos++;
		}
		
		my $seqRefOri = $refSeq;  ## Make a copy of the refSeq here and flip the in-frame stop codons with NNN
		my($refSeqNew,$list2CodonsFlipped) = flipStopCodons($refSeq,$phase5PrimeFinal,$phase3PrimeFinal);  ## Flip the stop codons to NNN
		$refSeq = $refSeqNew;

		## Ascertain the tag from the tags of the first and the last exon of the group
		my $tagFirst = $tagsList[0];
		my $tagLast  = $tagsList[$#tagsList];
		
		my $tagFinal = "";
		$tagFinal = "F" if ($tagFirst eq "F" && $tagLast eq "M");
		$tagFinal = "L" if ($tagFirst eq "M" && $tagLast eq "L");
		$tagFinal = "S" if ($tagFirst eq "F" && $tagLast eq "L");
		$tagFinal = "M" if ($tagFirst eq "M" && $tagLast eq "M");
			
#### Now run mafExtract, get the maf. From this maf, get the coordinates of the aligning locus for different species. From these coordinates, get the sequence plus/minus the length of the exon
#### and write it to a fasta file
		
		my $maf = `mktemp /dev/shm/XXXXX.file`; chomp $maf;
		push(@tempFilesList,$maf);
		
		my $mafExtractCall = "mafExtract $bbIndex -region=$chr:$exonStartComp-$exonStopComp stdout|mafSpeciesSubset stdin speciesList=$reference,$species species.lst=NULL $maf";
		system($mafExtractCall) == 0 || safeDie("Error running mafExtract for '$line'. The index was '$bbIndex'\n",\@tempFilesList);
		print "Mafextract run successful for '$line', the call was '$mafExtractCall'\n" if ($verbose);
		
		## Now extract the coordinates of the query species
		my ($ref2cdsStartQuery,$ref2cdsStopQuery,$ref2strandQuery,$ref2chrQuery,$ref2moreThanOneStrand,$ref2moreThanOneChr,$ref2eLine) = getCoordinatesFromMaf($maf,$reference,@tempFilesList);
		
		## Now get the values for the species of interest from the hashes that getCoordinatesFromMaf has returned
		my $cdsStartQuery     = $ref2cdsStartQuery->{$species};
		my $cdsStopQuery      = $ref2cdsStopQuery->{$species};
		my $strandQuery       = $ref2strandQuery->{$species};
		my $chrQuery          = $ref2chrQuery->{$species};
		my $moreThanOneStrand = $ref2moreThanOneStrand->{$species};
		my $moreThanOneChr    = $ref2moreThanOneChr->{$species};
		my $eLine             = $ref2eLine->{$species};
		
		next if ($cdsStartQuery > $cdsStopQuery); ## This indicates an inversion and for this case, I exclude such exons. Realignment_Exons.pl will catch this and write the message to the log file for such exons
		my $lengthExonQuery = $cdsStopQuery - $cdsStartQuery;
		
		print "### Aligning locus for the query species --> $cdsStartQuery,$cdsStopQuery,$strandQuery,$chrQuery,$moreThanOneStrand,$moreThanOneChr,$eLine\n\n" if ($verbose);
		## Find the rightBDBs
		my $fivePrimeBDB	= "$intronLengthBDBDir/$species.5Prime.BDB";
		my $threePrimeBDB	= "$intronLengthBDBDir/$species.3Prime.BDB";
		
		#ENST00000373995	chr10	52603747-52603882	rn5
		my $refIDFirst = "$gene#$isoform#$chr#$firstExon#$species";
		my $refIDLast  = "$gene#$isoform#$chr#$lastExon#$species";
		
		my $intronLength5Prime = readBDB($fivePrimeBDB,$refIDFirst);
		my $intronLength3Prime = readBDB($threePrimeBDB,$refIDLast);
		
		## Special cases where the intronLength cannot be obtained from the look-up table
		$intronLength5Prime = $lengthExon if ($tagFinal eq "F");  ## For the first exon, you can not have a length in 5' direction from the look-up table
		$intronLength3Prime = $lengthExon if ($tagFinal eq "L");  ## Same for the last exon
		$intronLength5Prime = $intronLength3Prime = $lengthExon if ($tagFinal eq "S");  ## And for a single exon gene
		$intronLength5Prime = $lengthExon if ($intronLength5Prime eq ""); ## if the exon in question is exon#2 and the first exon is deleted, we cannot get any value for intronLength here.
		$intronLength3Prime = $lengthExon if ($intronLength3Prime eq ""); ## if the exon in question is exon#12 and the last exon (13th) is deleted, we cannot get any value for intronLength here.
		
		### Now set the default flankingSequence length = 50
		my $flankDef = 50;
		$flankDef = $flankDef + ($lengthExon - $lengthExonQuery) if ($lengthExon > $lengthExonQuery);
		
		my $flank5Prime = my $flank3Prime = $flankDef; ## By default, we take the length of the exon as the length of the flanking sequence (very conservative)
		$flank5Prime = $intronLength5Prime  if ($flankDef > $intronLength5Prime);  ## If the exon is longer than the flank length, then we just use the flank length
		$flank3Prime = $intronLength3Prime  if ($flankDef > $intronLength3Prime);
		
		my $startString = my $stopString = "";
		my $nAtStart = my $nAtEnd = "";
		my $cdsQueryStart = my $cdsQueryStop = "";

		if ($strandQuery eq "+"){  ## These tell you the actual coordinates from where the sequence should be extracted
			$cdsQueryStart = $cdsStartQuery - $flank5Prime;	
			$cdsQueryStop  = $cdsStopQuery + $flank3Prime;
		}else{
			$cdsQueryStart = $cdsStartQuery - $flank3Prime;
			$cdsQueryStop  = $cdsStopQuery + $flank5Prime;
		}
		
		my $chromSize = getChromSize($species,$chrQuery,\@tempFilesList);
		my $shortEnd = "";
		
		if ($cdsQueryStart < 0){
			$startString = createString(abs($cdsQueryStart),"N");
			$nAtStart = abs($cdsQueryStart);
			$cdsQueryStart = 1;
			$shortEnd = "T";
		}
		
		if ($cdsQueryStop > $chromSize){
			my $diff = $cdsQueryStop - $chromSize;
			$stopString = createString($diff,"N");
			$cdsQueryStop = $chromSize;
			$shortEnd = "T";
			$nAtEnd = $diff;
		}
	
		## Now run 2BitToFa to extract the sequence
		my $twoBitFileQuery = "$ENV{'genomePath'}/gbdb-HL/$species/$species.2bit"; ## Do not worry about the the qualityMasked file here, we qualityMask the final maf later.
		my $seqQueryPrint   = getSequenceFrom2Bit($twoBitFileQuery,$chrQuery,$cdsQueryStart,$cdsQueryStop,\@tempFilesList);
		
		## if the sequence has Ns at the beginning or the end (assembly gaps), strip them right here and adjust the coordinates accordingly
		my $seqQueryCopy1 = $seqQueryPrint;
		my $seqQueryCopy2 = $seqQueryPrint;
		
		$seqQueryCopy1 =~s/^N+//;
		$cdsQueryStart = $cdsQueryStart + (length($seqQueryPrint)  - length($seqQueryCopy1));
		$seqQueryCopy2 =~s/N+$//;
		$cdsQueryStop = $cdsQueryStop - (length($seqQueryPrint)  - length($seqQueryCopy2));
		
		$seqQueryPrint =~s/^N+//;
		$seqQueryPrint =~s/N+$//;
		$seqQueryPrint = revComp($seqQueryPrint) if ($strand ne $strandQuery);
		
		print "Fetched fasta sequence for '$species' using the following coordinates: $species#$cdsQueryStart#$cdsQueryStop#$chrQuery#$strandQuery\n" if ($verbose);
		## Fix the strand orientation
		my $queryID = "$species#$cdsQueryStart#$cdsQueryStop#$chrQuery#$strandQuery";

		$seqQueryPrint = $startString.$seqQueryPrint.$stopString;
		$seqQueryPrint = uc($seqQueryPrint);
		
		my $shortEnd5P = my $shortEnd3P = 0;
		$shortEnd5P = $nAtStart if ($nAtStart ne "");
		$shortEnd3P = $nAtEnd if ($nAtEnd ne "");
		
		## Create a fasta file that will be used as input to CESAR
		my $fasta = `mktemp /dev/shm/XXXXX.file`; chomp $fasta;
		unshift(@tempFilesList,$fasta); ## Every temp file goes to this array
		
		open(FO,">$fasta") || safeDie("Cannot write to fasta file '$fasta'\n",\@tempFilesList);
		print FO ">referenceExon\n$refSeq\n" if ($program eq "CESAR");
		print FO ">referenceExon\n$refSeq\n####\n"  if ($program eq "CESAR2");  ## CESAR2 requires that after the end of the reference sequence, I have a string with 4 hashes.
		print FO ">$queryID\n$seqQueryPrint\n";
		close FO;
		
		print "Maf processing done and fasta file '$fasta' successfully created. Now running CESAR.py on '$fasta'\n"  if ($verbose);
		## At this point, I am done with the maf operations, all the relevant information from the maf has gone to my fasta file.
		## Now run our aligner on the fasta file with different arguments depending on the index/tag of the exon:
		my $out = `mktemp /dev/shm/XXXXX.file`; chomp $out;
		unshift(@tempFilesList,$out); ## Every temp file goes to this array
		my $call = "";
		####
		## Specify the acc/donor profiles based on the exon that you are looking at:
		my $accProfile = "extra/tables/human/acc_profile.txt";
		my $donorProfile = "extra/tables/human/do_profile.txt";
		
		if ($program eq "CESAR") {
		
			if ($tagFinal eq "F"){ 
				$call = "CESAR --is_first_exon IS_FIRST_EXON --clade $clade $fasta > $out";
				$call = "CESAR --is_first_exon IS_FIRST_EXON --has_AT_donor HAS_AT_DONOR --clade $clade $fasta > $out" if ($donorAT eq "T");
			}elsif ($tagFinal eq "L"){
				$call = "CESAR --is_last_exon IS_LAST_EXON --clade $clade $fasta > $out";
				$call = "CESAR --is_last_exon IS_LAST_EXON --has_AC_acceptor HAS_AC_ACCEPTOR --clade $clade $fasta > $out" if ($acceptorAC eq "T");
			}elsif ($tagFinal eq "S"){
				$call = "CESAR --is_first_exon IS_FIRST_EXON --is_last_exon IS_LAST_EXON --clade $clade $fasta > $out";
			}else{
				$call = "CESAR --clade $clade $fasta > $out";		
				$call = "CESAR --has_AT_donor HAS_AT_DONOR --clade $clade $fasta > $out" if ($donorAT eq "T");
				$call = "CESAR --has_AC_acceptor HAS_AC_ACCEPTOR --clade $clade $fasta > $out" if ($acceptorAC eq "T");
			}
		} else {
		
			if ($tagFinal eq "F"){ 
				$accProfile = "extra/tables/$clade/firstCodon_profile.txt";
			}elsif ($tagFinal eq "L"){
				$donorProfile = "extra/tables/$clade/lastCodon_profile.txt";
			}elsif ($tagFinal eq "S"){
				$accProfile = "extra/tables/$clade/firstCodon_profile.txt";
				$donorProfile = "extra/tables/$clade/lastCodon_profile.txt";
			}
				
			$donorProfile = "extra/tables/$clade/u12_donor_profile.txt" if ($donorAT eq "T");
			$accProfile = "extra/tables/$clade/u12_acc_profile.txt" if ($acceptorAC eq "T");
			if ($parameters eq "") {
				$call = "CESAR2 $fasta --max-memory 100GB -m extra/tables/$clade/eth_codon_sub.txt -p $accProfile $donorProfile > $out";
			} else {
				$call = "CESAR2 $fasta --max-memory 100GB -m extra/tables/$clade/eth_codon_sub.txt -p $accProfile $donorProfile --set $parameters > $out";
			}
		}
		
		my $statusAligner = system($call);
		safeDie("Error running the system call '$call'\n",\@tempFilesList) if ($statusAligner != 0);
		print "Aligment file '$out' produced. CESAR was called as '$call'\n" if ($verbose);
		
		### Now parse this alignment file and make different pairwise mafs.
		open(FIA,$out) || safeDie("Cannot open the output of CESAR --> '$out'\n",\@tempFilesList);
		my @alignmentArray = <FIA>;
		close FIA;
		
		my %pairwiseMaf = ();
		## Now iterate through this array
		for(my $k = 0; $k < scalar(@alignmentArray); $k++)
		{
			if ($alignmentArray[$k] =~/>referenceExon/ || $alignmentArray[$k] =~/>References/) {
				my $queryLine  = $alignmentArray[$k+2];
				$queryLine =~s/^>//;
				$queryLine =~s/\s+$//;
				my $species = (split /#/,$queryLine)[0];
				
				my $refSeq  	= $alignmentArray[$k+1];
				my $querySeq	= $alignmentArray[$k+3];
				$refSeq 	=~s/\s+$//;
				$querySeq	=~s/\s+$//;
				
				my $refSeqCopy = $refSeq;
				$refSeq =~s/^\s+//;
				## Flip back the small c to n
				my $ref2Hash_SeqCds2AlignmentCds = getAlignmentCdsHash(\@smallNPositions,$refSeq,0);
				$refSeq = flipBackN($refSeq,$ref2Hash_SeqCds2AlignmentCds);
				my $diff = length($refSeqCopy) - length($refSeq);
				$refSeq = createString($diff," ").$refSeq;
				
				$querySeq = replaceNAtEnds($querySeq,$shortEnd5P,"5P") if ($shortEnd5P ne "");
				$querySeq = replaceNAtEnds($querySeq,$shortEnd3P,"3P") if ($shortEnd3P ne "");
				
				$pairwiseMaf{$species} = createPairWiseMaf_ExonGroups($reference,$refSeq,"refSeqString",$querySeq,$queryLine,$flankingSeqFinal,$strand,$phase5PrimeFinal,$phase3PrimeFinal,$tagFinal,$strand,$seqRefOri,$list2CodonsFlipped,$species,\@tempFilesList);
				print "Pairwise maf $pairwiseMaf{$species} created for '$species'\n" if ($verbose);
				unshift(@tempFilesList,$pairwiseMaf{$species});
			}
		}
		
		## Time to split the composite exon alignment into individual exon alignments
		my @seqStartList = my @seqStopList = ();
		my $lExon = 20;
		
		for(my $i = 0; $i < scalar(@exonsList); $i++)
		{
			my $start = $lExon;  ## the beginning
			$start = $start + 3 if ($tagsList[$i] eq "F");  ## For the first exon
			my $stop  = $lengthExonList[$i] + $start;  ## the end which is the start point plus the length of the exon			
			$stop = $stop - 3 if ($tagsList[$i] eq "F" || $tagsList[$i] eq "L");
			
			$lExon = $lExon + $lengthExonList[$i];  ## the cumulative exon length
			push(@seqStartList,$start);
			push(@seqStopList,$stop);
			#print "Start is '$start' while stop is '$stop'\n";
		}
		
		my $alignedSeqRef = my $alignedSeqQuery = my $sLinePartQuery = "";  ## Now extract the aligned reference and the query sequence
		open(FITM,$pairwiseMaf{$species}) || die "Error opening pairWise maf for '$species'\n";
		while (my $seqLine = <FITM>)
		{
			$seqLine =~s/\s+$//;
			my @tmp = split(/\s+/,$seqLine)	;
			
			$alignedSeqRef   = $tmp[$#tmp] if ($seqLine =~/refSeqString/);
			if ($seqLine =~/^s/){
				$alignedSeqQuery = pop(@tmp);
				$sLinePartQuery = $seqLine;
			}
		}
		close FITM;

		$alignedSeqRef = reverse($alignedSeqRef) if ($strandRef eq "-");  ## This has to be reversed again.
		$alignedSeqQuery = reverse($alignedSeqQuery) if ($strandRef eq "-");
		
		my($dm1,$dm2,$startCDS,$lengthQueryTotal,$dm4,$sizeQuery,$dm5) = getSLineData($sLinePartQuery);  ## All I care about is the startCDS for the query and the chromosome size
			
		my $ref2StartList = getAlignmentCdsHash(\@seqStartList,$alignedSeqRef,0);
		my $ref2StopList  = getAlignmentCdsHash(\@seqStopList,$alignedSeqRef,0);
		
		my %startListCDS = %$ref2StartList;
		my %stopListCDS  = %$ref2StopList; 
		
		print "Aligned sequence Reference is '$alignedSeqRef'\n" if ($verbose);
		print "Aligned sequence Query is     '$alignedSeqQuery'\n" if ($verbose);
		
		# For the first element
		my $firstElement = $seqStartList[0];
		$startListCDS{$firstElement} = 20 if ($startListCDS{$firstElement} != 20 && $tagsList[0] ne "F");
		$startListCDS{$firstElement} = 23 if ($startListCDS{$firstElement} != 23 && $tagsList[0] eq "F");
		
		## Correct start and stopList to account for gaps that happen right at the junction of two exons
		for(my $k = 0; $k < scalar(@exonsList)-1; $k++)
		{
			my $stopPtCurrent  = $stopListCDS{$seqStopList[$k]};
			my $startPtNext    = $startListCDS{$seqStartList[$k+1]};               
			
			$stopListCDS{$seqStopList[$k]} = $startPtNext if ($stopPtCurrent != $startPtNext);
		}
		
		## For the last element
		my $lastElement = $seqStopList[$#seqStopList];
		$stopListCDS{$lastElement} = length($alignedSeqRef) - 20 if ($stopListCDS{$lastElement} != length($alignedSeqRef) - 20);
		$stopListCDS{$lastElement} = length($alignedSeqRef) - 23 if ($stopListCDS{$lastElement} != length($alignedSeqRef) - 23 && $tagsList[$#tagsList] eq "L");
		
		my $cumLengthSub = 0;
		my $copySeq = $alignedSeqQuery;
		$copySeq =~s/-//g;
		my $lPrint = length($copySeq);
		#print "HERE YOU ARE !!!\n\n#### ==> lPrint is '$lPrint'\n";	

		for(my $e = 0; $e < scalar(@exonsList); $e++)
		{
			my $startPt = $startListCDS{$seqStartList[$e]};
            		my $stopPt  = $stopListCDS{$seqStopList[$e]};
		
			print "The sequence will be extracted from '$startPt' (start) until '$stopPt' (stop)\n" if ($verbose);
			my $tag = $tagsList[$e];
			
			## Reference sequence stuff
			my $refSeq = "";
			$refSeq = uc(substr($alignedSeqRef,$startPt,($stopPt-$startPt)));
			$refSeq = reverse($refSeq) if ($strandRef eq "-");
			print "### The aligned sequence for the reference for exon '$exonsList[$e]' is '$refSeq' ###\nThe sLine part is '$sLineFirstPartList[$e]'\n" if ($verbose);
			
			my($acceptorSeqR,$donorSeqR) = (split /#/,$flankingSeqList[$e])[0,1];
			$acceptorSeqR = revComp($acceptorSeqR) if ($strandRef eq "-");
			$donorSeqR    = revComp($donorSeqR) if ($strandRef eq "-");
			my $refSeqLine 		= $sLineFirstPartList[$e]." ".uc($acceptorSeqR.$refSeq.$donorSeqR);
			
			## Query sequence stuff
			my $acceptorSeqQ = my $donorSeqQ = my $querySeq = "";
			if ($e == 0){
				#print "----> Case1\n";
				$querySeq = substr($alignedSeqQuery,$startPt,($stopPt-$startPt));
				$donorSeqQ    = createString(20,"-");
				
				$acceptorSeqQ = substr($alignedSeqQuery,0,20);
				$acceptorSeqQ = substr($alignedSeqQuery,0,23) if ($tag eq "F");
			}elsif ($e == scalar(@exonsList) - 1){
				#print "-----> Case2\n";
				#print "Start point is '$startPt'\n";
				$querySeq = substr($alignedSeqQuery,$startPt,($stopPt-$startPt));
				$acceptorSeqQ = createString(20,"-");
				
				$donorSeqQ = substr($alignedSeqQuery,(length($alignedSeqQuery)-20));
				$donorSeqQ = substr($alignedSeqQuery,(length($alignedSeqQuery)-23)) if ($tagsList[$e] eq "L");
			}else{
				#print "-----> Case3\n";
				$querySeq = substr($alignedSeqQuery,$startPt,($stopPt-$startPt));
				$acceptorSeqQ = createString(20,"-");
				$donorSeqQ    = createString(20,"-");
			}
			
			if ($strandRef eq "-"){
				$querySeq 	= reverse($querySeq);
				$acceptorSeqQ  	= reverse($acceptorSeqQ);
				$donorSeqQ	= reverse($donorSeqQ);
			}
			
			print "\n### The aligned sequence for query for exon '$exonsList[$e]' is '$querySeq' ###\n" if ($verbose);
			
			my $seqPrintQuery = "";
			#print "These are the componenets\n'$acceptorSeqQ'\n'$querySeq'\n'$donorSeqQ'\n";
			$seqPrintQuery = uc($acceptorSeqQ.$querySeq.$donorSeqQ);
			$seqPrintQuery = uc($donorSeqQ.$querySeq.$acceptorSeqQ) if ($strand eq "-");
			
			my $seqPrintQueryCopy = $seqPrintQuery;
			$seqPrintQueryCopy =~s/-//g;
			
			my $startPrint = "";
			my $lengthQuery = length($seqPrintQueryCopy);
			if ($strandRef eq "+"){
				$startPrint = $startCDS + $cumLengthSub;
			}else{  ## The coordinates have to be adjusted if the gene is on the minus strand in the reference. Thankfully, it is not so hard.	
				$startPrint = $startCDS + $lengthQueryTotal - $lengthQuery - $cumLengthSub;
			}
			$cumLengthSub = $cumLengthSub + $lengthQuery;
		
			#print "SLine manipulation --> $species, $chrQuery, $startPrint, $lengthQuery, $strandQuery, $sizeQuery, $seqPrintQuery\n";	
			my $querySeqLine = "s $species.$chrQuery $startPrint $lengthQuery $strandQuery $sizeQuery $seqPrintQuery";
			#print "The sLine part is '$querySeqLine'\n" if ($verbose);
			
			my $tmpMaf = `mktemp /dev/shm/XXXXX.file`; chomp $tmpMaf;
			push(@tempFilesList,$tmpMaf);
			
			open(FOTM,">$tmpMaf") || die "Error writing to tempMaf '$tmpMaf' while doing the checkMafCorrectness test\n";
			print FOTM "$refSeqLine\n$querySeqLine\n";
			close FOTM;

			checkMafCorrectness($tmpMaf,$reference,\@tempFilesList);  ## Sanity check for mafs:
			my $key = $refIDList[$e];
			writeMafToBDB($tmpMaf,$key,$outBDB,\@tempFilesList);
			print "Exon $exonsList[$e] done\n" if ($verbose);
			`rm -rf $tmpMaf`;
		}
	}
	
	## As a conservative approach, still delete everything that is present in tmpFilesList even though there should not be anything
	if (! $inter){  ## But do not do that for interactive mode, we need to see the temp files when 'inter' has been specified
		my $tmpFilesListString = join(" ",@tempFilesList);
		`rm -rf $tmpFilesListString`;
	}
}
close FI;
