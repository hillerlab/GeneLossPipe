#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;
use MyKentFunctions;
use MyBDB;

use lib "/home/vsharma/cesarScripts"; #use lib "/scratch/users/vsharma/cleaning_145Way";
use realignCESAR;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my $verbose = my $inter = 0;

my $in = my $bbIndex = my $reference = my $speciesListFile = my $outBDB = my $clade = my $logFile = my $deletedIntronsBDB = my $intronTableDir = my $program = my $parameters = "";
my $masking = 0;  ## If this is specified, then qualityMask the files:

sub usage
{
	die "usage: $0 -in inFile -bbIndex bigBedIndexFile -ref RefereceSpecies -queryList QuerySpeciesList -outBDB outFile.bdb -clade human|mouse|drosophila|celegans 
	-deletedIntronsBDB BDB file containing alignments for exons where the down|upstream intron is deleted -intronTableDir  Directory containing intron-look-up tables 
	-program (CESAR|CESAR2) -parameters [Optional: Only if you do not want to use the default transition probabilities. Works only with cCESAR]
	-log logFile -inter[Optional - useful for bug fixes, keeps temporary files] -v -masking (if you want to quality mask files)\n";
}
GetOptions ("v|verbose"  => \$verbose, "in=s" => \$in, "bbIndex=s" => \$bbIndex, "ref=s" => \$reference, "queryList=s" => \$speciesListFile, "outBDB=s" => \$outBDB, 
"clade=s" => \$clade, "log=s" => \$logFile, "deletedIntronsBDB=s" => \$deletedIntronsBDB, "intronTableDir=s" => \$intronTableDir, "program=s" => \$program, "parameters:s" => \$parameters,
"inter" => \$inter, "masking" => \$masking) || usage();
usage() if ($in eq "");

## Check different parameters:
die "Big-bed index does not exist|is not specified\n" if (! -e $bbIndex || $bbIndex eq "");
die "Reference species not specified\n" if ($reference eq "");
die "QueryList file does not exist|not specified\n" if ($speciesListFile eq "" || ! -e $speciesListFile);
die "Output BDB file is not specified\n" if ($outBDB eq "");
die "Intron look-up tables directory not specified|does not exist\n" if (! -d $intronTableDir || $intronTableDir eq "");
die "deletedIntron BDB directory not specified| does not exist\n" if (! -d $deletedIntronsBDB || $deletedIntronsBDB eq "");
die "The parameter value for program '$program' is not recognized\n" if ($program ne "CESAR" && $program ne "CESAR2");

## All parameters are OK. Time for business:
my @tempFilesList = ();	## This array contains a list of all temp files.

my @speciesList = (); 
open(FI,$speciesListFile) || die "Cannot open speciesList file '$speciesListFile'\n";	## Create an array of species from speciesList
while (my $line = <FI>)
{
	$line =~s/\s+$//;
	push(@speciesList,$line);
}
close FI;
unshift(@speciesList,$reference);		## Add reference species to the speciesList if not already present
my $speciesString = join(",",@speciesList);

my $speciesListFileAll = `mktemp /dev/shm/XXXXXX.file`; chomp $speciesListFileAll;
push(@tempFilesList,$speciesListFileAll);
open(FO,">$speciesListFileAll");
foreach my $species(@speciesList)
{
	print FO "$species\n";
}
close FO;

## Check if all species have 2 bit files, intron-lookUp tables as well as bdbFiles that contain alignments for those exons which have a deleted up/downstream intron.
my $qualityMaskedFilesList = `mktemp /dev/shm/XXXXX.file`; chomp $qualityMaskedFilesList;  ## This file stores the location of the qualityMasked files, used later by mafQualityMasking.pl
push(@tempFilesList,$qualityMaskedFilesList);

open(FO,">$qualityMaskedFilesList") || die "Cannot write to qualityMaskedFilesList '$qualityMaskedFilesList'\n";
foreach my $species(@speciesList)
{
	$species =~s/\s+$//;
	my $queryFile  = "$ENV{'genomePath'}/gbdb-HL/$species/$species.2bit";
	my $twoBitFileQueryMasked = "$ENV{'genomePath'}/gbdb-HL/$species/$species.qualityStrict.2bit";
   	$twoBitFileQueryMasked    = "$ENV{'genomePath'}/gbdb-HL/$species/$species.quality.2bit" if(! -e $twoBitFileQueryMasked);	## Just use the species.2bit file if no qualityStrict2bit file is available
	die "The species '$species' does not have a qualityStrict.2bit or a quality.2bit file\n" if (! -e $twoBitFileQueryMasked);
	
	print FO "$species\t$twoBitFileQueryMasked\n";
	die "2 bit file i.e. '$queryFile' not present for '$species'\n" if (! -e $queryFile);
	if ($species ne $reference) {
		my $intronTable1 = "$intronTableDir/$species.5Prime.BDB";
		my $intronTable2 = "$intronTableDir/$species.3Prime.BDB";
		my $deletedIntronExonsAlignment = "$deletedIntronsBDB/$species.realigned.BDB";
		
		die "Intron look-up tables missing for '$species'\n" if (! -e $intronTable1 || ! -e $intronTable2);
		#print "WARNING!! The BDB file containing alignments for those exons whose introns are deleted is missing for '$species' \n" if (! -e $deletedIntronExonsAlignment);
	}
}
close FO;

@speciesList = grep {$_ ne $reference} @speciesList;	## Now remove the reference species from the list
my %ignoreExonsHash = (); ## Store all exons that should be ignored in a hash
foreach my $species(@speciesList)
{
	my $file = "$intronTableDir/$species.ignoreExons";
	if (-e $file) {
		open(FI,$file) || die "Error opening ignore exons file '$file'\n";
		while (my $line = <FI>)
		{
			my($gene,$chr,$cds) = (split /\t/,$line)[0,2,3];
			my $key = "$gene\_$chr\_$cds";
			$ignoreExonsHash{$species}{$key} = "T";
		}
		close FI;
	} else {
		#print "WARNING!! No ignore exons file found for '$species'\n";
	}
}

my $twoBitFileRef = "$ENV{'genomePath'}/gbdb-HL/$reference/$reference.2bit";

my %incompCds = ();

## pre-process the input file
my $tmpIn = `mktemp`; chomp $tmpIn;

open(FI,$in) || die "Error opening input exons list '$in'\n";
open(FO,">$tmpIn");
while (my $line = <FI>)
{
	$line =~s/\s+$//;
	my @tmp = split(/\t/,$line);
	my $line1 = shift(@tmp);
	my $tList = shift(@tmp);
	my $line2 = join("\t",@tmp);
	
	my @transcriptsList = split(/,/,$tList);

	foreach my $t(@transcriptsList)
	{
		print FO "$line1\t$tList\t$t\t$line2\n" if (! exists $incompCds{$t});
	}
}
close FI;
close FO;

$in = $tmpIn;  # Use the processed file as the input;

open(FI,$in) || die "Error opening the input file '$in'\n";
open(FOLOG,">>$logFile") || die "Error!! Cannot write to the logFile '$logFile'\n";  ## You append to the logFile, because in case of a resubmitted job, you do not re-do the computation that you have already done for an exon
while (my $line = <FI>)
{
	$line =~s/\s+$//;
	print "Aligning the exon --> '$line'\n";  ## #IL13RA2	ENST00000243213,ENST00000371936	chrX	-	1-1	M	114245206-114245391  ## This is how an input line looks like
	

	my($gene,$transcriptIDs,$isoform,$chr,$strand,$phase,$tag,$exon) = (split /\t/, $line)[0,1,2,3,4,5,6,7];
	my @transcriptsList = split(",",$transcriptIDs);
	
	my $keyOri = "$gene\_$chr\_$exon"; ## Key for the ignoreExons hash;
	my ($phase5Prime,$phase3Prime) = (split /-/,$phase)[0,1];	
	my ($exonStart,$exonStop) = (split "-",$exon)[0,1];
	
	## Check if exon in the reference is not at the boundaries of a scaffold/chromsome. If yes, skip this exon. For example:
	##Gene19590	scaffold_3496	-	0-0	S	2-2078
	my $tmpStart = $exonStart - 25;
	my $tmpStop  = $exonStop + 25;
	my $chrSize = getChromSize($reference,$chr,\@tempFilesList);
	if ($tmpStart < 0 || $tmpStop > $chrSize){
		print FOLOG "Skipped $line\tThis exon lies at the beginning/end of the scaffold\nStart is '$tmpStart', stop is '$tmpStop' while chrSize is '$chrSize'\n";
		next;
	}
	
	my $refID = "$reference#$transcriptIDs#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag";             ## This id will serve as the key in the BDB file
	my $refIDIsoform = "$reference#$isoform#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag";      ## This id stores the alignment obtained by running RealignmentPairs_Pipeline.pl
	my $refIDIsoform_Table = "$gene#$isoform#$chr#$exonStart-$exonStop";
	
	my $exonDone = 0;
	if (-e $outBDB){
		foreach my $isoform(@transcriptsList)
		{
			my $key = "$reference#$isoform#$exonStart#$exonStop#$chr#$strand#$phase5Prime#$phase3Prime#$tag";
			my $mafValue = readBDB($outBDB,$key);
			$exonDone++ if ($mafValue ne "");
		}
	}

	if ($exonDone != scalar(@transcriptsList)) {	## Do this only if the keys are not present in the bdb file already.
							## This helps in case of resubmitted jobs and prevents re-aligning exons that have been realigned once.
		if (abs($exonStart-$exonStop) <= 6) { ## Ignore those exons which are less than 6 base pairs
			my $maf = `mktemp /dev/shm/XXXXX.file`; chomp $maf;
			unshift(@tempFilesList,$maf); ## Every temp file goes to this array
				
			$exonStart = $exonStart - 20;
			$exonStop  = $exonStop + 20;
			my $statusME = system("mafExtract -region=$chr:$exonStart-$exonStop $bbIndex stdout|mafSpeciesSubset stdin speciesList=$speciesString species.lst=NULL $maf");  ## Call mafextract 
			safeDie("Error running mafExtract for '$exon'. The index was '$bbIndex'\n",\@tempFilesList) if ($statusME != 0);

			writeMafToBDB($maf,$refID,$outBDB,\@tempFilesList); ### And write this maf to bdb
			next;
		}
	
		($exonStart,$exonStop) = fixCoordinates($exonStart,$exonStop,$strand,$tag) if ($tag ne "M"); ## Fix these coordinates right here for non-middle exons
		my $lengthExon = $exonStop - $exonStart;
		
		## get sequence for the reference species:
		my $seqRefPrint = getSequenceFrom2Bit($twoBitFileRef,$chr,$exonStart,$exonStop,\@tempFilesList);
		$seqRefPrint = revComp($seqRefPrint) if ($strand eq "-");	## this is the sequence that is used for alignment, but its revComp is written to maf if the gene comes from minus strand
		$seqRefPrint = fixPhase($seqRefPrint,$phase5Prime,$phase3Prime,$program);
		
		## Temporary bug fix for CESAR since CESAR fails with those exons where the split codon base contains "n". Discovered this while running CESAR on HLtupMer3-30way
		my $refSeqWithN = $seqRefPrint;
		$seqRefPrint =~s/n/c/g;
		## Locate positions which had the 'n' character
		my @smallNPositions = ();
		my $nPos = 0;
		while ($nPos < length($refSeqWithN))
		{
			my $char = substr($refSeqWithN,$nPos,1);
			push(@smallNPositions,$nPos) if ($char eq "n");
			$nPos++;
		}
		
		my $seqRefOri = $seqRefPrint;	## Make a copy of Reference sequence here;
		### Replace any in-frame stop codon in the Reference sequence with NNN
		my ($seqRefPrintNew,$list2CodonsFlipped) = flipStopCodons($seqRefPrint,$phase5Prime,$phase3Prime);
		$seqRefPrint = $seqRefPrintNew;
		
		my $flankingSeqString = getFlankingSeq($exonStart,$exonStop,$chr,$strand,$twoBitFileRef,$tag,20,\@tempFilesList); ## Seq upstream and downstream of the exon
		### Check if the exon is flanked by a U12 intron
		my $acceptorAC = my $donorAT = "";
		my ($accSeq,$donorSeq) = (split "#",$flankingSeqString)[0,1];
		$acceptorAC = "T" if ($accSeq =~/ac$/);
		$donorAT    = "T" if ($donorSeq =~/^at/);
		
		my $refSeq_1 = getRefSeqForMaf($seqRefPrint,$refID,$flankingSeqString,$strand,$tag);		## This gives me the first part of the s-Line for the Reference species
		
#### Now run mafExtract, get the maf. From this maf, get the coordinates of the aligning locus for different species. From these coordinates, get the sequence plus/minus 50bp
#### and write it to a fasta file
		
		my $maf = `mktemp /dev/shm/XXXXX.file`; chomp $maf;
		unshift(@tempFilesList,$maf); ## Every temp file goes to this array
		
		my $statusME = system("mafExtract -region=$chr:$exonStart-$exonStop $bbIndex $maf");  ## Call mafextract	
		safeDie("Error running mafExtract for '$exon'. The index was '$bbIndex'\n",\@tempFilesList) if ($statusME != 0);
		print "Mafextract run successful for $refID\n" if ($verbose);
		
		## Get the number of lines in the maf:
		my $nMaf = `wc -l < $maf`; chomp $nMaf;
		
		my %eLinesHash = ();
		if ($nMaf == 1) { ## This could happen, when none of the species in the maf have something that aligns to the reference species. In such cases, mafExtract returns nothing
						## An example: mafExtract /projects/genome/gbdb-HL/ce10/maf/ce10_7wayCleaned/multiCleaned.bb -region=chrX:4646151-4646369 stdout
			writeOnlyReferenceSequence($reference,$flankingSeqString,$seqRefPrint,$refSeq_1,$strand,$refID,$outBDB,$seqRefOri,$list2CodonsFlipped,$phase5Prime,$phase3Prime,$masking,\%eLinesHash,\@speciesList,\@tempFilesList);
			next;	
		}
		
		### Otherwise create a fasta file that will be input to CESAR
		my $fasta = `mktemp /dev/shm/XXXXX.file`; chomp $fasta;
		unshift(@tempFilesList,$fasta); ## Every temp file goes to this array
		open(FO,">$fasta") || safeDie("Cannot write to fasta file '$fasta'\n",\@tempFilesList);
		print FO ">$refID\n$seqRefPrint\n" if ($program eq "CESAR");
		print FO ">$refID\n$seqRefPrint\n####\n" if ($program eq "CESAR2");  ## CESAR2 requires that after the end of the reference sequence, I have a string with 4 hashes.
		my $lengthRefSequence = length($seqRefPrint);
		
		my $nSeqFasta = 1;
		
		my %shortEnd = my %alignmentFromCombinedExonHash = ();
		my %shortEnd5P = my %shortEnd3P = ();
		my %flank5PrimeZero = my %flank3PrimeZero = ();
		print "Now filtering maf for individual species\n" if ($verbose);
		
		#my $sleepTime = int(rand(10));  ## Sleep for a few seconds, this lets me ensure that not all jobs are accessing the same binary and the same data (particularly twoBitToFa and the twoBit files) at the same time. Only need for furiosa.
		#sleep $sleepTime;
		
		my $allPairWiseAlignment = `mktemp /dev/shm/XXXXX.file`; chomp $allPairWiseAlignment;
		unshift(@tempFilesList,$allPairWiseAlignment);
		
		open(FOAP,">$allPairWiseAlignment") || safeDie("Error running to temp allPairWiseFile '$allPairWiseAlignment'",\@tempFilesList);
		foreach my $species(@speciesList)
		{
			$species =~s/\s+$//;
			next if ($ignoreExonsHash{$species}{$keyOri});
			my $deletedIntronExonsAlignment = "$deletedIntronsBDB/$species.realigned.BDB";
			
			if (-e $deletedIntronExonsAlignment) {
				my $alignmentIntronDel = readBDB($deletedIntronExonsAlignment,$refIDIsoform);
				if ($alignmentIntronDel ne ""){
					print "--> Extracting alignment for '$refIDIsoform' for '$species' from pre-existing grouped alignment\n" if ($verbose);
					print FOAP "##maf version=1 scoring=tba.v8\n";
					print FOAP "$alignmentIntronDel";
					print FOLOG "$line\t$species\tAlignment extracted from the output of RealignmentPipeline_Pairs.pl\n";
					next;
				}
			} else {
				#print "WARNING!!! The deletedIntronExonAlignment file '$deletedIntronExonsAlignment'  was not found. Did you forget to run Realignment_ExonPairs.pl?\n";	
			}
			
			my($ref2cdsStartQuery,$ref2cdsStopQuery,$ref2strandQuery,$ref2chrQuery,$ref2moreThanOneStrand,$ref2moreThanOneChr,$ref2eLinePrint,$ref2SpeciesList) = getCoordinatesFromMaf($maf,$reference,\@tempFilesList);
			my %speciesPresentHash = %$ref2SpeciesList;
			
			if (exists $speciesPresentHash{$species}) { ## Check if the query species is present in the filtered maf
				## Next get the coordinates
				my $cdsStartQuery 		= $ref2cdsStartQuery->{$species};
				my $cdsStopQuery 		= $ref2cdsStopQuery->{$species};
				my $strandQuery 		= $ref2strandQuery->{$species};
				my $chrQuery 			= $ref2chrQuery->{$species};
				my $moreThanOneStrand 	= $ref2moreThanOneStrand->{$species};
				my $moreThanOneChr 		= $ref2moreThanOneChr->{$species};
				my $eLine 				= $ref2eLinePrint->{$species};
				
				print "### Aligning locus for the query species --> $cdsStartQuery,$cdsStopQuery,$strandQuery,$chrQuery,$moreThanOneStrand,$moreThanOneChr,$eLine\n\n" if ($verbose);
				
				if ($cdsStartQuery eq "NA" && $cdsStopQuery eq "NA") {	## Filter1: This will only happen if the $mafOut contains only the e-lines 
					$eLinesHash{$species} = $eLine;
					print FOLOG "$line\t$species\tELines\n";
					next;
				} elsif ($moreThanOneStrand eq "T" || $moreThanOneChr eq "T" ){ ## Filter2: Do nothing for these exons. We will treat such loci as missing data for the moment
					print FOLOG "$line\t$species\tThe exon comes from morethanOneStand '$moreThanOneStrand' or moreThanOneChr '$moreThanOneChr'\n";
					next;
				} elsif ($cdsStartQuery > $cdsStopQuery){		## Filter3: This could happen when an exon undergoes translocation|inversion
					print FOLOG "$line\t$species\tThe cdsStart '$cdsStartQuery' is greater than cdsStop '$cdsStopQuery'\n";
					next;
				} elsif (abs($cdsStartQuery - $cdsStopQuery) - length($seqRefPrint) > 10000){	## Filter4: This basically means that I do not allow an insertion of more than ~10,000bp in an exon
														## If that is the case, I ignore the exon --> Fair enough
					print FOLOG "$line\t$species\tThe exon has an insertion greater than 1000 bp'\n";
					next;
				}
				
				my $lengthExonQuery = $cdsStopQuery - $cdsStartQuery;
				my $refIDIsoform_Species = "$refIDIsoform_Table#$species";
				my $fivePrimeBDB  = "$intronTableDir/$species.5Prime.BDB";
				my $threePrimeBDB = "$intronTableDir/$species.3Prime.BDB"; 
				
				my $intronLength5Prime = readBDB($fivePrimeBDB,$refIDIsoform_Species);
				my $intronLength3Prime = readBDB($threePrimeBDB,$refIDIsoform_Species);
				
				$intronLength5Prime = 50 if ($tag eq "F");  ## For the first exon, you can not have a length in 5' direction from the look-up table
				$intronLength3Prime = 50 if ($tag eq "L");  ## Same for the last exon
				$intronLength5Prime = $intronLength3Prime = 50 if ($tag eq "S");  ## And for a single exon gene
				$intronLength5Prime = 50 if ($intronLength5Prime eq ""); ## if the exon in question is exon#2 and the first exon is deleted, we cannot get any value for intronLength here.
				$intronLength3Prime = 50 if ($intronLength3Prime eq ""); ## if the exon in question is exon#12 and the last exon (13th) is deleted, we cannot get any value for intronLength here.
				
				### Now set the default flankingSequence length = 50
				my $flankDef = 50;
				$flankDef = $flankDef + ($lengthExon - $lengthExonQuery) if ($lengthExon > $lengthExonQuery);
				
				my $flank5Prime = my $flank3Prime = $flankDef; ## By default, we take the length of the exon as the length of the flanking sequence (very conservative)
				$flank5Prime = $intronLength5Prime  if ($flankDef > $intronLength5Prime);  ## If the exon is longer than the flank length, then we just use the flank length
				$flank3Prime = $intronLength3Prime  if ($flankDef > $intronLength3Prime);
		
				print "Intron length5P is '$intronLength5Prime', Intron length3P is '$intronLength3Prime', Flank length5P is '$flank5Prime' and Flank length3P is '$flank3Prime' for '$species'\n" if ($verbose);
		
				my $startString = my $stopString = "";
				my $nAtStart = my $nAtEnd = "";
				my $cdsQueryStart = my $cdsQueryStop = "";

				if ($strandQuery eq "+") {  ## These tell you the actual coordinates from where the sequence should be extracted
					$cdsQueryStart = $cdsStartQuery - $flank5Prime;	
					$cdsQueryStop  = $cdsStopQuery + $flank3Prime;
				} else {
					$cdsQueryStart = $cdsStartQuery - $flank3Prime;
					$cdsQueryStop  = $cdsStopQuery + $flank5Prime;
				}
		
				my $chromSize = getChromSize($species,$chrQuery);
				my $shortEnd = "";
		
				if ($cdsQueryStart < 0) {
					$startString = createString(abs($cdsQueryStart),"N");
					$nAtStart = abs($cdsQueryStart);
					$cdsQueryStart = 1;
					$shortEnd = "T";
				}
		
				if ($cdsQueryStop > $chromSize) {
					my $diff = $cdsQueryStop - $chromSize;
					$stopString = createString($diff,"N");
					$cdsQueryStop = $chromSize;
					$shortEnd = "T";
					$nAtEnd = $diff;
				}
				
				## Now run 2BitToFa to extract the sequence
				my $twoBitFileQuery = "$ENV{'genomePath'}/gbdb-HL/$species/$species.2bit"; ## Do not worry about the the qualityMasked file here, we qualityMask the final maf later.
				#$twoBitFileQuery = "orcOrc1.quality.2bit" if ($species eq "orcOrc1");  ## Delete this later
				my $seqQueryPrint   = getSequenceFrom2Bit($twoBitFileQuery,$chrQuery,$cdsQueryStart,$cdsQueryStop,\@tempFilesList);
				
				## if the sequence has Ns at the beginning or the end (assembly gaps), strip them right here and adjust the coordinates accordingly
				my $seqQueryCopy1 = $seqQueryPrint;
				my $seqQueryCopy2 = $seqQueryPrint;
				
				$seqQueryCopy1 =~s/^N+//;
				$cdsQueryStart = $cdsQueryStart + (length($seqQueryPrint)  - length($seqQueryCopy1)); ## adjusting coordinates
				$seqQueryCopy2 =~s/N+$//;
				$cdsQueryStop = $cdsQueryStop - (length($seqQueryPrint)  - length($seqQueryCopy2)); ## adjusting coordinates
				
				$seqQueryPrint =~s/^N+//;
				$seqQueryPrint =~s/N+$//;
				$seqQueryPrint = revComp($seqQueryPrint) if ($strand ne $strandQuery); 
				
				print "Fetched fasta sequence for '$species' using the following coordinates: $species#$cdsQueryStart#$cdsQueryStop#$chrQuery#$strandQuery\n" if ($verbose);
				## Fix the strand orientation
				my $queryID = "$species#$cdsQueryStart#$cdsQueryStop#$chrQuery#$strandQuery";

				$seqQueryPrint = $startString.$seqQueryPrint.$stopString;
				$seqQueryPrint = uc($seqQueryPrint);
				my $lengthQuerySequence = length($seqQueryPrint);
				
				$shortEnd5P{$species} = $nAtStart if ($nAtStart ne "");
				$shortEnd3P{$species} = $nAtEnd if ($nAtEnd ne "");
				
				print FO ">$queryID\n$seqQueryPrint\n";
				$nSeqFasta++;
			}
		}
		close FO;	
		`rm -rf $maf`;
		close FOAP;
		
		if ($nSeqFasta == 1) { ## This could happen, when there was say sequence from 1 aligning species in the maf, so in total you have 2 sequences in the fasta - the reference and 1 query
							 ## However, this sequence was omitted because of Filter4, (see above - an insertion greater than 10,000 bp) in this query sequence
			writeOnlyReferenceSequence($reference,$flankingSeqString,$seqRefPrint,$refSeq_1,$strand,$refID,$outBDB,$seqRefOri,$list2CodonsFlipped,$phase5Prime,$phase3Prime,$masking,\%eLinesHash,\@speciesList,\@tempFilesList);
			next;
		}
		
		print "Maf processing done and fasta file '$fasta' successfully created. Now running CESAR.py on '$fasta'\n"  if ($verbose);
		## At this point, I am done with the maf operations, all the relevant information from the maf has gone to my fasta file.
		## Now run our aligner on the fasta file with different arguments depending on the index/tag of the exon:
		my $out = `mktemp /dev/shm/XXXXX.file`; chomp $out;
		unshift(@tempFilesList,$out); ## Every temp file goes to this array
	
		my $call = "";
		## Specify the acc/donor profiles based on the exon that you are looking at:
		my $accProfile = "extra/tables/$clade/acc_profile.txt";
		my $donorProfile = "extra/tables/$clade/do_profile.txt";
		
		if ($program eq "CESAR") {
		
			if ($tag eq "F"){ 
				$call = "CESAR --is_first_exon IS_FIRST_EXON --clade $clade $fasta > $out";
				$call = "CESAR --is_first_exon IS_FIRST_EXON --has_AT_donor HAS_AT_DONOR --clade $clade $fasta > $out" if ($donorAT eq "T");
			}elsif ($tag eq "L"){
				$call = "CESAR --is_last_exon IS_LAST_EXON --clade $clade $fasta > $out";
				$call = "CESAR --is_last_exon IS_LAST_EXON --has_AC_acceptor HAS_AC_ACCEPTOR --clade $clade $fasta > $out" if ($acceptorAC eq "T");
			}elsif ($tag eq "S"){
				$call = "CESAR --is_first_exon IS_FIRST_EXON --is_last_exon IS_LAST_EXON --clade $clade $fasta > $out";
			}else{
				$call = "CESAR --clade $clade $fasta > $out";		
				$call = "CESAR --has_AT_donor HAS_AT_DONOR --clade $clade $fasta > $out" if ($donorAT eq "T");
				$call = "CESAR --has_AC_acceptor HAS_AC_ACCEPTOR --clade $clade $fasta > $out" if ($acceptorAC eq "T");
			}
		} else {
		
			if ($tag eq "F"){ 
				$accProfile = "extra/tables/$clade/firstCodon_profile.txt";
			}elsif ($tag eq "L"){
				$donorProfile = "extra/tables/$clade/lastCodon_profile.txt";
			}elsif ($tag eq "S"){
				$accProfile = "extra/tables/$clade/firstCodon_profile.txt";
				$donorProfile = "extra/tables/$clade/lastCodon_profile.txt";
			}
				
			$donorProfile = "extra/tables/$clade/u12_donor_profile.txt" if ($donorAT eq "T");
			$accProfile = "extra/tables/$clade/u12_acceptor_profile.txt" if ($acceptorAC eq "T");
			$call = "CESAR2 $fasta --max-memory 300GB -m extra/tables/$clade/eth_codon_sub.txt -p $accProfile $donorProfile --set i3_i1_acc=0.2 i3_i1_do=0.2 cd_acc=0.2 cd_do=0.2 > $out";
			$call = "CESAR2 $fasta --max-memory 300GB -m extra/tables/$clade/eth_codon_sub.txt -p $accProfile $donorProfile --set $parameters --set i3_i1_acc=0.2 i3_i1_do=0.2 cd_acc=0.2 cd_do=0.2 > $out" if ($parameters ne "");  ## Append the parameter string to the call
		}
		
		#print "$call\n"; exit;
		
		my $statusAligner = system($call);
		safeDie("Error running the system call '$call'\n",\@tempFilesList) if ($statusAligner != 0);
		print "Aligment file '$out' produced\n" if ($verbose);
		
		### Now parse this alignment file and make different pairwise mafs.
		open(FIA,$out) || safeDie("Cannot open the output of CESAR --> '$out'\n",\@tempFilesList);
		my @alignmentArray = <FIA>;
		close FIA;
		
		## Now iterate through this array
		for(my $k = 0; $k < scalar(@alignmentArray); $k++)
		{
			if ($alignmentArray[$k] =~/>referenceExon/ || $alignmentArray[$k] =~/>$reference/)
			{
				my $queryLine  = $alignmentArray[$k+2];
				$queryLine =~s/^>//;
				$queryLine =~s/\s+$//;
				my $species = (split /#/,$queryLine)[0];
				
				my $refSeq  	= $alignmentArray[$k+1];
				$refSeq 	=~s/\s+$//;
				my $refSeqCopy = $refSeq;
				$refSeq =~s/^\s+//;
				
				## Flip back the small c to n
				my $ref2Hash_SeqCds2AlignmentCds = getAlignmentCdsHash(\@smallNPositions,$refSeq,0);
				$refSeq = flipBackN($refSeq,$ref2Hash_SeqCds2AlignmentCds);
				my $diff = length($refSeqCopy) - length($refSeq);
				$refSeq = createString($diff," ").$refSeq;
				
				my $querySeq	= $alignmentArray[$k+3];
				$refSeq 	=~s/\s+$//;
				$querySeq	=~s/\s+$//;
				
				$querySeq = replaceNAtEnds($querySeq,$shortEnd5P{$species},"5P") if (exists $shortEnd5P{$species});
				$querySeq = replaceNAtEnds($querySeq,$shortEnd3P{$species},"3P") if (exists $shortEnd3P{$species});
				
				createPairWiseMaf($reference,$refSeq,$refSeq_1,$querySeq,$queryLine,$flankingSeqString,$strand,$phase5Prime,$phase3Prime,$tag,$strand,$seqRefOri,$list2CodonsFlipped,$species,$allPairWiseAlignment,\@tempFilesList);
				print "Pairwise maf created for '$species' and written to '$allPairWiseAlignment'\n" if ($verbose);
			}
		}
		
		if ($inter){
			print "Fasta file is $fasta\nAlignment file is $out\n";
		}else{
			`rm -rf $fasta $out`;	## Delete the fasta and the alignment file
		}
		
		#print "$allPairWiseAlignment\n"; exit;
		createMultiMaf($reference,$allPairWiseAlignment,\%eLinesHash,\@speciesList,$outBDB,$refID,$masking,\@tempFilesList);	### This function creates a multimaf from different pairwise mafs
	}
	
	## As a conservative approach, still delete everything that is present in tmpFilesList even though there should not be anything
	if (! $inter)  ## But do not do that for interactive mode, we need to see the temp files
	{
		my $tmpFilesListString = join(" ",@tempFilesList);
		`rm -rf $tmpFilesListString`;
	}
}
close FI;
close FOLOG;

`rm -f $qualityMaskedFilesList`;
`rm -f $speciesListFileAll`;



