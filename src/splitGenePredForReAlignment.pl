#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;

my $usage = "perl $0 inFile outFile outDir lengthCutOffPerJob numberOfExonsPerJob\n####\n
PURPOSE: $0 reads the genePred file (myGenePred format) where each line is an isoform. And then merges different isoforms 
so that only one exon per line is written. All the unique exons from different isoforms of a gene are printed. No redundant 
information/exons that are shared between different isoforms are printed. Subsequently, the exons are binned and the big 
genepred is split into smaller genePreds, these smaller genePreds are used as an input to CESAR\n#####\n";
die $usage if (scalar(@ARGV) != 5);

my $in  			         = $ARGV[0];
my $out 			         = $ARGV[1];
my $outDir      		     = $ARGV[2];
my $lengthCutOffPerJob       = $ARGV[3];
my $maxNumberOfExonsPerJob   = $ARGV[4];

`mkdir -p $outDir`;

## Check input file;
die "Input file '$in' does not exist\n" if (! -e $in);

my %hashAllIsoforms = my %hashChrStrand = ();
my %phaseHash = my %numberHash = ();

open(FI,$in);
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$isoform,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,3,8];
	
	my @exonsList = split(",",$exons);
	@exonsList = reverse(@exonsList) if ($strand eq "-");
	
	### get incoming phase and outgoing phase for each exon here
	my $ref2PhaseList = getExonPhaseList(\@exonsList);
	my @phaseList = @$ref2PhaseList;
	
	## Associate each exon to its phase and its number in the transcript
	for(my $i = 0;$i < scalar(@phaseList); $i++) {
		my $j = $i+1;	
		my $number = "M";
		
		if (scalar(@exonsList) == 1) {	
			$number = "S";
		} elsif ($j == 1) {	
			$number = "F";
		} elsif ($j == scalar(@exonsList)) {
			$number = "L";
		}
		
		my $exon  = $exonsList[$i];
		my $phase = $phaseList[$i];
		
		$phaseHash{$isoform}{$exon}  = $phase; 		## This is important, the same exon could have different phase
		$numberHash{$isoform}{$exon} = $number;		## And even different numbers despite having exactly the same coordinates
	}
	
	$hashAllIsoforms{$gene}{$isoform} = $exons;
	$hashChrStrand{$gene} = "$chr\t$strand";
}
close FI;

my %hashExonsList = ();
my %hashGeneIsoformList = ();

## Now create triplet for each exon which includes the exon in question, plus the exon upstream and the exon downstream. And store everything in a hash
open(FI,$in);
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$isoform,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,3,8];
	
	my @exonsList = split(",",$exons);
	@exonsList = reverse(@exonsList) if ($strand eq "-");
	
	if (scalar(@exonsList) == 1) {
		my $exonUp = "F";
		my $exonDown  = "L";
		my $phaseUp = my $phaseDown = my $numberUp = my $numberDown = "";
		my $string = "$gene\t$isoform\t$phaseUp\t$numberUp\t$exonUp\t$phaseHash{$isoform}{$exonsList[0]}\t$numberHash{$isoform}{$exonsList[0]}\t$exonsList[0]\t$phaseDown\t$numberDown\t$exonDown";
		$hashExonsList{$string} = "T";
		$hashGeneIsoformList{$gene}{$isoform}{$string} = "T";
	} else {
		for(my $j = 0; $j < scalar(@exonsList); $j++) {
			my $exonUp = my $exonDown = my $phaseUp = my $phaseDown = my $numberUp = my $numberDown = "";
			if ($j == 0)	{	$exonUp = "F"; $phaseUp = "NA"; $numberUp = "NA";	}
			else			{	$exonUp = $exonsList[$j-1]; $phaseUp = $phaseHash{$isoform}{$exonUp}; $numberUp = $numberHash{$isoform}{$exonUp}; 	}
			
			if ($j == scalar(@exonsList)-1) {	
				$exonDown = "L";
				$phaseDown = "NA";
				$numberDown = "NA";
			} else {
				$exonDown = $exonsList[$j+1];
				$phaseDown = $phaseHash{$isoform}{$exonDown};
				$numberDown = $numberHash{$isoform}{$exonDown};
			}
			
			my $string = "$gene\t$isoform\t$phaseUp\t$numberUp\t$exonUp\t$phaseHash{$isoform}{$exonsList[$j]}\t$numberHash{$isoform}{$exonsList[$j]}\t$exonsList[$j]\t$phaseDown\t$numberDown\t$exonDown";
			$hashExonsList{$string} = "T";
			$hashGeneIsoformList{$gene}{$isoform}{$string} = "T";  ## A 3D hash which stores the gene name, the isoform name and the exon info (exon info is info about the exon plus the up and the downstream exon)
		}
	}
}
close FI;

# Now print the unique exons list to the exonUnion file. For exons that are common to different transcripts, group the transcripts and print them together
# For exons that are specific to a particular isoform, no grouping (obviously)

open(FO,">$out") || die "Cannot write to outFile '$out'\n";
foreach my $gene(keys(%hashGeneIsoformList)) {
	my %uniqueExons = my @isoformsList = my %isoformHash = ();
	foreach my $isoform(keys (%{$hashGeneIsoformList{$gene}})) {
		push(@isoformsList,$isoform);  ## For every exon create a list of the isoforms that have this exon
		foreach my $exon(keys (%{$hashGeneIsoformList{$gene}{$isoform}})) {
			my @tmp = split(/\t/,$exon);
			splice(@tmp,0,2);  ## Chop off the first 2 fields, the gene name and the isoform
			my $exonString = join("\t",@tmp);
			$uniqueExons{$exonString} = "T"; ## Store all unique exons in this has
			$isoformHash{$isoform}{$exonString} = "T";  ## This hash helps me to figure which exon is contained in which isoform, see below
		}
	}
	
	foreach my $exon(keys(%uniqueExons)) {  ## Loop over all unique exons
		my @exonList = ();
		foreach my $isoform(@isoformsList) {
			push(@exonList,$isoform) if (exists $isoformHash{$isoform}{$exon}); ## Store the isoform in the exonList array if this unique exon is present in the isoform
		}
		my $allIsoforms = join(",",@exonList);  ## Combine all isoforms to form a string
		my($phase,$tag,$cds) = (split /\t/,$exon)[3,4,5]; ## The useful stuff about the exon in question/
		
		print FO "$gene\t$allIsoforms\t$hashChrStrand{$gene}\t$phase\t$tag\t$cds\n";
	}
}
close FO;

### Sort the output file:
my %lengthHash = ();
open(FI,$out);
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my $exonCds = (split /\t/,$line)[6];
	my ($exonStart,$exonStop) = (split /-/, $exonCds)[0,1];
	
	my $length = abs($exonStart - $exonStop);
	$lengthHash{$line} = $length;
}
close FI;

my $fileTmp = `mktemp|tr -d "\n"`;

open(FO,">$fileTmp");
foreach my $keys(sort {$lengthHash{$a} <=> $lengthHash{$b}} keys(%lengthHash) ) {
	print FO "$keys\n";
}
close FO;
`mv $fileTmp $out`;

## The following operation splits the larger genePred into smaller genePreds
my $file = 1;
my $lengthTotal = 0;
my %hashList = my @list = my @lengthList = ();
my $k = 0;

open(FI,$out);
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$isoform,$chr,$strand,$phaseString,$numberString,$uniqueExons) = (split /\t/,$line);
	
	my $length = getLength($uniqueExons);
	$lengthTotal = $lengthTotal + $length;
	
	## The lengthTotal threshold has been determined with some hit and trial attempts. The original choice was 7500 based on a test with 46way alignment. 36 exons whose cumulative length was 10000
	## was bundled as a job. In short queue (with no-resubmission allowed), only 27 exons finished. So I made this threshold to 0.75 times 10000. However, some tests with this number
	## did not work either, so I made this to 5000 first and 4000 later -->  with 4000, all jobs that are supposed to finish in short queue indeed finish in the short queue
	## NOTE: More tests need to be done with hg19-100way alignment or say mm10-60way alignment.
	
	if ( ($lengthTotal > $lengthCutOffPerJob) || ($k >= $maxNumberOfExonsPerJob) ) {	 ## Irrespective of the exons length, do not have more than 70 exons per job
											
		my $string = join("\n",@list); 
		my $queue = getQueue(\@lengthList);	
		
		$hashList{"$file\_$queue"} = $string;
		
		@list = @lengthList = ();
		$k = $lengthTotal = 0;
		$file++;
		
		$lengthTotal = $lengthTotal + $length;
	}
	
	$list[$k] = "$line";
	$lengthList[$k] = $length;
	$k++;
}
close FI;

## The last values have not been assigned to a string because the if statement was not encountered
my $string = join("\n",@list);
my $queue = getQueue(\@lengthList);

$hashList{"$file\_$queue"} = $string;

foreach my $keys( keys(%hashList)) {
	my $string = $hashList{$keys};
	
	my $outFile = "$outDir/$keys";
	open(FO,">$outFile");
	print FO "$string\n";
	close FO;
}

###########################################################
###########################################################
## Sub-routines

sub getExonPhaseList {		## This subroutine gives the phases (5' and 3') of each exon.
	my $ref2ExonsList = shift;
	my @exonsList =  @$ref2ExonsList;
	
	my $i = 0;
	my @exonLength = ();
	my $lengthCum = 0;
	
	my @exonPhaseList = ();
	
	foreach my $exon(@exonsList) {
		my ($start,$stop) = (split /-/,$exon)[0,1];
		my $length = abs($start-$stop);
		$lengthCum = $lengthCum + $length;
		
		$exonLength[$i] = $lengthCum;
		$i++;
	}
	
	for (my $i = 0;$i < scalar(@exonLength); $i++) {
		my $exonLengthCum = $exonLength[$i];
		
		my $phase5Prime = my $phase3Prime = "";
		
		if ($i != 0) {
			$phase5Prime = 3 - ($exonLength[$i-1]%3);
			$phase3Prime = ($exonLength[$i]%3);
		} else {
			$phase5Prime = 0;
			$phase3Prime = ($exonLength[$i]%3);
		}
		
		$phase5Prime = 0 if ($phase5Prime == 3);
		$phase3Prime = 0 if ($phase3Prime == 3);
		
		my $phaseList = "$phase5Prime-$phase3Prime";
		$exonPhaseList[$i] = $phaseList;
	}
	
	return \@exonPhaseList;
}

sub getLength {
	my $exons = shift;
	my @exonsList = split(",",$exons);
	
	my $length = 0;
	foreach my $ex(@exonsList) {
		my ($exStart,$exEnd) = (split /-/,$ex)[0,1];
		$length = $length + abs($exStart - $exEnd);
	}
	
	return $length;
}

sub getQueue {	## This function tells me which queue should I use:
	my $ref2List = shift;
	
	my @lengthList = @$ref2List;
	my @lengthsSorted = sort{$b <=> $a} (@lengthList);		## Arrange in descending order
	
	my $sumFirst5 = 0;		## Look at the first 5 elements of each smaller genePred:
	## if the number of the elements is less than 5 (happens with a very long exon), then simply look at all the elements
	
	my $k = 5;
	$k = scalar(@lengthsSorted) if (scalar(@lengthsSorted) < 5);
	
	for(my $i = 0; $i <= ($k - 1); $i++) {
		$sumFirst5 = $sumFirst5 + $lengthsSorted[$i];
	}
	
	my $queue = "short";	## The default is the short queue
	
	if ($sumFirst5 >= 5000) {	## If the sum of first 5 exons is greater than 5000, let it run in the long queue
		$queue = "long";
	} elsif ( ($sumFirst5 >= 2000) && ($sumFirst5 < 5000) ) { ### if it is between 2000-5000, let it run in the medium queue
		$queue = "medium";
	}
	return $queue;
}
