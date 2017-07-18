#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;
use MyKentFunctions;
use MyBDB;

my $geneLossPipe = $ENV{'GeneLossPipeCode'};
use lib "$ENV{'GeneLossPipeCode'}";
use TabularFileFunctions;
use useful_functions_VS;
my $configFile = "$ENV{'GeneLossPipeCode'}/GeneLoss.config";

$| = 1;		# == fflush(stdout)

sub usage
{
	die "usage: $0 [-v|verebose] -input[M] MyGenePredFormatFile -ref[M] RefSpecies -out[M] OutputDirectory -log[M] logFile -index[M] (could either be an indexed maf(with a .bb extension) or a bdb file corresponding to your maf (with a .bdb extension)) -alignment[O] (46way/100way/60way) -species[O] SpeciesList(comma separated) -pMammals[O] -sameSS[O] -excludeGenes[O] (a list of genes that should be excluded, one gene per line) \nParameters with an '[M]' against them are mandatory. Parameters with an '[O]' against them are optional\n";
}

my $verbose = my $sameSSFlag = 0;
my $input = my $alignment = my $RefSpecies = my $data_path = my $speciesL = my $pMammals = my $strict = my $logFile = my $index = my $excludeGenes = "";
GetOptions ("v|verbose"  => \$verbose, "input=s" => \$input, "ref=s" =>\$RefSpecies, "alignment=s" =>\$alignment, "out=s" => \$data_path, "species:s" => \$speciesL, "pMammals" => \$pMammals, "sameSS" => \$sameSSFlag, "strictQM" => \$strict, "log=s" => \$logFile, "index=s" => \$index, "excludeGenes:s" => \$excludeGenes) || usage();
usage() if ($input eq "");

### Input is always a file in myGenePred format.
## Check if the input file is in the right format
open(FI,$input) || die "Error opening the input file '$input'\n";
while (my $line = <FI>)
{
	$line =~s/\s+$//;
	my ($chr,$strand) = (split /\t/,$line)[2,3];
	
	print "Warning: This '$chr' does not look like a valid chromosome value\n" if ($chr !~/chr/ && $chr !~/scaffold/);
	die "This '$strand' is not a valid strand value\n" if ( ($strand ne "+") && ($strand ne "-") );
}
close FI;

## Checking for the remaning parameters

die "Output directory not specified\n" if ($data_path eq "");
die "logFile not specified\n" if ($logFile eq "");
die "Index '$index' not specified|does not exist|does not look either a bdb or a big-bed file\n" if ( ($index eq "") || (! -e $index) || ( ($index !~/\.bdb/i) && ($index !~/\.bb/) ) );

## Get species list:
my @includeSpeciesList = ();
if ($pMammals)		## Filter out only the placental mammals
{
	my $speciesTerm = "$RefSpecies\_pMammals\_$alignment";
	die "Alignment '$alignment' not recognised\n" if  ($alignment !~ /way/ );
	
	open(FI,$configFile) || die "Error opening config file which is $configFile\n";
	while (my $line  = <FI>)
	{
		$line =~s/\s+$//;
		my ($species,$term) = (split /\t/,$line)[0,1];
		push(@includeSpeciesList,$species) if ($term eq $speciesTerm);
	}
	close FI;
}elsif ($speciesL ne ""){
	@includeSpeciesList = split(",",$speciesL);
}else{ ## All the species that "should be" present in the maf becomes my includeSpeciesList then.
	
	my $speciesTerm = "$RefSpecies\_$alignment";
	die "Alignment '$alignment' not recognised\n" if  ($alignment !~ /way/ );
	
	open(FI,$configFile) || die "Error opening config file which is $configFile\n";
	while (my $line  = <FI>)
	{
		$line =~s/\s+$//;
		my ($species,$term) = (split /\t/,$line)[0,1];
		push(@includeSpeciesList,$species) if ($term eq $speciesTerm);
	}
	close FI;
}
## Everything looks OK, time for business

my %excludeGenesHash = ();
if ($excludeGenes ne ""){
	open(FI,$excludeGenes) || die "Error opening excludeGenesList '$excludeGenes'\n";
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		$excludeGenesHash{$line} = "T";
	}
	close FI;
}

my @tmpFilesList = (); ## This array stores all temporary files.

## Create all temporary files/directories only once. And over-write them at every iteration of the loop.

my $tmpFileInput = `mktemp /dev/shm/XXXXXXXX.file`; chomp $tmpFileInput;
my $tempInputDir = `mktemp -d /dev/shm/XXXXX.file`; chomp $tempInputDir;
my $mafFile = `mktemp /dev/shm/XXXXX.file`; chomp $mafFile;
my $mafFileOut = `mktemp /dev/shm/XXXXX.file`; chomp $mafFileOut;
my $tmpBed = `mktemp /dev/shm/XXXXX.file`; chomp $tmpBed;
my $tmpMaf = `mktemp /dev/shm/XXXXX.file`; chomp $tmpMaf;
push(@tmpFilesList,$tmpFileInput,$tempInputDir,$mafFile,$mafFileOut,$tmpBed,$tmpMaf);

open(FIL,$input);
open(FOL,">>$logFile") || die "Cannot write to logFile\n";
while (my $line = <FIL>)
{
	my $sleepTime = int(rand(5));  ## Sleep for a randomTime upto 5 seconds, this helps on furiosa
	sleep $sleepTime;

	$line =~s/\s+$//;
	my @tmp = split(/\t/,$line);
	
	open(FOT,">$tmpFileInput");
	print FOT "$line\n";
	close FOT;
	
	my ($geneName,$accession,$chr,$strand_input,$posStart,$posStop,$exons) = (split /\t/,$line)[0,1,2,3,6,7,8];
	next if (exists $excludeGenesHash{$accession});   
	
	if ($strand_input eq "+"){	
		$posStart++;	
	}else{	
		$posStart = $tmp[7];	
		$posStop  = $tmp[6]+1;
	}
	
	@includeSpeciesList = grep {$_ ne $RefSpecies} @includeSpeciesList;  ## remove reference species list from the list
	unshift(@includeSpeciesList,$RefSpecies); ### In case the list does not have the reference species, add it at the beginning of the list
	my %includeSpeciesHash  = map{ $_ => 1} @includeSpeciesList;	## create a hash where every species in the includeSpeciesList is a key
	my $includeSpeciesString = join(",",@includeSpeciesList);		## And also a comma separated string, useful for mafSpeciesSubset
	
	### First create the bed and the nature files: I need the bedFile only if the input is a mafIndex, I need the nature file for both kinds of inputs though
	
	
	my $statusBNF = system("perl $geneLossPipe/createBNFiles.pl -in $tmpFileInput -out $tempInputDir -length 20");
	safeDie(\@tmpFilesList,"Cannot run createBNFiles on $line\n") if ($statusBNF != 0);
	
	my $exon_ct_file = "$tempInputDir/$accession.nature";
	my $bedFile 	 = "$tempInputDir/$accession.bed";
	safeDie(\@tmpFilesList,"bedFile with exonic coordinates plus/minus 20bp '$bedFile' is absent\n") if (! -e $bedFile);
	safeDie(\@tmpFilesList,"Genomic coordinate annotation file $exon_ct_file is absent\n") if (! -e $exon_ct_file);
	
	### Get maf for all the exons from the either the bdb file (most likely) or from the bbIndex file (less likely):
	
	if ($index =~/\.bdb$/i){	## i.e the index is a BDB file
		my @exonsList = split(/,/,$exons);
		
		### Get 5' and 3'phase for exon - we need this because the key for the exon in the BDB file needs this information
		my @exonsListCopy = @exonsList;
		@exonsListCopy = reverse(@exonsListCopy) if ($strand_input eq "-");  ## This is needed because the gene is on minus strand
		my $ref2PhaseList = getExonPhaseList(\@exonsListCopy);
		my @phaseList = @$ref2PhaseList;
		@phaseList = reverse(@phaseList) if ($strand_input eq "-");	## Now reverse this back so that the phases are in accordance with the different exons as given in the input file
		
		open(FO,">$mafFile") || die "Error !! Cannot write to mafFile '$mafFile'\n";
		my $noe = scalar(@exonsList);	
		my $exonMissingFlag = "";

		for (my $k = 0; $k < scalar(@exonsList); $k++)
		{
			## Get tag for the exon here:
			my $tag = "M"; ## by default every exon is a middle exon
			if (scalar(@exonsList) == 1){
				$tag = "S";
			}else{
				if ($strand_input eq "-"){
					$tag = "L" if ($k == 0);
					$tag = "F" if ($k == scalar(@exonsList) - 1);
				}else{
					$tag = "F" if ($k == 0);
					$tag = "L" if ($k == scalar(@exonsList) - 1);
				}
			}

			my ($exonStart,$exonStop) = (split /-/,$exonsList[$k])[0,1];
			my ($phase5Prime,$phase3Prime) = (split /-/,$phaseList[$k])[0,1];
			my $lengthRefExpected = 40 + abs($exonStart - $exonStop);
		
			my $key = "$RefSpecies#$accession#$exonStart#$exonStop#$chr#$strand_input#$phase5Prime#$phase3Prime#$tag";
			my $maf = readBDB($index,$key);

			my $mafCopy = $maf;
			$mafCopy =~s/\s+$//;
	
			if ($maf eq "" || $mafCopy eq "##maf version=1 scoring=mafExtract"){  ## For certain alignments, there is nothing for the reference exon either. Exclude such genes
				print "Exon '$exonsList[$k]' is missing. Cannot find the key '$key' in the BDB file '$index'\n" if ($verbose);
				$exonMissingFlag = "T";
				last;
			}

			## Also check if the length of the reference exon is not the exonLength+40 bp (since I ask for a 20bp flanking sequence on either side)
			my $lengthRef = 0;
			my @mafLines = split(/\n/,$maf)	;
			foreach my $mLine(@mafLines)
			{
				next if ($mLine !~/^s/);
				my ($species,$dm1,$dm2,$length,$dm3,$dm4,$dm5) = getSLineData($mLine);
				$lengthRef = $lengthRef + $length if ($species eq $RefSpecies);
			}
			
			if ($lengthRef != $lengthRefExpected) {
				print "$maf\n";
				print "Something is wrong with the exon alignment for '$exonsList[$k]'\n" if ($verbose);
				$exonMissingFlag = "T";
				last;
			}
			print FO "$maf\n\n";
		}
		close FO;
		
		if ($exonMissingFlag eq "T"){
			print FOL "Skipping $line\t Not all exons realigned for this gene\n";
			next;
		}
	
		my $statusMSS = system("mafSpeciesSubset $mafFile speciesList=$includeSpeciesString species.lst=NULL $mafFileOut");		## statusMSS --> statusMafSpeciesSubset
		safeDie(\@tmpFilesList,"Error running mafSpeciesSubset for '$line' bdb maf file\n") if ($statusMSS != 0);
		`mv $mafFileOut $mafFile`;
	
		my $Ref2SpeciesList = GetSpeciesList($mafFile);		## This gives me a list of species present in the maf.
		my %speciesPresentInMaf = map{ $_ => 1} @$Ref2SpeciesList;
		## Write to logFile for every species that was present in the includeSpeciesString but could not be found in the mafFile
		foreach my $species(keys(%includeSpeciesHash))
		{
			if (! exists $speciesPresentInMaf{$species}){
				print FOL "$accession\t$geneName\t$species\tIgnored. Nothing in the maf for this gene\n";
				delete $includeSpeciesHash{$species}; ### And delete them from the includeSpeciesHash
			}
		}
		### Filter out those species where the gene does not come from the same chromosome/strand
		if ($sameSSFlag){
			my ($ref2ExcludeSpecies) = checkForSameScaffoldStrand($mafFile,\%includeSpeciesHash);
			{
				my %excludeSpeciesList = %$ref2ExcludeSpecies;
				foreach my $species(keys(%excludeSpeciesList))
				{
					print FOL "$accession\t$geneName\t$species\t$excludeSpeciesList{$species}\n";
					delete $includeSpeciesHash{$species}; ### And delete them from the includeSpeciesHash
				}
			}
			my $includeSpeciesListUpdated = join(",",keys(%includeSpeciesHash));
		
			my $statusMSS = system("mafSpeciesSubset $mafFile speciesList=$includeSpeciesListUpdated species.lst=NULL $mafFileOut");		## statusMSS --> statusMafSpeciesSubset
			`mv $mafFileOut $mafFile`;
		}
	}
	
	else{
		
		my @exonsList = split(",",$exons);
		open(FOB,">$tmpBed");		## Create a bed file with ONLY coding exons
		foreach my $exon(@exonsList)
		{
			my @tmpExon = split(/-/,$exon);
			print FOB "$chr\t$tmpExon[0]\t$tmpExon[1]\n";
		}
		close FOB;
		## And then run mafExtract + mafSpeciesSubset on this bed file
		my $call = "set -o pipefail; mafExtract $index -regionList=$tmpBed stdout|mafSpeciesSubset stdin speciesList=$includeSpeciesString species.lst=NULL $tmpMaf";
		system("$call")	== 0 || safeDie(\@tmpFilesList,"Error running '$call'\n");
		
		my $Ref2SpeciesList = GetSpeciesList($tmpMaf);		## This gives me a list of species present in the maf.
		my %speciesPresentInMaf = map{ $_ => 1} @$Ref2SpeciesList;
		## Write to logFile for every species that was present in the includeSpeciesString but could not be found in the mafFile
		foreach my $species(keys(%includeSpeciesHash))
		{
			if (! exists $speciesPresentInMaf{$species}){
				print FOL "$accession\t$geneName\t$species\tIgnored. Nothing in the maf for this gene\n";
				delete $includeSpeciesHash{$species}; ### And delete them from the includeSpeciesHash
			}
		}
		### Filter out those species where the gene does not come from the same chromosome/strand
		if ($sameSSFlag){
			my ($ref2ExcludeSpecies) = checkForSameScaffoldStrand($tmpMaf,\%includeSpeciesHash);
			{
				my %excludeSpeciesList = %$ref2ExcludeSpecies;
				foreach my $species(keys(%excludeSpeciesList))
				{
					print FOL "$accession\t$geneName\t$species\t$excludeSpeciesList{$species}\n";
					delete $includeSpeciesHash{$species}; ### And delete them from the includeSpeciesHash
				}
			}
		}
		my $includeSpeciesStringUpdated = join(",",keys(%includeSpeciesHash));
		
		## Now run mafExtract again with a different bedFile, here we use a flanking region of 20bp on either side of every exon, we also filter this maf for includeSpeciesList
		my $callFinal = "set -o pipefail; mafExtract $index stdout -regionList=$bedFile|mafSpeciesSubset stdin -speciesList=$includeSpeciesStringUpdated species.lst=NULL $mafFile"; #### Run mafExtract to get a maf file from the input bed file.
		system("$callFinal") == 0 || safeDie(\@tmpFilesList,"Error running '$call'\n");
		
		print "mafExtract run was successful for $accession\n\n" if ($verbose);
	}
	
	## Now I have a maf and it only contains the species that I am interested in. All the useless stuff has gone away: Time to generate the tabular files:
	my $out_dir = "$data_path/$accession";
	`mkdir -p $out_dir`;
	
	my %genome_cds_nature_hash = ();	## This hash stores the position based information for every genomic coordinate in the reference species, e.g. xxx is coding_exon_1, xxy 1_A_1 etc.
	my %codon_index_hash = ();			## This hash stores the codon index for every genomic coordinate in the reference species e.g. xxx is C1, xxy is C2 etc.
	open(FI,"$exon_ct_file") || die "Error opening Genomic coordinates annotation file\n";
	while (my $line=<FI>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
		
		$genome_cds_nature_hash{$tmp[1]} = $tmp[0];
		$codon_index_hash{$tmp[1]} = $tmp[2];
	}
	close FI;

	my @species_list = keys(%includeSpeciesHash);	### A list of species that I want in my maf
	### Before you do the entire processing... check if the accession is already done and dusted. In that case, do not do anything.
	
	my $flagInComplete = "T";
	my $tabFileAllSpecies = "$out_dir.txt";

	$flagInComplete = checkForCompleteness(\@species_list,$tabFileAllSpecies) if (-e $tabFileAllSpecies);
	
	if ($flagInComplete eq "T"){	## If indeed the tabFileDir is incomplete, then do this operation
			
		GetHashFromMafAndPrintIt($mafFile,$RefSpecies,$out_dir,$strand_input,\%genome_cds_nature_hash,\%codon_index_hash);
		addReadingFrame_TabularFiles($out_dir,\@species_list,$RefSpecies);
		formatFiles($out_dir,\@species_list);
		
		## Combine all individual files to make one composite tabular file
		open(FO,">$data_path/$accession.txt") || die "Error writing to composite tabular file\n";		
		foreach my $species(@species_list)
		{
			open(FIS,"$out_dir/$species.txt") || die "Error opening species file for '$species' in $out_dir\n";
			while(my $l = <FIS>)
			{
				print FO "$l";
			}
			close FIS;
		}
		close FO;
		`rm -rf $out_dir`;
	}else{
		print "MESSAGE: Nothing to be done for this $accession --> All the tabular files are already present\n";
	}
	
	print "####\n\nAll operations with $accession are finished\n####\n\n" if ($verbose);
	`rm $exon_ct_file $bedFile`;  ## These need to be deleted because they are created inside a temp directory.
}
close FIL;
close FOL;

my $tmpFilesListString = join(" ",@tmpFilesList);	## As a conservative measure, delete all temp files though there should not be any left by now
`rm -rf $tmpFilesListString`;

#### Sub-routines for maf-processing:

sub checkForCompleteness
{
	my ($ref2SpeciesList,$tabFile) = @_;
	my @species_list = @$ref2SpeciesList;
	
	my $flagInComplete = "F";

	my $ct = 0;
	open(FITA,$tabFile) || die "Error opening tabFile '$tabFile' in checkForCompleteness function\n";
	while (my $line = <FITA>)
	{
		$line =~s/\s+$//;
		$ct++ if ($line =~/FileDone/);
	}
	close FITA;
	$flagInComplete = "T" if ($ct != scalar(@species_list));
	return $flagInComplete;
}

sub addReadingFrame_TabularFiles
{
	my($outDir,$ref2SpeciesList,$ref) = @_;
	my @speciesList = @$ref2SpeciesList;
	
	foreach my $species(@speciesList)
	{
		my $refFile	= "$outDir/$ref.txt";
		my $speciesFile = "$outDir/$species.txt";

		open(FIS,$speciesFile) || die "Error opening speciesFile for '$species' in addReadingFrame_TabularFiles\n";
		my @fileArray = <FIS>;
		close FIS;
		
		my $i = 1;
		my $noe = get_noe($refFile);
		my %readingFrame = my %readingFrameEnd = ();
		

		my $ref2ExonBoundaryHash = getStartStopPos($speciesFile,"NA");  ## I want the entire hash back and not just the boundaries for one particular exon
		my %exonBoundaryHash = %$ref2ExonBoundaryHash;  ## Gives me the boundary for each exon in the tabFile
		
		while($i <= $noe)  ## Find the reading frame for each exon, plus the exon boundaries in the tabFile
		{
			my ($rf,$rfLast,$flagUpstreamInsertion,$flagDFBI) = getReadingFrame($refFile,$i,$speciesFile);  ## Get reading frame for each exon, the start of the exon and the end of the exon
			
			$readingFrame{$i} = $rf;
			$readingFrameEnd{$i} = $rfLast;
			
			if ($flagUpstreamInsertion eq "T" || $flagDFBI eq "T"){ ## in that case, the reading frame is undetermined and I have to look at the upstream exon
				if ($i == 1){
					$readingFrame{$i} = 1;
				}else{
					my $rfUpstreamExon = $readingFrameEnd{$i - 1};
				
					my $rfCurrent = $rfUpstreamExon +1;
					$rfCurrent = $rfCurrent%3;
					$rfCurrent = 1 if ($rfCurrent == 0);
					$readingFrame{$i} = $rfCurrent;
				}
			}
			$i++;
		}
		
		my %nonExonBoundaryHash = ();
		foreach my $keys(sort {$a <=> $b} keys(%exonBoundaryHash))   ## Here I get the boundaries for the non-Exonic regions in the tabFile i.e the first 20 lines, the intronic regions between two exons
		{
			my $nonExonBoundary = "";
			if ($keys == 1){
				my($startC,$stopC) = (split /-/,$exonBoundaryHash{$keys})[0,1];    ## startC is startCurrent, stopC is stopCurrent
				$startC--;
				$nonExonBoundary = "0-$startC";
			}else{
				my($startC,$stopC) = (split /-/,$exonBoundaryHash{$keys-1})[0,1];    ## startC is startCurrent, stopC is stopCurrent
				my($startN,$stopN) = (split /-/,$exonBoundaryHash{$keys})[0,1];  ## startN is startNext
				$stopC++;
				$startN--;
				$nonExonBoundary = "$stopC-$startN";
			}
			$nonExonBoundaryHash{$keys} = $nonExonBoundary;
		}
	
		## Last 20 lines
		my($dm,$startL20) = (split /-/,$exonBoundaryHash{$noe})[0,1];
		$startL20++;
		my $stopL20 = `wc -l < $speciesFile`; chomp $stopL20;
		$stopL20--;
		
		open(FO,">$speciesFile");
		foreach my $exonNumber(sort {$a <=> $b} keys(%exonBoundaryHash))
		{
			my $nonExonicRange = $nonExonBoundaryHash{$exonNumber};
			my($startNE,$stopNE) = (split /-/,$nonExonicRange)[0,1];
	
			my @chunkNonExonic = @fileArray[$startNE..$stopNE];
			foreach my $line(@chunkNonExonic)
			{
				chomp $line;
				print FO "$line\t-\n";  ## No codon position for nonExonic regions
			}	

			my $exonicRange = $exonBoundaryHash{$exonNumber};
			my($startE,$stopE) = (split /-/,$exonicRange)[0,1];
		
			my $rf = $readingFrame{$exonNumber};
			## Enforce the ancestral/reference reading frame on the query species if the exonic sequence in the query begins with a deletion
			my @chunkExonic = @fileArray[$startE..$stopE];
			my $deletionsAtStart = 0;
			foreach my $line(@chunkExonic) ## Now indicate reading frame/codon positions for every base in the exonic region
			{
				chomp $line;
				my $base = (split/\t/,$line)[11];
				if ($base eq "-"){
					$deletionsAtStart++;
				}else{
					last;
				}
			}
			$rf = $rf + $deletionsAtStart;  ## Updated
			
			foreach my $line(@chunkExonic) ## Now indicate reading frame/codon positions for every base in the exonic region
			{
				chomp $line;
				my $base = (split /\t/,$line)[11];

				my $codonPosPrint = "";
				if ($base ne "-"){		## Treat every character as a real nucleotide unless it is a deletion ("-").
					if ($base eq "?"){
						$codonPosPrint = "-"; ## Do not "print a real codon position in case of missing data even though I consider it as some sequence there.
					}else{
						$codonPosPrint = ($rf % 3);
						$codonPosPrint = 3 if ($codonPosPrint == 0);
						$codonPosPrint = "C".$codonPosPrint;
					}
					$rf++;
				}else{
					$codonPosPrint = "-";
				}
				print FO "$line\t$codonPosPrint\n";
			}
		}
		
		## Print last 20 lines
		my @chunkNonExonicLast = @fileArray[$startL20..$stopL20];
		foreach my $line(@chunkNonExonicLast)
		{
			chomp $line;
			print FO "$line\t-\n";	
		}
		close FO; ## Thats it, done
	}
}

sub formatFiles
{
	my($outDir,$ref2SpeciesList) = @_;
	my @speciesList =  @$ref2SpeciesList;
	
	foreach my $species(@speciesList)
	{
		open(FIS,"$outDir/$species.txt") || die "Error opening tabFile for species '$species' in formatFiles function\n";
		my @fileArray = <FIS>;
		close FIS;
		
		open(FO,">$outDir/$species.txt") || die "Error writing to tabFile for species '$species' in formatFiles function\n";
		#print FO "FileBegin-$species\n";
		foreach my $line(@fileArray)
		{
			$line =~s/\s+$//;
			my($lNumber,$ref,$nature,$rfRef,$strandRef,$cdsRef,$ntRef,$chrRef,$query,$strandQuery,$cdsQuery,$ntQuery,$chrQuery,$rfQuery) = (split /\t/,$line);

			## Convert $ntRef and $ntQuery to upper case, just in case there are in lower case
			$ntRef 	 = uc($ntRef);
			$ntQuery = uc($ntQuery);			
			print FO "$lNumber\t$ref\t$cdsRef\t$chrRef\t$strandRef\t$nature\t$rfRef\t$ntRef\tX\t$query\t$cdsQuery\t$chrQuery\t$strandQuery\t$ntQuery\tX\t$rfQuery\n";  ## No protein translation, just put X there
		}
		#print FO "FileDone-$species\n";
		print FO "FileDone\n";
		close FO;
	}
}

															###############################
sub getReadingFrame
{
	my ($refFile,$exon,$speciesFile) = @_;
	
	my $flagUpstreamInsertion = "F";
	my $exonicTerm = "Coding_exon\_$exon";  ### This gives me the reading frame for the exon
	
	my $readingFrame = my $readingFrameLast = "";
	open(FIR,$refFile) || die "Error opening reference file '$refFile' in get_reading_frame function\n";
	while (my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my($nature,$codonPosition,$base) = (split /\t/,$line)[2,3,11];
		if ($nature eq $exonicTerm && $base ne "-"){
			$readingFrame = $codonPosition;
			$readingFrame =~s/C//;
			last;
		}
	}
	close FIR;
	
	open(FIR,$refFile) || die "Error opening reference file '$refFile' in get_reading_frame function\n";
	while (my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my($nature,$codonPosition,$base) = (split /\t/,$line)[2,3,11];
		if ($nature eq $exonicTerm && $base ne "-"){
			$readingFrameLast = $codonPosition;
			$readingFrameLast =~s/C//;
		}
	}
	close FIR;
	
	## Now check if there is an insertion upstream of the start of the exon. This needs to be dealt with because there is no reading frame that is defined for this position. And of course, this should not be checked for the first exon
	## Firstly, there is no line upstream of the first coding position of the first exon. And even if there was one, that would not correspond to the Acceptor positions
	
	if ($exon != 1)
	{
		my $ct = 0;
		open(FIR,$refFile) || die "Error opening reference file '$refFile' in get_reading_frame function\n";
		while (my $line = <FIR>)
		{
			$line =~s/\s+$//;
			my $nature = (split /\t/,$line)[2];
			last if ($nature eq $exonicTerm);
			$ct++;
		}
		close FIR;
		
		open(FIR,$refFile) || die "Error opening reference file '$refFile' in get_reading_frame function\n";
		my @fileArray = <FIR>;
		close FIR;
		
		my $lineOfInterest = $fileArray[$ct-1];
		my $thirdFieldLOI = (split /\t/,$lineOfInterest)[2];
		$flagUpstreamInsertion = "T" if ($thirdFieldLOI !~/_A_1/); ## If no, set the flag to T.		
	}
	
	my $flagDFBI = "F";   ## Deletion followed by Insertion
	
	if ($readingFrame eq ""){ 
		my $exonicChunk = get_exonic_chunk($exon,$refFile,$speciesFile);
		my @exonChunk = split(/\n/,$exonicChunk);
		my $line = $exonChunk[0];
		my $firstBaseQuery = (split /\t/,$line)[11];
		
		if ($firstBaseQuery eq "-"){   ## the first base in the query is indeed a deletion
			foreach my $line(@exonChunk)
			{
				my($ntRef,$ntQuery) = (split /\t/,$line)[6,11];
				if ($ntRef eq "-" && $ntQuery ne "-")
				{
					$flagDFBI = "T";
					last;
				}
			}
		}
	}
	
	return ($readingFrame,$readingFrameLast,$flagUpstreamInsertion,$flagDFBI);
}

sub checkForSameScaffoldStrand
{
	my($tmpMaf,$ref2IncludeSpeciesHash) = @_;
	
	my %includeSpeciesHash = %$ref2IncludeSpeciesHash; 
	my %chrHash = my %strandHash = my %speciesExclude = ();
	
	open(FI,"$tmpMaf");		## Now read the maf, and collect information about the strand and chromosome for every species in the maf (from the "s" lines)
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		
		if ($line =~/^s/){
			my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line);
			$chrHash{$species}{$chr} 	= "T";
			$strandHash{$species}{$strand}	= "T";
		}
	}
	close FI;
	
	foreach my $species(keys(%includeSpeciesHash))
	{
		my @chrArray = my @strandArray = ();
			
		foreach my $chr (keys %{$chrHash{$species}})
		{
			push(@chrArray,$chr);
		}
			
		foreach my $strand (keys %{$strandHash{$species}})
		{
			push(@strandArray,$strand);
		}
			
		my $chrList = join(",",@chrArray);
			
		if (scalar(@chrArray) > 1 || scalar(@strandArray) > 1){
			my @logMessage = ();
			my $m = 0;
				
			push(@logMessage,"Ignored. The gene in the maf comes from more than one chromosome (i.e. $chrList)") if (scalar(@chrArray) > 1);
			push(@logMessage,"Ignored. The gene in the maf comes from two strands") if (scalar(@strandArray) > 1);
				
			my $message = join(", ",@logMessage);
			$speciesExclude{$species} = $message;
		}
	}
	return (\%speciesExclude);
}

sub safeDie
{
	my ($ref2List,$message);
	
	my @tempFilesList = @$ref2List;
	my $tempList = join(" ",@tempFilesList);
	`rm -rf $tempList`;
	die $message;
}

sub getExonPhaseList		## This subroutine gives the phases (5' and 3') of each exon.
{
	my $ref2ExonsList = shift;
	my @exonsList =  @$ref2ExonsList;
	
	my $i = 0;
	my @exonLength = ();
	my $lengthCum = 0;
	
	my @exonPhaseList = ();
	
	foreach my $exon(@exonsList)
	{
		my ($start,$stop) = (split /-/,$exon)[0,1];
		my $length = abs($start-$stop);
		$lengthCum = $lengthCum + $length;
		
		$exonLength[$i] = $lengthCum;
		$i++;
	}
	
	for (my $i = 0;$i < scalar(@exonLength); $i++)
	{
		my $exonLengthCum = $exonLength[$i];
		
		my $phase5Prime = my $phase3Prime = "";
		
		if ($i != 0){
			$phase5Prime = 3 - ($exonLength[$i-1]%3);
			$phase3Prime = ($exonLength[$i]%3);
		}else{
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


sub get_exonic_chunk
{
	my ($exonIndex,$refFile,$speciesFile) = @_;				
	
	if (-e $speciesFile){	
		open(FIR,$speciesFile) || die "Error opening species file '$speciesFile' in get_exonic_chunk function\n";
		my @fileArray = <FIR>;
		close FIR;  
		
		my $exonic_term = "Coding_exon\_$exonIndex";
		my $tmp_result = `grep -w "$exonic_term" $refFile`;
		die "Improper exon index '$exonIndex'. This exon is not present in the file $refFile \n" if ($tmp_result eq "");
		
		my $term1 = "$exonIndex\_A_1";
		my $term2 = "$exonIndex\_D_1";
		
		## Get the number of exons
		my $noe = get_noe($refFile);
		$term1 = "Coding_exon_1" if ($exonIndex == 1);
		$term2 = "Coding_exon\_$noe" if ($exonIndex == $noe);
		
		## Get start and stop positions:
		my ($startPt,$stopPt) = getStartStopPos($speciesFile,$exonIndex);
		
		my @exonicSeqArray = @fileArray[$startPt..$stopPt];
		my $exonicSeq = join("",@exonicSeqArray);
		
		return $exonicSeq;
	}else{
		die "The gene is missing since '$speciesFile' is not present\n";
	}
}

sub get_noe
{
	my $refFile = shift;

	my $noe = "";
	open(FIR,"$refFile") || die "Error opening input file '$refFile' in get_noe function'\n";
	while(my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my $nature = (split /\t/,$line)[2];
		$noe = $nature if ($nature =~/Coding_exon/);
	}
	close FIR;
	$noe =~s/Coding_exon_//;
	
	return $noe;
}

sub getStartStopPos  ## Given a tabular file, this function returns a hash (when the second argument is "NA") where the key is the exon number while the value is a string of the type "start-stop"
{					 ## start indicates the line where the exon starts in the tabFile while stop indicates the line where the exon ends in the tabFile
					 ## If the second argument is an exon number, the function returns the start/stop for this exon number:
	my ($file,$exonNumber) = @_;
	
	my $noe = get_noe($file);
	die "'$exonNumber' greater than the total number of exons '$noe' in '$file'\n" if ($exonNumber ne "NA" && $exonNumber > $noe);
	
	my @termsList = my @termsListPair =  ();
	push(@termsList,"Coding_exon_1","1_D_1");
	if ($noe > 1){	
		my $start = 2;
		while ($start <= $noe)
		{
			my $term1 = "$start\_A_1";
			my $term2 = "$start\_D_1";
			push(@termsList,$term1,$term2);
			
			$start++;
		}
	}
	
	my $term2 = "Coding_exon\_$noe";
	push(@termsList,$term2);
		
	@termsList = grep{$_ ne "$noe\_D_1"} @termsList;  ## throw out LastExon_D_1
	
	## Now make pairs
	for(my $i = 0; $i < scalar(@termsList)-1; $i = $i + 2)
	{
		my $p1 = $termsList[$i];
		my $p2 = $termsList[$i+1];
		push(@termsListPair,"$p1-$p2");
	}	
	
	## Get the line start and line stops:
	my $ct = 0;
	my %hashAllTerms = ();
	
	open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
	while (my $line = <FIR>)
	{
		my $nature = (split /\t/,$line)[2];
		$hashAllTerms{$nature} = $ct;
		
		$ct++;
	}
	close FIR;
	
	##  Exception for first exon
	my $ctStart = 0;
	open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
	while (my $line = <FIR>)
	{
		my $nature = (split /\t/,$line)[2];
		last if ($nature eq "Coding_exon_1");
		$ctStart++;
	}
	close FIR;
	
	my %exonBoundariesHash = ();
	my $i = 1;
	foreach my $pair(@termsListPair)
	{
		my($start,$stop) = (split /-/,$pair);
		
		my $begin = "";
		if ($i == 1){
				$begin = $ctStart;
		}else{
				$begin = $hashAllTerms{$start};
		}
		
		my $end   = $hashAllTerms{$stop};
		$begin++ if ($start =~/_A_/);
		$end-- if ($stop =~/_D_/);
		
		$exonBoundariesHash{$i} = "$begin-$end";
		$i++;
	}
	
	if ($exonNumber eq "NA"){
			return \%exonBoundariesHash;
	}else{
		my($start,$stop) = (split /-/, $exonBoundariesHash{$exonNumber})[0,1];
		return($start,$stop);
	}
}
