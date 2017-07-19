#!/sw/bin/perl
## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;

my $usage = "perl $0 ReferenceSpecies QuerySpecies geneChainFilteringDir outputDirectoryWhite outputDirectoryBlack\n";
my $purpose = "$0 is the preprocessing script for mafFilterForWhiteListedChains.pl. Creates bedFiles from whiteListed chains which are used by the main script\n";
die $usage.$purpose if (scalar(@ARGV) != 5);

my $ref 		  = $ARGV[0];
my $species 		  = $ARGV[1];
my $chainsDir	          = $ARGV[2];
my $outDirWhite		  = $ARGV[3];
my $outDirBlack		  = $ARGV[4];

`mkdir -p $outDirWhite $outDirBlack`;
my $chainFileWhite = "$chainsDir/ChainPseudoGeneParalogFiltering/$ref.$species.whiteListed.chain.gz";
die "WhiteListed ChainsFile '$chainFileWhite' does not exist. Exiting\n" if (! -e $chainFileWhite);

## collate the blackListed chains
my $bl1 = "$chainsDir/GeneChainFiltering/$ref.$species.nonGenespanning.chain.gz";
my $bl2 = "$chainsDir/ChainPseudoGeneParalogFiltering/$ref.$species.blacklistedPseudogene.chain.gz";

die "The nonGeneSpanning blackListed chains file does not exist\n" if (! -e $bl1);
die "The pseudogene/paralog like blackListed chains file does not exist\n" if (! -e $bl2);

my $tmpDir = `mktemp -d`; chomp $tmpDir;
my $call = "zcat $bl1 $bl2 |  gzip -c > $tmpDir/$ref.$species.blackListedAll.chain.gz";
system($call) == 0 || die "Error running '$call'\n";

createBedFiles($chainFileWhite,$outDirWhite);  
createBedFiles("$tmpDir/$ref.$species.blackListedAll.chain.gz",$outDirBlack);
`rm -rf $tmpDir`;

sub createBedFiles
{
	my ($chainFile,$outDir) = @_;
	
	### Run chain to alignedBed on the chainFile first
	my $tmpBed = `mktemp /dev/shm/XXXXX.bed`; chomp $tmpBed;
	my $call = "ChainToAlignBed.perl $chainFile $tmpBed -mergeTargetBlocks -mergeTargetBlocks_Query";
	system($call) == 0 || die "Error running '$call'\n";

	## Now sort the file:
	my $callBedSort = "bedSort $tmpBed $tmpBed";
	system($callBedSort) == 0 || die "Error running bedSort\n";

	`mkdir -p $outDir/$species`;
	## Now split the sorted bedFile into chrom-split files:
	my $callBedSplit = "bedSplitOnChromUnlimited $tmpBed $outDir/$species/";
	system($callBedSplit) == 0 || die "Error running bedSplit\n";

	## Sort every bedFile, it helps later
	my $listBedFiles = `ls $outDir/$species`;
	my @listFiles = split(/\n/,$listBedFiles);
	foreach my $f(@listFiles)
	{
		my $bedFileIn = "$outDir/$species/$f";
		my $callBedSort = "bedSort $bedFileIn $bedFileIn";
		system($callBedSort) == 0 || die "Error running '$callBedSort'\n";
	}

	`rm -rf $tmpBed`; ## Clean up
}
