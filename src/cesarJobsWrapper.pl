#!/sw/bin/perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my $parameters = "";
my $masking = 0;
GetOptions ("parameters:s" => \$parameters, "masking" => \$masking);

my $usage = "USAGE: $0 bbIndex referenceSpecies speciesList outputDir intronTableDir clade splitGenesDir logDir jobsListName program[CESAR|CESAR2] parameters[Optional]\n";
my $purpose = "PURPOSE: $0 is a wrapper around RealignmentPipeline_Exons.pl. It creates the jobs file that need to be submitted via para\n";
die $usage.$purpose if (scalar(@ARGV) < 10);

my $index 		= $ARGV[0];
my $ref   		= $ARGV[1];
my $queryList 		= $ARGV[2];
my $outDir	 	= $ARGV[3];
my $intronTable 	= $ARGV[4];
my $clade		= $ARGV[5];
my $dir 		= $ARGV[6];
my $logDir 		= $ARGV[7];
my $jobsPrefix		= $ARGV[8];
my $program             = $ARGV[9];

## check paramaters
die "bbIndex '$index' does not exist\n" if (! -e $index);
die "SpeciesList is not specified\n" if (! -e $queryList);
die "outputDirectory not specified\n" if ($outDir eq "");
die "intronTablesDir '$intronTable' does not exist\n" if (! -d $intronTable);
die "logDir is not specified\n" if ($logDir eq "");
die "splitExonsList directory does not exist\n" if(! -d $dir);
die "jobsPrefix is not specified\n" if ($jobsPrefix eq "");
die "Clade not specified\n" if ($clade eq "");
die "Reference species not specified\n" if ($ref eq "");
die "Program '$program' is not recognized\n" if ($program ne "CESAR" && $program ne "CESAR2");

`mkdir -p $logDir $outDir`;
## All OK, time for business
open(FO1,">$jobsPrefix.short") || die "Error writing to shortJobs file, the prefix is '$jobsPrefix'\n";
open(FO2,">$jobsPrefix.medium") || die "Error writing to mediumJobs file, the prefix is '$jobsPrefix'\n";;
open(FO3,">$jobsPrefix.long") || die "Error writing to longJobs file, the prefix is '$jobsPrefix'\n";;

my $ls = `ls $dir`;
my @list = split(/\n/,$ls);

my $ct = 1;
my $nShort = my $nMed = my $nLong = 0;
foreach my $file(@list)
{
	my $fileIn = "$dir/$file";
	my $outBDB = "$outDir/realigned.$ct.BDB";
	my $logFile = "$logDir/$ct.log";	
	
	my $line = "RealignmentPipeline_Exons.pl -in $fileIn -bbIndex $index -ref $ref -queryList $queryList -outBDB $outBDB -clade human  -deletedIntronsBDB $intronTable -intronTableDir $intronTable  -log $logFile -program $program";
	my $lineFormat = formatLine($line);
	
	if ($file =~/short/){
		print FO1 "$lineFormat\n";
		$nShort++;
	} elsif ($file =~/medium/){
		print FO2 "$lineFormat\n";
		$nMed++;
	}elsif ($file =~/long/){
		print FO3 "$lineFormat\n";
		$nLong++;
	}else{
		die "file type not recognised for '$file'\n";
	}
	$ct++;
}
close FO1;
close FO2;
close FO3;

deleteIfEmpty("$jobsPrefix.short") if ($nShort == 0);  ## Never gonna happen
deleteIfEmpty("$jobsPrefix.medium") if ($nMed == 0);   ## Unlikely to happen
deleteIfEmpty("$jobsPrefix.long") if ($nLong == 0);    ## May happen	

print "ALL DONE\n $nShort jobs for short queue, $nMed jobs for medium queue and $nLong jobs for long queue produced in '$jobsPrefix' jobs file\n";

sub formatLine
{
	my $line = shift;
	$line = $line." -parameters $parameters" if ($parameters ne "");
	$line = $line." -masking"if ($masking);
	return $line;
}

sub deleteIfEmpty
{
	my $file = shift;
	print "'$file' deleted because there are no jobs in this file\n";
}
