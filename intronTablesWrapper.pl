#!/sw/bin/perl
## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

use strict;
use warnings;
use BerkeleyDB;

my $usage = "perl $0 speciesList allTranscriptsList (the full path needs to be specified for this file) ReferenceSpecies bigBedIndex clade (human|mouse etc) jobsPrefix program[CESAR|CESAR2]\n";
my $purpose = "$0 is a wrapper around getIntronTables.pl\n";
die $usage.$purpose if (scalar(@ARGV != 7));

my $speciesFile 	    = $ARGV[0];
my $allTranscriptsFile  = $ARGV[1]; ## The absolute path
my $ref			        = $ARGV[2];
my $index		        = $ARGV[3];
my $clade		        = $ARGV[4];
my $jobsPrefix		    = $ARGV[5]; ## jobsPrefix file
my $program             = $ARGV[6];

## Check parameters
die "Species file '$speciesFile' does not exist\n" if (! -e $speciesFile);
die "The allTranscripts file '$allTranscriptsFile' does not exist\n" if (! -e $allTranscriptsFile);
die "The bigBed index '$index' does not exist\n" if (! -e $index);
die "The value '$program' is not valid for the program\n" if ($program ne "CESAR" && $program ne "CESAR2");

## All paramaters OK
my $exonsUnionFile     = "$allTranscriptsFile\_exonUnion";
my $splitGenesDir      = "splitGenesDir";
`mkdir -p $splitGenesDir`;

## Step 0, run the splitGenePredForReAlignment.pl tool to split the allTranscriptsGenePred to smaller files such that each file is a job
my $callSGP = "splitGenePredForReAlignment.pl $allTranscriptsFile $exonsUnionFile $splitGenesDir 4000 50";
system($callSGP) == 0 || die "Error running '$callSGP'. Aborting at the first step\n";

## Check if intronTables and intronTables_Split directory exist. If yes, then die with an error message:
die "intronTables and|or intronTables_Split already directory exist, delete the existing directories\n" if (-d "intronTables" || -d "intronTables_Split");

## Step1, split the allTranscriptsGenePred file
my $pwd = `pwd`; chomp $pwd;
my $tempDir = "splitGenesList";
`mkdir -p $tempDir`;

chdir $tempDir;
my $callSplit = "splitFile $allTranscriptsFile 2000 file";  ### 2000 transcripts per job is a good number for the jobs to finish in short queue
system($callSplit) == 0 || die "Error running '$callSplit'\n";
chdir $pwd;

my $lsFiles = `ls $tempDir`;
my @listFiles = split(/\n/,$lsFiles);

`mkdir -p intronTables intronTables_Split exonGroupsDir exonGroupsRealigned`;
open(FO,">jobsIntronTables");
foreach my $file(@listFiles) {
	print FO "getIntronTables.pl $tempDir/$file $ref $speciesFile $index intronTables_Split\n";
}

## Step 2, push the jobsIntronTables jobs
my $paraCall = "para make jobsIntronTables jobsIntronTables";
system($paraCall) == 0 || die "Error running jobsIntronTables, '$paraCall' failed\n";

print "#####\n\ngetIntronTables.pl run successfull. Now merging the output of different split files\n#####\n\n";

## Step 3, merge the individual outputs of jobsIntronTables
open(FIS,$speciesFile) || die "Error opening speciesList file '$speciesFile'\n";
open(FOS,">exonGroups") || die "Error writing to exonGroups file\n";
open(FOR,">realignmentExonGroups") || die "Error writing to realignmentExonGroups file\n";

while (my $species = <FIS>) {
	$species =~s/\s+$//;
	my @list5P = my @list3P = my @list5PDel = my @list3PDel = my @listExclude = ();

	print "Here you are, the species is '$species'\n";
	
	foreach my $file(@listFiles) {
		my $f1 = "intronTables_Split/".$file.".".$species.".fivePrimeTable"; push(@list5P,$f1);
		my $f2 = "intronTables_Split/".$file.".".$species.".fivePrimeTableDel"; push(@list5PDel,$f2);
		my $f3 = "intronTables_Split/".$file.".".$species.".threePrimeTable"; push(@list3P,$f3);
		my $f4 = "intronTables_Split/".$file.".".$species.".threePrimeTableDel"; push(@list3PDel,$f4);
		my $f5 = "intronTables_Split/".$file.".".$species.".exclude"; push(@listExclude,$f5);
	}
	
	mergeFunction(\@list5P,"intronTables/$species.fivePrimeTable");
	mergeFunction(\@list5PDel,"intronTables/$species.fivePrimeTableDel");
	mergeFunction(\@list3P,"intronTables/$species.threePrimeTable");
	mergeFunction(\@list3PDel,"intronTables/$species.threePrimeTableDel");
	mergeFunction(\@listExclude,"intronTables/$species.ignoreExons");
	
	## Convert the tables to BDB files
	fileToBDB("intronTables/$species.fivePrimeTable","intronTables/$species.5Prime.BDB");
	fileToBDB("intronTables/$species.threePrimeTable","intronTables/$species.3Prime.BDB");
	
	## Also remove the tables, I only need the BDB files
	`rm intronTables/$species.fivePrimeTable intronTables/$species.threePrimeTable`;
	
	my $fivePrimeTableDel  = "intronTables/$species.fivePrimeTableDel";
	my $threePrimeTableDel = "intronTables/$species.threePrimeTableDel";
	
	my $n = 0;
	open(FIST,$fivePrimeTableDel) || die "Error opening 5PrimeDelTable  '$fivePrimeTableDel'\n";
	while (my $l = <FIST>) {
		$n++;
	}
	close FIST;
	
	if ($n > 0) {
		print FOS "groupExons.pl $allTranscriptsFile intronTables/$species.fivePrimeTableDel intronTables/$species.threePrimeTableDel exonGroupsDir/$species.groups $exonsUnionFile $species intronTables/$species.ignoreExons\n";
		print FOR "RealignmentPipeline_ExonGroups.pl -in exonGroupsDir/$species.groups -bbIndex $index -ref $ref -intronLookUpTables intronTables -outBDB intronTables/$species.realigned.BDB -clade $clade -program $program\n" ;
	} else {
		`rm $fivePrimeTableDel $threePrimeTableDel`;	
	}
}
close FOS;
close FOR;

## Step 4, group exons
`chmod +x exonGroups`;
system("exonGroups") == 0 || die "Error running exonGroups\n";

## Step 5, push Realignment_ExonGroups.pl jobs
if ($program eq "CESAR") { 
	my $realignmentCall = "para make realignmentExonGroups realignmentExonGroups -q medium";
	system($realignmentCall) == 0 || die "Error running '$realignmentCall'\n";
} else {
	`chmod +x realignmentExonGroups`;  ## For CESRA2, there is no need to submit jobs to cluster
	system("realignmentExonGroups") == 0 || die "Error running realignmentExonGroups\n";
}

## Remove unnecessary stuff
`rm -rf $tempDir jobsIntronTables exonGroups intronTables_Split`;

`mkdir -p realignedExons realignedExonsLog`;
## Step 6, final step --> create jobs file for exons that need to be submitted via para
my $createJobsFileCall = "cesarJobsWrapper.pl $index $ref $speciesFile realignedExons intronTables $clade $splitGenesDir realignedExonsLog $jobsPrefix $program";
system($createJobsFileCall) == 0 || die "Error running '$createJobsFileCall'\n";

## Done

### Sub-routines
sub mergeFunction {
	my ($ref2List,$outFile) = @_;
	my @listFiles = @$ref2List;
	my $allFiles = join(" ",@listFiles);
	
	`cat $allFiles > $outFile`;
}


sub fileToBDB {
	my($inFile,$bdbFile) = @_;
	
	my %hash = ();
	tie %hash, "BerkeleyDB::Btree", -Filename => "$bdbFile", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";

	open(FI,$inFile) || die "Error opening inFile '$inFile' in fileToBDB file\n";
	while (my $line = <FI>) {
		$line =~s/\s+$//;
		my($gene,$acc,$chr,$exonCDS,$species,$length) = (split /\t/,$line)[0,1,2,3,4,5];
		my $key = "$gene#$acc#$chr#$exonCDS#$species";
		$hash{$key} = $length;
	}
	close FI;
	untie %hash;
}
