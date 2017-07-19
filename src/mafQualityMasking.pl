#!/sw/bin/perl

### Script to cleanse maf blocks for each chromosome. 
### By cleanisng, I mean replace each base of low quality with an "N", so that subsequent results regarding these bases which are converted to "N"s are treated with caution.
### Virag Sharma, March 2013.

use strict;
use warnings;
use FileHandle;
use BerkeleyDB;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;
use MyKentFunctions;

my $genomePath = $ENV{'genomePath'};
die "ERROR: directory specified in environment variable genomePath (set to $genomePath) does not exist\n" if (! -d $genomePath);
my $path_assembly = "$genomePath/gbdb-HL"; ### Directory which contains the 2 bit files (for genome)

$| = 1;		# == fflush(stdout)

# input options/parameters
my $verbose = 0;
my $quiet = 0;
my $mafIn =""; 	## the input MAF file.
my $mafOut ="";   	## the output file.
my $listFiles = "";     ## the path to the file that contains the path to the user specified masked files. 

sub usage 
{
	die "usage: $0 -input Input_maf_file -output Output_file -listFiles (optional, only if you want a non-default quality masking. Specify a tab separated file which has 2 fields --> the assembly and the path to the User Specified Masked File) -quiet (does not print anything to the screen) \n";
}

GetOptions ("v|verbose"  => \$verbose, "q|quiet" => \$quiet, "input=s" => \$mafIn, "output=s" => \$mafOut, "listFiles:s" => \$listFiles) || usage();

usage() if ($mafIn eq "");
### Now check if all the input parameters have been correctly specified.
die "Input maf file i.e. $mafIn does not exist\n" if (! -e $mafIn);
die "Output maf file is not specified\n" if ($mafOut eq "");

## Everything seems ok
## Check if the file has nothing but the header lines and no maf blocks (this is the case with a few chromosomes like: chrUn_gl000226.maf ->hg19_46way/100way), for this get the number of lines in the maf
my $getNumberOfLines = "wc -l $mafIn|awk '{ print \$1}'|tr -d \"\n\"";
my $nol = `$getNumberOfLines`;
die "ERROR: $getNumberOfLines call failed, cannot get the number of lines in maf\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);

if ($nol == 1){		## If yes, then copy this file to the output directory and do nothing
	`cp $mafIn $mafOut`;
	exit; ## And just exit.
}

### Run mafSpeciesList to get a list of species present in the input maf block.	
my @species_array=();			### create a species_array which contains the list of all the species present in the maf block.
@species_array = `mafSpeciesList $mafIn stdout`; 
chomp(@species_array); 
die "ERROR: mafSpeciesList must have failed because the list of species is empty\n" if (scalar @species_array < 1); 
$" = " ";	# separate array elements by a space in the output
print "List of species in the maf: @species_array\n" if ($verbose);

my %hash_2bit_file = ();
my %assembly_quality_ua = ();    ## %assembly_quality_ua [Assembly_quality_Unavailable] is a hash where the key is the species and the value is either "T" 
                                 ## (for species with real quality.2bit files) or "F" (for species with fake quality.2bit files).
                                                        
## For the moment, let us believe that no species has a quality.2bit file. In that case, all the quality.2bit files will be faked.
## Each key/species in the hash %assembly_quality_ua will have the value "T".

foreach my $species(@species_array)
{
	my $sp_2bit_file = "$path_assembly/$species/$species.2bit";
	$hash_2bit_file{$species} = $sp_2bit_file;
	$assembly_quality_ua{$species} = "T";
}

### Now if you want to specify user specified 2 bit files (produced by a different quality threshold), then do the following:
if ($listFiles ne "")
{
	open(FIL,$listFiles) || die "Error !! Cannot open the listFiles '$listFiles' specified\n";
	while (my $line = <FIL>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
		
		my $species = $tmp[0];
		my $twoBitFile_userSpecified = $tmp[1];
		
		if (! -e $twoBitFile_userSpecified)
		{
			die "The user specified twoBit file --> $twoBitFile_userSpecified for $species does not exist\n";
		}
		else
		{
			$hash_2bit_file{$species} = $twoBitFile_userSpecified;
			$assembly_quality_ua{$species} = "F";
		}
	}
	close FIL;
}

else	## check if the simpler quality.2bit files exist (default where every base that is below the quality threshold is masked, nothing else)
{
	foreach my $species(@species_array)
	{
		my $sp_2bit_file = "$path_assembly/$species/$species.quality.2bit";
		
		$hash_2bit_file{$species} = $sp_2bit_file if (-e $sp_2bit_file);
		$assembly_quality_ua{$species} = "F";
	}
}

## Check that an exisiting 2bit file has been assigned to each assembly. If it is not a quality.2bit file, then it should be the "fake" genome.2bit file
my @filesUnavailable = ();
foreach my $species(@species_array)
{
	my $twoBitFile = $hash_2bit_file{$species};
	push(@filesUnavailable,$twoBitFile) if (! -e $twoBitFile);
	
}

if (scalar(@filesUnavailable) != 0){
	my $string = join(",",@filesUnavailable);
	die "No 2 bit files (neither whole genome, nor any quality masked) are found for these assemblies --> $string\n";
}

### Now figure out which files are RealQualityMasked files and which files are links, like hg19.quality.2bit -> hg19.NoQualityScore.2bit
my %hashFakeFile = ();
foreach my $assembly(keys(%hash_2bit_file))
{
	my $twoBitFileQM = $hash_2bit_file{$assembly};
	$hashFakeFile{$assembly} = "T" if (-l $twoBitFileQM); ## Dump those assemblies in hashFakeFile where the 2bit file is a link.
}

### If all the 2bit files are found, proceed further.
#### Now the main body of code which actually cleanses the maf blocks
print "All 2 bit files are found. Now ready to cleanse maf\n" if (! $quiet);

my $temp_dir = create_temp("d");			### First of all, create a temporary directory, create two sub-directories within this temp diretory as well.
`mkdir $temp_dir/bed_files`;
`mkdir $temp_dir/fasta_files`;

my %sp2File = my %sp_bed_file = my %sp_fasta_file = ();

print "open a tmp bed file $temp_dir/bed_files/\$species.bed for every species ..." if ($verbose);
foreach my $species (@species_array)
{
	if (! exists $hashFakeFile{$species}){
		my $fh = FileHandle->new();        # need to create a new filehandle.
		$sp2File{$species}=$fh;
	
		my $bed_file = "$temp_dir/bed_files/$species.bed"; 
		$sp_bed_file{$species} = $bed_file;
		$fh->open(">$bed_file") || die "ERROR: cannot write to $bed_file\n";
	}else{
		my $fh = FileHandle->new();        # need to create a new filehandle.
		$sp2File{$species}=$fh;
	
		my $fasta_file = "$temp_dir/fasta_files/$species.fa"; 
		$sp_fasta_file{$species} = $fasta_file;
		$fh->open(">$fasta_file") || die "ERROR: cannot write to $fasta_file\n";
	}
}
print "DONE\n" if ($verbose);

my $deletedSeqFasta = create_temp("f");
open(FOD,">$deletedSeqFasta");

#### Now create temporary bed files per species which would contain the coordinates for each block from the input maf. 1 bed file per species/assembly.
#### We read the maf file for the first time here. 
print "Read the $mafIn and get the coordinates for every aligning seq for all species...\n" if ($verbose);
open(FI,"$mafIn") || die "Error opening input file which is $mafIn\n";
while (my $line=<FI>)
{
	$line=~s/\s+$//;
	
	if ($line =~/^s/){
		my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line); 
		my $name = "$species\_$chr\_$start\_$size\_$strand";		#### The name is the "unique key" here

		my $fh = $sp2File{$species};
		
		if ($size != 0){
			if (! exists $hashFakeFile{$species}){
				if ($strand eq "-"){
					my $stop = $srcSize-$start;
					$start = $stop-$size;
			
					print $fh "$chr\t$start\t$stop\t$name\n";
					print "\t$species --> write $chr\t$start\t$stop\t$name  for - strand aligning s line: $line\n" if ($verbose);
				}else{
					my $stop = $start+$size;
					print $fh "$chr\t$start\t$stop\t$name\n";
					print "\t$species --> write $chr\t$start\t$stop\t$name  for + strand aligning s line: $line\n" if ($verbose);
				}
			}else{
				$seq = uc($seq); ## Capitalize the sequence, Note that if you had run 2BitToFa on a fake quality 2bit file, then it would be all Caps. Now since you do notrun 2BitToFa, you need to do allCaps operation here.
				print $fh ">$name\n$seq\n";
			}
		}else{
			print FOD ">$name\n$seq\n";
		}
	}
}
close FI;
close FOD; 
print "writing the coordinates (for assemblies with a Quality.2bit file) and sequences (for assemblies with a fake Quality.2bit file) - DONE\n" if ($verbose); 

print "Now close the file handles and get the seqs with twoBitToFa ..." if ($verbose);
foreach my $species (keys %sp2File) 			### close all file handles, also run twoBitToFa at the same time to extract the sequence for the coordinates present in the bed file.
{
	my $fh = FileHandle->new();
	$fh = $sp2File{$species};
	$fh->close;

	if (! exists $hashFakeFile{$species}){
		my $bed_file 			= "$temp_dir/bed_files/$species.bed"; 
		my $fasta_file 			= "$temp_dir/fasta_files/$species.fa";	
		my $quality_masked_file = $hash_2bit_file{$species};
	
		my $status = system("bedSort $bed_file stdout| twoBitToFa $quality_masked_file -bed=stdin $fasta_file"); ## First run bedSort on the bed file, subsequently run twoBitToFa on the sorted bed file.
																											## a sorted bed file is processed faster by twoBitToFa than an unsorted one.		
		die "Error running twoBitToFa on bed file for $species" if ($status != 0);
		print "twoBitToFa run successful for $species \n" if ($verbose);
	}
}
print "DONE\n" if ($verbose);

print "Now concatenate all species.fasta files and create a BDB hash ...\n" if ($verbose);
my $all_fasta_combined = create_temp("f");

my $status = system("cat $temp_dir/fasta_files/*.fa $deletedSeqFasta > $all_fasta_combined");	### Create a BIG fasta file from all fasta files.
die "Error : cannot combine all fasta files present in $temp_dir/fasta_files/ " if ($status != 0);
print "All fasta file $all_fasta_combined created\n" if ($verbose);

###### Read the concatenated fasta file and convert it to a hash which is stored in a BerkeleyDB file	######	
my $bdb_file = create_temp("f");

my %hash = ();
tie %hash, "BerkeleyDB::Btree",    -Filename => "$bdb_file", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";

open(FIS,"$all_fasta_combined") || die "Error opening concatenated $all_fasta_combined \n";
my $seq = my $id = "";

while (my $line=<FIS>)
{ 
	$line =~s/\s+$//;

	if ($line =~/>/){
		if ($seq ne ""){
			$hash{$id} = $seq;	# store the previous seq and ID in the hash
			print "\tstore in $bdb_file: $id -> $seq\n" if ($verbose);
			$id = $seq = "";
		}

		$id = $line;			# ID of current sequence
		$id =~s/>//;
	}else{
		$seq = $seq.$line;
	}
}
close FIS;

$hash{$id} = $seq;	# store the last sequence
print "\tstore in $bdb_file: $id -> $seq\n" if ($verbose);
untie %hash;
print "Writing $bdb_file DONE\n" if ($verbose);	
## Fasta file read and BDB file created.

`rm -r $temp_dir`; ### Get rid of the temporary directory, not needed anymore.
`rm $all_fasta_combined $deletedSeqFasta`; ### Also delete this combined fasta file, the id-sequence pairs have already been stored in a hash.
		
## Now read the BDB file created in the previous step into a hash called %hash_read.
		
my %hash_read = ();
tie %hash_read, "BerkeleyDB::Btree", -Filename => "$bdb_file", -Flags => DB_RDONLY or die "Cannot open BerkeleyDB file   $! $BerkeleyDB::Error\n";
		
print "All processing done. Now writing to the output file $mafOut\n" if (! $quiet);

## At this stage, read the entire file again ($mafIn) and replace the old sequence with the new/quality masked sequence.
print "Now reading $mafIn again and masking every sequence ...\n" if ($verbose);
open(FIM,$mafIn) || die "Error opening chromosome maf file\n";
open(FO,">$mafOut") || die "Error opening output file which is $mafOut \n";

while (my $seq_line=<FIM>)
{
	
	if ($seq_line =~/^s/) {
		$seq_line =~s/\s+$//;
		my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($seq_line); 
		
		my $id_old = "$species\_$chr\_$start\_$size\_$strand";		## this is the key that is used in the hash created earlier.

		die "ERROR: this key $id_old does not exist in the $bdb_file hash\n" if (! exists $hash_read{$id_old});
		my $seq_new = $hash_read{$id_old};	     			## and here we have the "new sequence".
		my $seq_old = $seq;
		
		if ($strand eq "-"){
			$seq_new = revComp($seq_new);
		}
		
		my $seq_new_gaps = "";
		
		##### Finally, the entire thing comes down to this point
		
		if ( ($assembly_quality_ua{$species} eq "T") or ($seq_new !~/[acgtn]/) ){		## Situation 1: if the species has no quality scores available (or fake quality.2bit files) or if the sequence has no lower case characters, just use the old sequence but capitalize everything
			$seq_new_gaps = uc($seq_old);
		}else{ 											## Situation 2: the sequence has some lower case characters, substitute such characters with "N"
			$seq_new 	  =~tr/atcgn/NNNNN/;
			$seq_new_gaps = sub_introduce_gaps($seq_new,$seq);
			print "\t$species: replace this seq by this masked seq: \n\t\t$seq\n\t\t$seq_new_gaps\n" if ($verbose);
		}
		
		my $pos = rindex($seq_line, $seq);			    ## start searching at the end and get the position from where the sequence begins in the $seq_line. Remember "rindex" returns the last occurrence of the pattern in the string while "index" returns the first.
		my $allButSeq = substr($seq_line, 0, $pos);    ## get the contents of the line upstream of the position from where the sequence begins. Store it in variable $allButSeq
		my $print_line = $allButSeq.$seq_new_gaps;		## concatenate $allButSeq and $seq_new_gaps.
		
		print FO "$print_line\n";
		print "\t$species: replace this entire line by this: \n\t\t$seq_line\n\t\t$print_line\n" if ($verbose);
	}else{
		print FO "$seq_line";
	}
}
close FIM;
close FO;

`rm $bdb_file`;  ## Do not see any need to DIE here since this is the last step.
print "$mafIn successfully masked and written to $mafOut\n" if (! $quiet);
	
## End of the main body of the code. 
## Sub-routines follow.

sub sub_introduce_gaps			## Subroutine to introduce gaps in the new sequence (quality masked) based on the position of gaps in the old (but aligned with gaps) sequence.
{
	my ($seq_masked,$seq_aligned) = @_;		## the new sequence and the old aligned sequence

	my @seq_masked_array  = split(//,$seq_masked);
	my @seq_aligned_array = split(//,$seq_aligned);
	
	my $masked_aligned_seq = "";	
	my $ug_index=0;  # $ug_index refers to the index in the ungapped sequence.
	
	for(my $j=0; $j<@seq_aligned_array; $j++)
	{
		if ($seq_aligned_array[$j] eq "-"){
			$masked_aligned_seq .="-";
		}else{
			$masked_aligned_seq .=$seq_masked_array[$ug_index];
			$ug_index++;
		}
	}
	return $masked_aligned_seq;
}
