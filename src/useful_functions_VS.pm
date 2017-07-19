#!/sw/bin/perl

#### Virag Sharma, July 2013.
#### This package contains certain trivial but nevertheless very useful functions.

package useful_functions_VS;
use strict;
use warnings;
use Scalar::Util::Numeric qw(isint);

use Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(seq_multiple_3 find_stop_codon find_stop_codon_all find_particular_codon change_cds split_codons seq2codon 
find_pattern fasta_hash randomizeFile RemoveHeader nol delete_if_empty SortAndMerge GetSpeciesList CreateBedFile createNRList 
ParseParameters getRandomInsert createOutDir getHashMax getHashMin create_temp revComp getSLineData getELineData getILineData readBDB secureWriteBDB);

our $genomePath = $ENV{'genomePath'};

####################################################################################################################
####      Add N's at the end of the sequence so that its length becomes a multiple of 3. This is useful 	   #####
####################################################################################################################

sub seq_multiple_3 {	
	my $sequence = shift;
	my $rem = length($sequence)%3;
	
	if ($rem == 1) { $sequence=$sequence."NN"; }
	elsif ($rem == 2) { $sequence=$sequence."N";	}
	
	return $sequence;
}

########################################################################
#### 	Find the first/last/all stop codon(s) in an input sequence	####
########################################################################

sub find_stop_codon {
	my($seq,$argument) = @_; 	### Takes two arguments: the sequence and the label for the stop codon - first, last or all.
	
	$seq = seq_multiple_3($seq);
	
	my $i = my $s = 0;
	my $stop = "";
	my @stop_list = ();
	
	while ($i < length($seq)-2)
	{
		my $codon = substr($seq,$i,1).substr($seq,$i+1,1).substr($seq,$i+2,1);
		
		if ( ($codon eq "TAA") or ($codon eq "TAG") or ($codon eq "TGA") )
		{
			$stop = $i;
			
			if ($argument eq "first") { last; }
			elsif ($argument eq "all") {$stop_list[$s] = $i; $s++}
		}
		
		$i = $i+3;
	}
	
	return $stop if ( ($argument eq "first") or ($argument eq "last") );
	return (\@stop_list);
}

#########################################################################
####   Find all stop codons in an input sequence (in all frames) 	 ####
#########################################################################

sub find_stop_codon_all {
	my $seq = shift;			
	
	my @stop_list = ();
	my $s = my $j = 0;

	while ($j <=2) {
		my $sequence_frame = substr($seq,$j);
		$sequence_frame = seq_multiple_3($sequence_frame);
		
		my $i = 0;
		while ($i < length($sequence_frame)-2) {
			my $codon=substr($sequence_frame,$i,1).substr($sequence_frame,$i+1,1).substr($sequence_frame,$i+2,1);
			
			if ( ($codon eq "TAA") or ($codon eq "TAG") or ($codon eq "TGA") ) {
				if ($j == 0)
					{ $stop_list[$s] = $i;
				} elsif ($j == 1)
					{ $stop_list[$s] = $i+1;
				} else  {
					$stop_list[$s] = $i+2;
				}
				
				$s++
			}
			$i = $i+3;
		}
		$j++;
	}
	
	@stop_list = sort {$a <=> $b} @stop_list;				## Just to ensure that the stop codons are listed in ascending order/sorted, it helps.
	return (\@stop_list);
}

########################################################################################
#### 	Find the first instance of a user specified codon in an input sequence 		####
########################################################################################

sub find_particular_codon {
	my ($seq,$check_codon) = @_;	
	$seq = seq_multiple_3($seq);
	
	my $i = 0;
	my $stop = "";
	
	while ($i < length($seq)-2) {
		my $codon=substr($seq,$i,1).substr($seq,$i+1,1).substr($seq,$i+2,1);
		if  ($codon eq $check_codon) {
			$stop = $i;
			last;
		}
		
		$i = $i+3;
	}
	
	return $stop;
}

###########################################################################################################
####   	Change the coordinates given a starting position, the genomic coordinates and the strand       ####
###########################################################################################################

sub change_cds {
	my ($stop_pos,$genomic_cds,$strand) = @_;
	my $report_cds = "";
	
	die "Strand not defined for $strand\n" if ( ($strand ne "+") and ($strand ne "-") );
	die "Genomic cds $genomic_cds is less than the specified position $stop_pos\n" if ( ($genomic_cds < $stop_pos) and ($strand eq "-") );
	
	if ($strand eq "+") {
		$report_cds = $genomic_cds + $stop_pos;
	} else {
		$report_cds = $genomic_cds - $stop_pos;
	}
	
	return $report_cds;
}

############################################################################################################
####   Given a nucleotide sequence, split it into codons and separate the codons with a white space 	####
############################################################################################################

sub split_codons {
	my $seq = shift;
	my $seq_new = seq_multiple_3($seq);
	
	my $add_nt = length($seq_new)-length($seq);
	
	my $i = my $c = 0;
	my @codon_list = ();
	
	while ($i < length($seq_new)) {
		my $codon = substr($seq_new,$i,1).substr($seq_new,$i+1,1).substr($seq_new,$i+2,1);
		$codon_list[$c] = $codon;
		
		$c++;
		$i = $i+3;
	}

	my $seq_codons = join(" ",@codon_list);
	my $seq_codons_trim = substr($seq_codons,0,(length($seq_new)-$add_nt));
	
	return $seq_codons_trim;
}

####################################################################################################################
####   Very similar to split_codons. However, it returns a reference to the array that contains all the codons 	####
####################################################################################################################

sub seq2codon {
	my $seq = shift;
	my $seq_new = seq_multiple_3($seq);
	
	my $add_nt = length($seq_new) - length($seq);
	
	my $i = my $c = 0;
	my @codon_list = ();
	
	while ($i < length($seq_new)) {
		my $codon = substr($seq_new,$i,1).substr($seq_new,$i+1,1).substr($seq_new,$i+2,1);
		$codon_list[$c] = $codon;
		
		$c++;
		$i = $i+3;
	}

	return \@codon_list;
}


########################################################################################################################################################################################
####		Given a pattern and a sequence, the function returns a reference to an array than contains a list of all the locations of the pattern in the input sequence				####
########################################################################################################################################################################################


sub find_pattern {				### Very useful function. Input arguments -> a sequence and a pattern/k-mer that you want to find in the input sequence.								### Output -> a reference to an array that contains the locations of the pattern in the sequence.
								### Very nostalgic
	my ($input_seq,$pattern) = @_;
	my $result = index($input_seq,$pattern);

	my @return_list = ();
	my $rl = 0;

	if ($result >-1) {
		my $r_chop = 0;
		my $r      = $result;
		
		do {
			$return_list[$rl] = $r;
			$rl++;
			
			$r_chop     = $r_chop + $result + length($pattern);                    
			my $input_seq_tr = substr($input_seq,($result + length($pattern)) );

			$result	= index($input_seq_tr,$pattern);
			$r	= $result+$r_chop;
			
			$input_seq	= $input_seq_tr;
		}
		while ($result !=-1);
	}
	
	return(\@return_list);
}

####################################################################################################################################
####	Given a fasta file, the function returns a reference to the hash which contains the key-value pairs for the fasta file	####
####################################################################################################################################

sub fasta_hash {
    my $file_input = shift;
    
    my $id = my $seq = "";
    my %seq_hash = ();

    open(FIS,"$file_input") || die "Error opening input file '$file_input' in fasta_hash function\n";
    while (my $line=<FIS>) {
		chomp $line;

		if ($line =~/>/) {
			if ($seq ne "") {
				$seq_hash{$id} = $seq;
				$id = $seq = "";
			}

			$id =$line;
			$id =~s/>//;
		} else {
			$seq = $seq.$line;
		}
    }
    close FIS;

    $seq_hash{$id} = $seq;
    return(\%seq_hash);
}

#######################################################################
##### This function finds the number of lines in the input file   #####
#######################################################################

sub nol
{
	my $file = shift;
	my $nol = `wc -l < $file`; chomp $nol;
	return $nol;
}

###########################################################################################
#####    Given a file/directory, this function checks if the input is empty i.e. the  #####
#####    file has no lines or the directory has no files. If yes, then delete it.     #####
###########################################################################################

sub delete_if_empty {
        my ($input,$type) = @_;
        die "Input type $input not defined for delete_if_empty function\n" if ( ($type ne "f") && ($type ne "d") );
        
        if ($type eq "f") {
                my $nol = nol($input);                                                             ## Get the number of lines for this file
		`rm $input` if ($nol == 0);                               ## Delete this file if there is actually nothing in this file
        } elsif ($type eq "d") {
                my $no_files = `ls $input|wc -l|tr -d "\n"`;                      ## Get the number of files present in this directory
                `rm -r $input` if ($no_files == 0);                               ## Delete if the number of files == 0
        }
}

##################################################################################
##### This function takes a bed file as an input. Then it sorts the BedFile, #####
##### then merges the continuous elements in the sorted file. Finally, it    #####
##### renames the new sorted plus merged file to the input file.             #####
##################################################################################

sub SortAndMerge         {                                       
	my $file = shift;
	die "The input file $file does not exist\n" if ($file eq "" || ! -e $file);

	my $sortedBed  = create_temp("f");
	my $mergedBed  = create_temp("f");
	my $call = "bedSort $file $sortedBed";
	system($call) == 0 || die "Error running '$call' in SortAndMerge function\n";
	
	my $callMerge = "mergeBed -i $sortedBed > $mergedBed";					### merge the continuous elements in the bed file e.g. 100-200 and 200-210 becomes 100-210. 
	system($callMerge) == 0 || die "Error running '$callMerge' in SortAndMerge function\n";
																			### This is important because if there is an insertion in between 200 and 201, it will be 
																			### captured only if the element is presented as 100-210 and not as two separate elements.
	`mv $mergedBed $file`;
	`rm $sortedBed`;		## Delete this temp file.
}

sub GetSpeciesList {   ### This function returns a reference to an array which contains all species present in the maf.
        my $MafIn = shift;              ### The only argument is the MAF block.
        my @species_array = ();
        @species_array = `mafSpeciesList $MafIn stdout`;
        chomp(@species_array);

        die "ERROR: mafSpeciesList must have failed because the list of species is empty\n" if (scalar @species_array < 1);
        $" = " ";       ### separate array elements by a space in the output

        return (\@species_array);
}

sub CreateBedFile {  ## Given a chromosome, a start coordinate and a stop coordinate, this function writes these fiels to a bedFile and returns it.

	my($chr,$start,$stop) = @_;

	my $fileTmp = create_temp("f");
	open(FO,">$fileTmp") || die "Error writing to tempFile in CreateBedFile function";
	print FO "$chr\t$start\t$stop\n";
	close FO;

	return $fileTmp;
}

sub createNRList {	## Given an array that has some values that could be present multiple times, the function returns an array where each value is present only once
					## set intFlag to "T" if the input list has integers/numbers and then the function returns a sorted list.
	
	my ($ref2List,$intFlag) = @_;
	my @list = @$ref2List;
	
	my %tmpHash = map{$_ => 1}@list;
	my @returnList = keys(%tmpHash);
	
	@returnList = sort {$a <=> $b} @returnList if ($intFlag eq "T");
	return \@returnList;
}

sub ParseParameters {		## Parse parameters for the different modules of the GeneLoss pipeline
	my($fmt,$input,$pMammals,$speciesList,$reference,$alignment) = @_;

	my $GeneLossConfig = "$ENV{'GeneLossPipeCode'}/GeneLoss.config";
	
	## Get a list of accessions
	my @accList = ();
	if ($fmt eq "lst") {
		open(FIL,"$input");
		@accList = <FIL>;
		close FIL;
	} else {	
		push(@accList,$input);
	}

	## Now get a list of species

	my @placental_mammals = ();
	my $p = 0;
	
	if ($speciesList ne ""){
		@placental_mammals = split(",",$speciesList);
	}else{
		my $searchTerm = "$reference\_$alignment";
		$searchTerm = "$reference\_pMammals\_$alignment" if ($pMammals eq "T");
		
		open(FI,"$GeneLossConfig") || die "Error opening geneLossConfig '$GeneLossConfig' file\n";
		while (my $line = <FI>) {
			$line =~s/\s+$//;
			my ($species,$speciesTerm) = (split /\t/,$line)[0,1];
			push(@placental_mammals,$species) if ($speciesTerm eq $searchTerm);
		}
		close FI;
	}	
	
	@placental_mammals = grep {$_ ne $reference} @placental_mammals; #### Remove the reference species from the list of species you want to look into.
	
	return(\@accList,\@placental_mammals);
}

sub getRandomInsert {	## This function generates an Insert of a given length (its length is a multiple of 3 and it does not contain any in-frame stop codon)

	my ($length,$type) = @_;
	
	die "Argument '$type' not recognized by getRandomInsert function\n" if ( ($type eq "CDS") && ($type ne "non-CDS") );
	
	### Make sure that the length is a multiple of 3
	my $rem = $length%3;
	$length = $length - $rem;
	
	my @bases = qw(A T C G);
	my %stopCodon = ("TAA" => "T", "TAG" => "T", "TGA" => "T");
	my $string = "";
	
	my $i = 0;
	while ($i < $length) {
		my $codon = $bases[int(rand(4)-1)].$bases[int(rand(4)-1)].$bases[int(rand(4)-1)];
	
		if ($type eq "CDS") {		## Only worry about this if you want a random CDS insert where stop codons are disallowed.
			if (! exists $stopCodon{$codon}) {
				$string = $string.$codon;
			} else {
				$i = $i - 3;
			}
		}
	
		$i = $i+3;
	}
	
	return $string;
}

sub createOutDir {  ### Suppose you want to write to a file /scratch/users/vsharma/test/testOutput. This function will create the /scratch/users/vsharma/test/ directory
	my $file = shift;
	my $dirName = `dirname $file`;
	chomp $dirName;	
	`mkdir -p $dirName` if ($dirName ne ".");
}


sub getHashMax {  ## given a hash, this function returns the key and the associated value that is maximum among all the values in the hash
	my $ref2Hash = shift;
	my %hash = %$ref2Hash;
	
	my @keysList = keys(%hash);
	my $keyMax = $keysList[0];		## Assign this to the first value in the hash.
	
	my $valueMax = $hash{$keyMax};
	foreach my $keys(@keysList) {		## Now check which key has the maximum value
		my $value = $hash{$keys};
		if ($value > $valueMax) {
			$keyMax = $keys;
			$valueMax = $value;	
		}
	}
	return ($keyMax,$valueMax);
}

sub getHashMin {  ## given a hash, this function returns the key and the associated value that is minimum among all the values in the hash
	my $ref2Hash = shift;
	my %hash = %$ref2Hash;
	
	my @keysList = keys(%hash);
	my $keyMin = $keysList[0];		## Assign this to the first value in the hash.
	
	my $valueMin = $hash{$keyMin};
	foreach my $keys(@keysList) {		## Now check which key has the minimum value
	
		my $value = $hash{$keys};
		if ($value < $valueMin) {
			$keyMin = $keys;
			$valueMin = $value;	
		}
	}
	return ($keyMin,$valueMin);
}

#############################################################################################################################################################
# creates a temporary file or directory, if passed with "f" argument -> create a temp file, if passed with "d" argument -> create a temp directory.
#############################################################################################################################################################
sub create_temp {
	my $arg_1=my $temp_dir="";
	
	$arg_1=shift;
	
	if ($arg_1 eq "d") {
		$temp_dir=`mktemp -d`;
	} elsif ($arg_1 eq "f") {
		$temp_dir=`mktemp`;
	} else {
		die "Argument not defined for create_temp function\n";
	}
	
	$temp_dir=~s/\s+$//;
	return ($temp_dir);
}

#########################################################################################################
# return reverse complement seq
#########################################################################################################
sub revComp {
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}


#############################################################
# split an e line like
# 		e anoCar1.scaffold_563               657656 18140 +    818251 I	
# into the fields
#############################################################
sub getELineData {
	my $line = shift; 
	
	if ($line =~ /^e (.*)/) {		# empty part
		# anoCar1.scaffold_563               657656 18140 +    818251 I
		my @f = split(/[ ]+/, $1); 
		if ($#f != 5) {
			die "ERROR in getELineData: cannot parse 6 elements from $line\n";
		} 
		my ($src, $start, $size, $strand, $srcSize, $status) = (@f)[0,1,2,3,4,5];
		my ($species) = (split(/[\.]/, $src))[0];
		my $chr = substr($src, length($species)+1, length($src));	# given "tupBel1.scaffold_142368.1-170510" $chr is "scaffold_142368.1-170510"

#		print "e line: ($blockNo, $species, $chr, $start, $size, $strand, $srcSize, $status)   [[$line]]\n" if ($verbose);

		return ($species, $chr, $start, $size, $strand, $srcSize, $status);
	}else{
		die "ERROR: call getELineData with no e line: $line\n";
	}
}

#############################################################
# split an i line like
# 		i ornAna1.Contig17774             C 0 I 49
# into the fields
#############################################################
sub getILineData {
	my $line = shift; 
	
	if ($line =~ /^i (.*)/) {		# information part
		# ornAna1.Contig17774             C 0 I 49
		my @f = split(/[ ]+/, $1); 
		if ($#f != 4) {
			die "ERROR in getILineData: cannot parse 5 elements from $line\n";
		} 
		my ($src, $statusUp, $countUp, $statusDown, $countDown) = (@f)[0,1,2,3,4];
		my ($species) = (split(/[\.]/, $src))[0];
		my $chr = substr($src, length($species)+1, length($src)); 	# given "tupBel1.scaffold_142368.1-170510" $chr is "scaffold_142368.1-170510"

#		print "i line: ($blockNo, $species, $chr, $statusUp, $countUp, $statusDown, $countDown)   [[$line]]\n" if ($verbose);

		return ($species, $chr, $statusUp, $countUp, $statusDown, $countDown);
	}else{
		die "ERROR: call getILineData with no i line: $line\n";
	}
}

#############################################################
# split an s line like
# 		s ornAna1.Contig17774                  7736 10 +     18306 CTGGG----GCTGT
# into the fields
#############################################################
sub getSLineData {
	my $line = shift; 
	
	if ($line =~ /^s (.*)/) {		# sequence part
		# ornAna1.Contig17774                  7736 10 +     18306 CTGGG----GCTGT
		my @f = split(/[ ]+/, $1); 
		if ($#f != 5) {
			die "ERROR in getSLineData: cannot parse 6 elements from $line\n";
		} 
		my ($src, $start, $size, $strand, $srcSize, $seq) = (@f)[0,1,2,3,4,5];
		# few sanity checks
		if (! isint($start) || ! isint($size) || ! isint($srcSize)) {
			die "ERROR in getSLineData: start/size/srcSize are not integers in $line\n";
		}
		if ($strand ne "+" && $strand ne "-") {
			die "ERROR in getSLineData: strand is neither + nor - in $line\n";
		}
		my ($species) = (split(/[\.]/, $src))[0];
		my $chr = substr($src, length($species)+1, length($src));	# given "tupBel1.scaffold_142368.1-170510" $chr is "scaffold_142368.1-170510"

#		print "s line: ($blockNo, $species, $chr, $start, $size, $strand, $srcSize, $seq)  [[$line]]\n" if ($verbose);

		return ($species, $chr, $start, $size, $strand, $srcSize, $seq);
	}else{
		die "ERROR: call getSLineData with no s line: $line\n";
	}
}


##################################################################
# get value from a BDB file for a given key
##################################################################
sub readBDB {
	my ($file, $key) = @_;

	my $Path = `pwd`;  chomp($Path);
	$key =~ s/$Path//g;		# get rid of the absolute path until chr*/CNE*/
	$key =~ s/\/\//\//g;		# get rid of // in e.g. chr2/CNE.2//CNE.2  
	$key = substr($key, 1, length($key)) if (substr($key, 0,1) eq "/");			# get rid of leading /

	my %h;
	tie %h, "BerkeleyDB::Btree", -Filename => "$file", -Flags => DB_RDONLY or die "Cannot open BerkeleyDB: $file   $! $BerkeleyDB::Error\n";
	my $value = "";
	if (exists $h{$key}) {
		$value = $h{$key};
	}
	untie %h;
	return $value;
}


##################################################################
# lock the BDB file, then write and unlock
##################################################################
sub secureWriteBDB {
	my ($file, $key, $value, $useLock, $dir) = @_;

	my $Path = `pwd`;  chomp($Path);
	$key =~ s/$Path//g;		# get rid of the absolute path until chr*/CNE*/
	$key =~ s/\/\//\//g;		# get rid of // in e.g. chr2/CNE.2//CNE.2  
	$key = substr($key, 1, length($key)) if (substr($key, 0,1) eq "/");			# get rid of leading /
	
	my %hash;
	tie %hash, "BerkeleyDB::Btree",	-Filename => "$file", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";
	$hash{$key} = $value;
	untie %hash;

	# now read the same key from the BDB and see if we get exactly the same. Otherwise die. 
	my %h;
	tie %h, "BerkeleyDB::Btree", -Filename => "$file", -Flags => DB_RDONLY or die "Cannot open BerkeleyDB: $file   $! $BerkeleyDB::Error\n";
	my $valueRead = "";
	if (exists $h{$key}) {
		$valueRead = $h{$key};
	}
	untie %h;
	die "ERROR in secureWriteBDB: reading $key from $file does not give the value that was stored\nINPUT $value\nOUTPUT $valueRead\n" if ($value ne $valueRead);
}
