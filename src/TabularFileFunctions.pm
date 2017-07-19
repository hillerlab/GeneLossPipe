#!/sw/bin/perl

package TabularFileFunctions;
use strict;
use warnings;
use Exporter;

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;
use MyKentFunctions;

use lib "$ENV{'GeneLossPipeCode'}";
use useful_functions_VS;
use FixInsertsFromILines;

our @ISA = ('Exporter');
our @EXPORT = qw(GetHashFromMafAndPrintIt GetMafLength GetBlockLengthCum GetBlockBoundaries PlaceInsertions FindDeletionsAndPlaceThem FixMissingData FixMissingDataPositions 
ReverseFile InvertStrand IncrementCds); 

sub GetHashFromMafAndPrintIt	### This is the most important function, does practically everything.
{	
	my ($mafFile,$reference,$outDir,$strand_input,$Ref2CdsAnnotationHash,$Ref2CodonIndexHash) = @_;		## Arguments --> 1) bedFile 2) Reference species 3) Output directory 4) StrandReference

	my %genome_cds_nature_hash = %$Ref2CdsAnnotationHash;
	my %codon_index_hash       = %$Ref2CodonIndexHash;
	
	##### Now get a list of species present in the maf.
	my $Ref2SpeciesList = GetSpeciesList($mafFile);
	my @species_list = @$Ref2SpeciesList;	
	
	my %species_hash = ();
	foreach my $species (@species_list)
	{
		$species_hash{$species} = "T;"
	}
	
	### Look at the insertions that come from the "i" lines. Insertion_iLine Business
	
	my $Ref2HashInsertions = FixInsertions($mafFile,$strand_input,$reference);
	my %tmpHash = %$Ref2HashInsertions;
	
	my %SpeciesWithInserts = ();
	
	foreach my $keys(keys %tmpHash)
	{
		$SpeciesWithInserts{$keys} = "T";
	}

	### Now process the  i lines where you see "M" as the LeftStatusCharacter.
	
	my $ref2MissingData = FixMissingData($mafFile,$reference);
	my %HashMissingData = %$ref2MissingData;

	## get a list of species which have this anomaly
	my %MissingDataSpecies = ();
	
	foreach my $species(keys %HashMissingData)
	{
		$MissingDataSpecies{$species} = "T";
	}
	####
	
	my $length_maf_file = GetMafLength($reference,$mafFile);
	
	my %hash_genomic_cds = my %hash_base = my %hash_chr = my %hash_strand = (); #### Initialize 4 2-dimensional hashes: one for genomic coordinates, one for base/nucleotide, one for chromosome and one for strand.
	my $i = 0;

	while ($i < $length_maf_file)				
	{
		foreach my $species(@species_list)
		{
			$hash_genomic_cds{$species}{$i} = "?";
			$hash_base{$species}{$i} = "?";
			$hash_chr{$species}{$i} = "?";
			$hash_strand{$species}{$i} = "NA";		## Just initialize it to "NA" instead of "?" to avoid confusion with the minus strand
		}
		
		$i++;
	}
	
	open(FI,$mafFile);
	my @mafArray = <FI>;
	close FI;
	
	my $Ref2BlockLength = GetBlockLengthCum($mafFile,$reference);	## Get an array where each element is the length of the maf block (cumulative). Do not need a hash here, since the lengths can be easiy identified from their index in the array.
	my @BlockLength = @$Ref2BlockLength;							## If the length of 1st block is 10, the first element will be 0, the second will be 0+10 = 10....and so on.

	my $Ref2BBHash = GetBlockBoundaries($mafFile);			## This gives the reference to the hash where the key is the block, the value is a string of the form start-stop
	my %BlockBoundaries = %$Ref2BBHash;						## where start is the line number from where this block starts, stop is the line number where this block ends.

	my @blocks = keys(%BlockBoundaries);
	my @blocksList = sort {$a <=> $b} @blocks;				## Sort the blocks.
	
	foreach my $block(@blocksList)
	{
		my($start_block,$stop_block) = (split /-/,$BlockBoundaries{$block})[0,1]; ## Get the line number where the block starts and the block ends
		
		my @Block_Array = @mafArray[$start_block .. $stop_block];	## Get a smaller array (corresponds to the maf block) called Block_Array which is a sub-array of the bigger array (the entire maf file)
		pop @Block_Array;
		my  $nol = scalar(@mafArray);
		
		foreach my $line(@Block_Array)
		{
			if   ( ($line =~/^s/) && ($line ne "") )		### IF the lines are sequence lines (starting with s).
			{
				my $b_ct = $BlockLength[$block - 1];
				
				my $start_index = "";
				my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line);
				my @seqList = split(//,$seq);
				
				$start_index = $start;		## Read the sequence line sequentially. And populate the hashes accordingly.
					
				if ($strand eq "+")			## If the strand is +ve, the processing is relatively straight-forward.
				{
					foreach my $base(@seqList)
					{
						if ($base ne "-")
						{
							$hash_genomic_cds{$species}{$b_ct} = $start_index;
							$hash_base{$species}{$b_ct}        = $base;
							$hash_chr{$species}{$b_ct}         = $chr;
							$hash_strand{$species}{$b_ct}      = $strand;
							$start_index++;
						}
					
						$b_ct++;
					}
				}
				else					## If the strand is -ve, the processing is still relatively straight-forward, except that I need to "correct" the cds.
				{
					my $correct_cds = $srcSize - $start;
					$start_index	= $correct_cds;
					
					foreach my $base(@seqList)
					{
						if ($base ne "-")
						{
							$hash_genomic_cds{$species}{$b_ct}	= $start_index;
							$hash_base{$species}{$b_ct}			= $base;
							$hash_chr{$species}{$b_ct}			= $chr;
							$hash_strand{$species}{$b_ct}		= $strand;

							$start_index--;
						}

						$b_ct++;
					}
				}
			}						#### IF block code ends.
		}
	}
	
	foreach my $species(@species_list)
	{
		my $OutFile = "$outDir/$species.txt";
		open(FO,">$OutFile") || die "Error opening output file which is $OutFile\n";
		
		print "Writing to $OutFile\n";
		
		$i=0;
		while ($i < $length_maf_file)
		{		

			my $ref_g_cds    = $hash_genomic_cds{$reference}{$i};
			my $ref_g_base   = $hash_base{$reference}{$i};
			my $ref_g_chr    = $hash_chr{$reference}{$i};
			my $ref_g_strand = $hash_strand{$reference}{$i};
			
			
			my $status = my $codon_index = "-";
			if ( ($ref_g_cds ne "?") and ($ref_g_cds ne "-") )		## $ref_g_cds will be "?" if there is an insertion in a query sequence and this insertion is visible in the s line itself
			{														## $ref_g_cds will be "-" if there is an insertion in a query sequence and this insertion is captured from the FixILineModule
				if (exists $genome_cds_nature_hash{$ref_g_cds} )
				{
					$status = $genome_cds_nature_hash{$ref_g_cds};
					if ($status !~/Coding_exon/)	{	$codon_index = "-"; }
					else							{	$codon_index = $codon_index_hash{$ref_g_cds }; }
				}
			}
			
			else
			{
				$status = $codon_index = "-";
			}
			
			my $g_cds    = $hash_genomic_cds{$species}{$i};
			my $g_base   = $hash_base{$species}{$i};
			my $g_chr    = $hash_chr{$species}{$i};
			my $g_strand = $hash_strand{$species}{$i};
			
			print FO "$reference\t$status\t$codon_index\t$ref_g_strand\t$ref_g_cds\t$ref_g_base\t$ref_g_chr\t$species\t$g_strand\t$g_cds\t$g_base\t$g_chr\n";
			
			$i++;
		}
		close FO;
		
		PlaceInsertions($OutFile,$species,$Ref2HashInsertions)	if (exists $SpeciesWithInserts{$species});
		
		## Some formatting issues - Get rid of some useless lines, also place line numbers before each line.
		
		my $tmpFile1 = create_temp("f");
		my $tmpFile2 = create_temp("f");
		my $tmpFile3 = create_temp("f");
		
		system("awk '{ if ( (\$5 != \"-\") || (\$10 != \"-\") ) print \$0 }' $OutFile > $tmpFile1"); ## Get rid of those lines which are insertions in some other species and there is no sequence in either query or the Reference species
		system("awk '{ if ( (\$5 != \"?\") || (\$10 != \"?\") ) print \$0 }' $tmpFile1 > $tmpFile2"); ## Get rid of those lines which are missing sequence in both species (this would never happen with a reference species though)
		system("awk '{print NR \"\t\" \$0}' $tmpFile2 > $tmpFile3");	## Add line numbers at the beginning of each line in the file
		
		system("rm $tmpFile1 $tmpFile2");
		system("mv $tmpFile3 $OutFile");
		
		##  Formatting done
		
		FindDeletionsAndPlaceThem($mafFile,$strand_input,$species,$OutFile,$reference);	## Place Deletions
		FixMissingDataPositions($OutFile,$species,$ref2MissingData) if (exists $MissingDataSpecies{$species});	## Fix Missing Data Positions	
		ReverseFile($OutFile) if ($strand_input eq "-");	## Reverse file if strand is minus
		IncrementCds($OutFile);								## Fix CDS.
	}

}

sub GetBlockBoundaries ## Given a maf file, this function returns a hash. The key is the block number, and the value is a string of the form "Start-Stop" where Start is the line number where a maf block starts, Stop is the line number
{		       		   ## where the maf block ends.
	
	my $input = shift;

	my $lineCt = 1;
	my %BlockBoundaries = ();
	my $blockStart = my $blockStop = "";
	my $blockCt = 0;

	open(FI,$input) || die "Error opening $input in GetBlockBoundaries function \n";
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
	
		if ($line =~/score/)
		{
			$blockStart = $lineCt;
			$blockCt++;
		}
		elsif ($line eq "")
		{
			$blockStop = $lineCt;
			$BlockBoundaries{$blockCt} = "$blockStart-$blockStop";
		}
	
		$lineCt++;
	}
	close FI;

	return \%BlockBoundaries;
}

sub GetBlockLengthCum ### Gets the cumluative length of the blocks from a maf.
{
	my ($maf,$RefSpecies) = @_;
	my @lengthList = ();

	my $ct = 0;			## First get the length of each block and store these values in an array called @lengthList
	open(FI,$maf) || die "Error opening maf file ($maf) in GetBlockLengthCum\n";
	while (my $line = <FI>)
	{
		if ( ($line =~/^s/) && ($line =~/$RefSpecies/) )
		{
			$line =~s/\s+$//;
			my ($species,$chr,$start,$size,$strand,$srcSize,$seq) = getSLineData($line);

			$lengthList[$ct] = length($seq);
			$ct++;
		}
	}

	close FI;

	my $i  = my $sum = 0;	## Now get the cumulative length, the length of the first element is 0, the length of the second element is 0 + length of first element....and so on.
	my @BlockCumLength = ();

	while ($i < $ct)
	{
		if ($i == 0)
		{
			$BlockCumLength[$i] = 0;
		}

		else
		{	
			$BlockCumLength[$i] = $sum;
		}

		$sum = $sum + $lengthList[$i];       
		$i++;
	}
	
	return \@BlockCumLength;
}

sub GetMafLength	## This function gets the length of the maf.
{
	my ($reference,$tmp_maf) = @_;

	my $length_maf_file = 0;		## this variable contains the length of the entire maf block in the maf file, which includes all insertions, not just the "actual" sequence.
	open(FIM,$tmp_maf);
	while (my $line = <FIM>)
	{
		if ($line =~/$reference/)
		{
			my ($species,$chr,$start,$size,$strand,$srcSize,$seq) = getSLineData($line);
			$length_maf_file = $length_maf_file + length($seq);
		}
	}
	close FIM;
	return $length_maf_file;
}

sub PlaceInsertions	## This function helps you place the insert lines (which are captured from the i lines) in the tabular file
{
	my ($InFile,$species,$Ref2Hash) = @_;
	my %InsertHash = %$Ref2Hash;		## What you have here is a hash, the key is the position where the insert should be placed, the value is the a string of lines separated by new line character that are to be inserted.

	my @InsertPos = ();
	my $ip = 0;

	foreach my $genomeCds (keys %{$InsertHash{$species}})
	{
		$InsertPos[$ip] = $genomeCds;
		$ip++;
	}

	@InsertPos = sort {$a <=> $b} @InsertPos;

	foreach my $genomeCds(@InsertPos)
    	{
		my $string = $InsertHash{$species}{$genomeCds};		## The string of lines that has to be inserted.
		my @tmp = split(/\n/,$string);			

		my $StringFile = create_temp("f");		## This is important to make the two arrays, @tmp and @FileArray compatible
		open(FOSF,">$StringFile");
		foreach my $lineString(@tmp)
		{
			print FOSF "$lineString\n";
		}
		close FOSF;

		open(FOSF,$StringFile);
		@tmp = <FOSF>;
		close FOSF;					## Done, the two arrays will be compatible now.

		open(FileToArray,$InFile) || die "Error opening $InFile";
		my @FileArray = <FileToArray>;
		close FileToArray;

		splice (@FileArray,$genomeCds,0,@tmp);					## This is it basically. Add @tmp to @FileArray at a specific position
		
		open(FO_ATF,">$InFile") or die "error updating file\n";				## Update the file now
		foreach my $element(@FileArray)
		{
			print FO_ATF $element;
		}
		close FO_ATF;
		system("rm $StringFile");
	}
}

sub FindDeletionsAndPlaceThem	## This function helps you find deletions and place them appropriately in the tabular file. As a result, all the "?"s are replaced with "-"s at such positions.
{
	my ($mafIn,$StrandRef,$speciesCheck,$InFile,$RefSpecies) = @_;
	my %DelHash = ();
	
	my $species = my $chr = my $start = my $size = my $strand = my $srcSize = my $seq = "";
	my $delCt = 0;
	
	open(FIM,"$mafIn") || die "Error opening maf file\n";
	while (my $line = <FIM>)
	{
		$line=~s/\s+$//;
		
		if ( ($line =~/^s/) and ($line =~/$RefSpecies/) )
		{		
			($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line); 
		}
		
		if ($line =~/^e/) 			## This is the check point, checks if the line is an "e" line.
		{
			my ($speciesE,$chrE,$startE,$sizeE,$strandE,$srcSizeE,$statusE) = getELineData($line);
			
			if ($speciesE eq $speciesCheck)
			{
			
				my $begin = my $end = "";			## If you have nothing but an "e" line for a species, it implies that there is a deletion in the query species
													## Fix the positions of "begin" and "end" first.
				if ($StrandRef eq "+")
				{
					$begin = $start;
					$end   = $start + $size;
				}
				else
				{
					$begin = $start;
					$end   = $begin + $size;
				}
			
				while ($begin < $end)			## And put everything in a hash called DelHash
				{									## The key is the genomic coordinate in the reference, the value is "-"
					$delCt++;
					$DelHash{$begin} = "-";
					$begin++;
				}
			}
		}
		
		elsif ($line =~/^s/) 	## Now this is a bit tricky, here I do not have a clean deletion but instead some parts of the sequence are deleted.
		{			## Has to be dealt carefully.
			
			my ($speciesS, $dm2, $dm3, $dm4, $dm5, $dm6, $seqQuery) = getSLineData($line); 	## dm means dummy, I do not care about these variables. All I need is the sequence for the query species.
			if ($speciesS eq $speciesCheck)
			{
				my @SeqQuery = split(//,$seqQuery);
			
				my @SeqRef = split(//,$seq);

				my $begin = $start;
			
				for (my $j =0; $j <@SeqRef; $j++)
				{
					if ( ($SeqQuery[$j] eq "-") && ($SeqRef[$j] ne "-")  )	## If the sequence in query is "-" but it is not the same in reference, then it is a deletion
					{														## And hence, put these coordinates in DelHash.	
						$delCt++;
						$DelHash{$begin} = "-";
						$begin++;	
					}
					elsif (($SeqRef[$j] ne "-"))							## else do not do anything, just incremement the value of begin.
					{
						$begin++;
					}
				}
			}
		}
	}
	close FIM;
	
	## At this point, I substitute the "?" with "-" because by now, I have identified the deletions.
	if ($delCt > 0)				## Do this only if there are some deletions in the query species
	{							
		my $FileTmp = create_temp("f");
	
		open(FID,"$InFile");
		open(FOD,">$FileTmp");
	
		while (my $line = <FID>)
		{
			$line =~s/\s+$//;
			my @tmp = split(/\t/,$line);
	
			if (exists $DelHash{$tmp[5]})
			{
				$tmp[10] = $tmp[11] = $tmp[12] = "-";
				$tmp[9] = "NA";
				$line = join("\t",@tmp);
			}
		
			print FOD "$line\n";
		}
	
		close FID;
		close FOD;
		system("mv $FileTmp $InFile");
	}
	
	## Now fix insertions in the query sequence.
	
	my $FileTmp = create_temp("f");
	
	open(FID,"$InFile");
	open(FOD,">$FileTmp");
	
	while (my $line = <FID>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
	
		if ( ($tmp[10] ne "-" ) and ($tmp[5] eq "?") )	## This implies that if the sequence in the query is not "-" but the sequence in reference is ?, it is actually a deletion in the Reference sequence.
		{
			$tmp[2] = $tmp[3] = $tmp[5] = $tmp[6] = $tmp[7] = "-";
			$tmp[4] = "NA";
			$line = join("\t",@tmp);
		}
		print FOD "$line\n";
	}
	
	close FID;
	close FOD;
	system("mv $FileTmp $InFile");
}

sub FixMissingData		## This function helps you identify missing data in the tabular file, where it could be mistaken for deletion.
{
	my ($maf,$Reference) = @_;
	
	open(FI,$maf) || die "Error opening $maf input\n";
	my @mafArray = <FI>;
	close FI;
	
	my $RefCds = my $lengthBlockRef = "";
	my %HashMissingData = ();	## Everything is stored in a 2 dimensional hash called HashMissingData --> the first dimension is the species, 
								## the second is the genomic-Coordinate, the value is species
								
	for(my $i = 0; $i < @mafArray; $i++)
	{
		my $line = $mafArray[$i];
		if ($line =~/$Reference/)	## Get those lines which contain the reference species and get RefCds and length of the block from them.	
		{
			my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line);
			$RefCds 	= $start;
			$lengthBlockRef = $size;
		}
		
		elsif ($line =~/^s/)		## Get those lines which are s lines but not from the reference
		{
			my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($line);
			
			my $iLine = $mafArray[$i+1];	## Get the i line for this s line.
			$iLine =~s/\s+$//;
			
			if ($iLine =~/^i/)	## check if this is indeed an i" line. For evolver-generated mafs
			{			
				my($speciesDM, $chrDM, $statusUp, $countUp, $statusDown, $countDown) = getILineData($iLine); 
				
				if ( ($statusUp eq "M") && ($seq =~/^-/) )	## Check if the statusUp is "M" and sequence begins with "-"
				{
					my @seqArray = split(//,$seq);
					my $cdsStart = $RefCds;

					foreach my $base(@seqArray)			## Take this sequence and replace all those positions starting from the beginning until you find a non"-" character
					{									## and replace with N
						if ($base eq "-")
						{
							$HashMissingData{$species}{$cdsStart} = "N";
							$cdsStart++;
						}
						else
						{
							last;
						}
					}
				}
			}
		}
		
		elsif ($line =~/^e/)	### Also fix those e-lines which are actually Missing Data and not entire deletions
		{						### Remember all e-lines are not deletions.
			
			my ($species,$chr,$start,$size,$strand,$srcSize,$status) = getELineData($line);
			
			if ($status eq "M")		## An example:	"e gorGor1.Supercontig_0015369        14840 92 +     70208 M"		## The M at the end, i.e the status tells us that this is missing data
			{
				my $start = $RefCds;
				my $stop  = $RefCds + $lengthBlockRef -1;
				
				while ($start <= $stop)
				{
					$HashMissingData{$species}{$start} = "N";
					$start++;
				}
			}
		}
	}

	return \%HashMissingData;
}

sub FixMissingDataPositions		## The missing data positions that are identified at the previous step, are now substituted with "?" in place of "-".
{
	my ($fileIn,$sp,$Ref2Hash) = @_;
	my %HashMissingDataPos = %$Ref2Hash;
	
	my $fileOut = create_temp("f");
	
	open(FI,"$fileIn");		## Simply read a tabular file which is created at the previous step.
	open(FO,">$fileOut");	
	
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		my @lineTmp = split(/\t/,$line);
		
		my $RefCds = $lineTmp[5];	## Get the genomic coordinate of the Reference species
									## If this value occurs in the hash, then it is 
		if (exists $HashMissingDataPos{$sp}{$RefCds})
		{
			print FO "$lineTmp[0]\t$lineTmp[1]\t$lineTmp[2]\t$lineTmp[3]\t$lineTmp[4]\t$lineTmp[5]\t$lineTmp[6]\t$lineTmp[7]\t$lineTmp[8]\tNA\t?\t?\t?\n";		## Make changes to the last 4 fields.
		}			## Note that the browswer will show these positions as Ns.
		else
		{
			print FO "$line\n";
		}
	}
	
	close FI;
	close FO;
	
	system("mv $fileOut $fileIn");
}

sub ReverseFile	## For a gene on minus strand, the file has to be reversed. Read from bottom to top, get the reverse complement of the sequence, also invert the strands
{
	my $fileIn = shift;
	
	my $fileOut = create_temp("f");
	my $fileInverted = create_temp("f");
	
	my $statusTAC = system("tac $fileIn > $fileInverted"); 	## tac is the linux utility that helps you read a file from bottom to top
	die "Cannot invert file $fileIn\n" if ($statusTAC != 0);
	
	open(FI,$fileInverted);
	open(FO,">$fileOut");
	
	#my %hashNT = ("A" => "Y", "C" => "Y", "G" => "Y", "T" => "Y");
	my @ntChars = qw(a A t T c C g G n N);
	my %hashNT = map{$_ => 1}@ntChars;
	
	my $i = 1;
	
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);
		
		shift @tmp;
		## Invert strands for both the reference and the query, if the strand value is not "NA"
		$tmp[3] = InvertStrand($tmp[3])	if ($tmp[3] ne "NA");		
        	$tmp[8] = InvertStrand($tmp[8]) if ($tmp[8] ne "NA");

		## Get the reverse-complement for both the reference and the query, if the nucleotide at the position is one of the 4 standard nts.
		$tmp[5]  = revComp($tmp[5])	if (exists $hashNT{$tmp[5]});
		$tmp[10] = revComp($tmp[10])	if (exists $hashNT{$tmp[10]});
		
		$line = join("\t",@tmp);
		
		print FO "$i\t$line\n";
		$i++;
	}
	
	close FI;
	close FO;
	
	system("mv $fileOut $fileIn");
	system("rm $fileInverted");
}

sub InvertStrand
{
	my $in = shift;
	my $out = "";

	$out = "+" if ($in eq "-");
	$out = "-" if ($in eq "+");
	return $out; 
}

sub IncrementCds
{
	my $fileIn = shift;
	my $fileOut = create_temp("f");

	## There are two rules here. The Reference Genome Cds has to be always incremented
	## The Query Genome Cds --> well, depends. If the strand is same for both the query genome and the reference genome, then this increment is needed, otherwise not.

	## Get the first line where there is sequence for both the genomes, get the strands and make a decision on whether you want to increment the GenomeCds for the query or not.
	my $lineCheck = `awk '{ if ( (\$5 != "NA") && (\$10 != "NA") ) print \$0 }' $fileIn |head -n 1`;
    my $flagIncrementQueryCds = "F";

	if ($lineCheck ne "")
	{
	    $lineCheck =~s/\s+$//;
		my @tmp = split(/\t/,$lineCheck);

		my $strandReference = $tmp[4];
		my $strandQuery	    = $tmp[9];
	    $flagIncrementQueryCds = "T" if ($strandReference eq $strandQuery);
	}
	
	## Set flagIncrementQueryCds to True, if the strand for both the query and the reference species is the same.
	
	open(FI,$fileIn);
	open(FO,">$fileOut");
	while (my $line = <FI>)
	{
		$line =~s/\s+$//;
		my @tmp = split(/\t/,$line);

		$tmp[5]++ if ($tmp[4] ne "NA");
		$tmp[10]++ if ( ($tmp[9] ne "NA") && ($flagIncrementQueryCds eq "T") );

		my $lineNew = join("\t",@tmp);
		print FO "$lineNew\n";
	}

	close FI;
	close FO;

	system("mv $fileOut $fileIn");
}
