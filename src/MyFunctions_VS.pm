#!/sw/bin/perl

### Virag Sharma, July 2013.
### This module contains a set of functions for parsing the pairwise tabular files that I create for the gene loss pipeline.
### Dependencies -> Uses functions from MH's MyFunctions.pm, MyMaxEnt.pm and VS's useful_functions_VS.pm

package MyFunctions_VS;
use strict;
use warnings;
use Exporter;

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyKentFunctions;
use MyFunctions;
use MyMaxEnt;

use lib "$ENV{'GeneLossPipeCode'}";
use useful_functions_VS;

our @ISA = ('Exporter');
our @EXPORT = qw(fetch_seq  get_exonic_seq get_anchoring_cds  get_chr_strand get_noe get_exonic_chunk get_sequence_context_Stop_Codon get_chr_length);

our $genomePath = $ENV{'genomePath'};

#################################################################################################
####     Get the nucleotide sequence for an assembly given the genomic coordinate 	 		#####
####     (the start point), the length, the strand, the species and the chromosome 		   	#####
#################################################################################################

sub fetch_seq
{					
	my($start,$stop,$strand,$species,$chromosome) = @_;
	
	my $seqNU = my $seqND = my $flag = ""; ## $seqNU --> seq of Ns upstream, $seqND --> seq of Ns downstream
    	my $chromSize = get_chr_length($chromosome,$species);
	
	if ($start >= $chromSize)
	{
		my $seq = GenerateString(abs($start - $stop));
		return ($seq,$flag);
	}

	if ($start < 0) 
	{
		$seqNU  = GenerateString(abs($start));		### Pad something on the upstream
		$start = 0;
	}

	if ($stop > $chromSize)		## Padding with Ns to extend the size scaffold does not really work. That is because I end up with coordinates that are not real. And if we were to run TwoBitToFa again on these genomic coordinates, then it will fail.
	{
		$stop = $chromSize;
	}
	
	die "Incorrect parameters to twoBitToFa while fetching sequence\n" if ($start < 0 || $stop < 0 || $start > $stop);
	
	my $sp2bit = "$genomePath/gbdb-HL/$species/$species.quality.2bit";
	$sp2bit =  "$genomePath/gbdb-HL/$species/$species.2bit" if (! -e $sp2bit);

	my $seq = `twoBitToFa $sp2bit:$chromosome:$start-$stop stdout|grep -v ">"|tr -d "\n"`;
	die("Cannot get the sequence from 2BitToFa") if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	$seq = revComp($seq) if ($strand eq "-");	

	$seq =~tr/atcgn/NNNNN/;		### This is important, since we are extracting the sequence from quality.2bit file...there will be some lower case characters. Substitute those characters with N.
	
	$seq = $seqNU.$seq.$seqND;
	return $seq;
}

################################################################
####	Extract the exonic sequence from a pairwise file	####
################################################################
sub get_exonic_seq
{
	
	my ($acc,$species,$exonIndex,$pep_flag,$refSpecies,$data_path,$RefFlag) = @_;				### $pep_flag ==><== { "E" => "Return Exon", "P" => "Return Protein", "B" => "Return Both" };
	
	die "Argument not defined for $pep_flag in get_exonic_seq function\n" if ( ($pep_flag ne "E") and ($pep_flag ne "P") and ($pep_flag ne "B") );
	my $file = "$data_path/$acc/$species.txt";
	
	if (-e $file)
	{	
		open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
		my @fileArray = <FIR>;
		close FIR;  
		
		my $refFile = "$data_path/$acc/$refSpecies.txt";
		my $exonic_term = "Coding_exon\_$exonIndex";
		
		my $tmp_result = `grep -w "$exonic_term" $refFile`;
		die "Improper exon index $exonIndex. This exon is not present in the file $refFile \n" if ($tmp_result eq "");
		
		my $term1 = "$exonIndex\_A_1";
		my $term2 = "$exonIndex\_D_1";
		
		## Get the number of exons
		my $noe = get_noe($acc,$refSpecies,$data_path);		
		$term1 = "Coding_exon_1" if ($exonIndex == 1);
		$term2 = "Coding_exon\_$noe" if ($exonIndex == $noe);
		
		## Get start and stop positions:
		my ($startPt,$stopPt) = getStartStopPos($file,$term1,$term2);
		my $exonic_seq = my $exonic_seq_pep = "";
		
		if ($RefFlag eq "T"){	
			while ($startPt <= $stopPt)
			{
				my $line = $fileArray[$startPt];
				my($nt,$pep) = (split /\t/,$line)[7,8];
				$exonic_seq = $exonic_seq.$nt;
				$exonic_seq_pep = $exonic_seq_pep.$pep;
				
				$startPt++;
			}
		}else{	
			while ($startPt <= $stopPt)
			{
				my $line = $fileArray[$startPt];
				my($nt,$pep) = (split /\t/,$line)[13,14];
				$exonic_seq = $exonic_seq.$nt;
				$exonic_seq_pep = $exonic_seq_pep.$pep;
				
				$startPt++;
			}
		}

		if ($pep_flag eq "E"){	
			return $exonic_seq;
		}elsif ($pep_flag eq "P"){	
			return $exonic_seq_pep;
		}elsif ($pep_flag eq "B"){
			return ($exonic_seq,$exonic_seq_pep);
		}
	}
	
	else
	{
		die "The gene $acc is not present in the species $species, the file '$file' is missing\n";
	}
}

################################################################################################################################################
####     Extract the anchoring coordinate- the point from where you start to extract the intronic sequence	which means that you extract	####
####	 the position of the last base that aligns to the human exon (when looking for splice site in downstream intron => 	donor  			####
####	 site here) or the first base that aligns to the human exon (when looking for splice site in upstream intron => acceptor site) 		####			
################################################################################################################################################

sub get_anchoring_cds
{
	my ($acc,$species,$exon_index,$type,$RefSpecies,$data_path) = @_;			## $type is either Donor (D) or Acceptor (A).

	my $file = "$data_path/$acc/$species.txt";
	
	my $Ref_file = "$data_path/$acc/$RefSpecies.txt";
	my $exonic_term = "Coding_exon\_$exon_index";
	
	my $tmp_result = `grep -w "$exonic_term" $Ref_file`;
	die "Improper exon index $exon_index. This exon is not present in the file $Ref_file\n" if ($tmp_result eq "");
	
	my $term_1 = "$exon_index\_A_1";
	my $term_2 = "$exon_index\_D_1";
	my $nol    = `wc -l $file|awk '{print \$1}'|tr -d "\n"`;
	
	my $line_start = my $line_stop = "";
	$line_start = `grep -n -w "$term_1" $file|awk '{split(\$0,a,":"); print a[1]}'|tr -d "\n"`;
	
	if ($line_start eq "") { $line_start = `grep -n "Coding_exon" $Ref_file|head -n 1|awk '{ split(\$0,a,":"); print a[1]}'|tr -d "\n"`; }
	else 		       { $line_start++;   }
	
	$line_stop = `grep -n -w "$term_2" $file|awk '{split(\$0,a,":"); print a[1]}'|tr -d "\n"`;
	
	if ($line_stop eq "") { $line_stop = `grep -n "Coding_exon" $Ref_file|tail -n 1|awk '{ split(\$0,a,":"); print a[1]}'|tr -d "\n"`; }
	else { $line_stop--; }
	
	my $pos = "";
	
	$pos = `sed -n $line_start,$line_stop\\p $file|awk '{if ( (\$14 != "-") && (\$14 != "?") ) print \$11}'|tail -n 1|tr -d "\n"` if ($type eq "D");
	$pos = `sed -n $line_start,$line_stop\\p $file|awk '{if ( (\$14 != "-") && (\$14 != "?") ) print \$11}'|head -n 1|tr -d "\n"` if ($type eq "A");
	
	return $pos;
}


####################################################################################################################
####		Get the strand and the chromosome, given the accession, the species and the exon index		####
####################################################################################################################

sub get_chr_strand
{
	my($accession,$species,$exon_index,$data_path) = @_;
	
	my $spFile   = "$data_path/$accession/$species.txt";
	my $grep_term = "Coding_exon\_$exon_index";
	
	my $chromStart = my $chromStop = my $strandStart = my $strandStop = "";
	
	open(FIR,$spFile) || die "Error opening species file in '$spFile' in get_chr_strand function\n";
	while (my $line = <FIR>)
	{
		my($exonTerm,$chr,$strand,$nt) = (split /\t/,$line)[5,11,12,13];
		
		if ($grep_term eq $exonTerm && $nt ne "-" && $nt ne "?")
		{
			$chromStart = $chr;
			$strandStart = $strand;
			last;	
		}
	}
	close FIR;
	
	open(FIR,$spFile) || die "Error opening species file in '$spFile' in get_chr_strand function\n";
	while (my $line = <FIR>)
	{
		my($exonTerm,$chr,$strand,$nt) = (split /\t/,$line)[5,11,12,13];
		
		if ($grep_term eq $exonTerm && $nt ne "-" && $nt ne "?")
		{
			$chromStart = $chr;
			$strandStart = $strand;
		}
	}
	close FIR;
	
	if ($chromStart ne $chromStop || $strandStart ne $strandStop){
		return("NA","NA");
	}else{
		return($chromStart,$strandStart);
	}
}


############################################################################
####	Given an accession, get the number of exons in this transcript	####
############################################################################

sub get_noe
{
	my ($accession,$refSpecies,$data_path) = @_;

	my $noe = "";
	open(FIR,"$data_path/$accession/$refSpecies.txt") || die "Error opening input file in get_noe function for '$accession' in '$data_path'\n";
	while(my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my $nature = (split /\t/,$line)[5];
		$noe = $nature if ($nature =~/Coding_exon/);
	}
	close FIR;
	$noe =~s/Coding_exon_//;
	
	return $noe;
}


############################################################################
####	Same as get_exonic_seq. Just write all other fields as well.	####
############################################################################

sub get_exonic_chunk
{
	my ($acc,$species,$exonIndex,$refSpecies,$data_path) = @_;				
	
	my $file = "$data_path/$acc/$species.txt";
		
	if (-e $file){	
		open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
		my @fileArray = <FIR>;
		close FIR;  
		
		my $term1 = "$exonIndex\_A_1";
		my $term2 = "$exonIndex\_D_1";
		
		## Get the number of exons
		my $noe = get_noe($acc,$refSpecies,$data_path);
		$term1 = "Coding_exon_1" if ($exonIndex == 1);
		$term2 = "Coding_exon\_$noe" if ($exonIndex == $noe);
		
		## Get start and stop positions:
		my ($startPt,$stopPt) = getStartStopPos($file,$term1,$term2);
		
		my @exonicSeqArray = @fileArray[$startPt..$stopPt];
		my $exonicSeq = join("",@exonicSeqArray);
		
		return $exonicSeq;
	}else{
		die "The gene $acc is not present in the species $species\n";
	}
}

########################################################################################################
####		Given the species and the chromosome number, find the size of this sequence				####
########################################################################################################

sub get_chr_length		## basically grep for the input chromosome in the relevant chrom.sizes file
{
	my ($input_chr,$species) = @_;
	
	my $sp_chr_sizes = "$genomePath/gbdb-HL/$species/chrom.sizes";
	my $size = `grep -w "$input_chr" $sp_chr_sizes|awk '{print \$2}'|tr -d "\n"`;
	
	return $size;	
}

sub GenerateString
{
	my $length = shift;

	my $i = 1;
	my $string = "";

	while ($i <= $length)
	{
		$string = $string."N";
		$i++;
	}

	return $string;
}	

## NOT exported, used internally.

sub getStartStopPos
{
	my($file,$term1,$term2) = @_;
	
	die "The input file '$file' does not exist (Dying in getStartStopPos function)\n" if (! -e $file);
	
	## Get the line start and line stops:
	my $ct = 0;
	my $startPt = my $stopPt = "";
	open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
	while (my $line = <FIR>)
	{
		my $nature = (split /\t/,$line)[5];
		
		$startPt = $ct if ($nature eq $term1);
		$stopPt  = $ct if ($nature eq $term2);
		$ct++;
	}
	close FIR;
	
	##  Exception for first exon
	if ($term1 eq "Coding_exon_1")
	{
		my $ct = 0;
		open(FIR,$file) || die "Error opening species file '$file' in get_exonic_seq function\n";
		while (my $line = <FIR>)
		{
			my $nature = (split /\t/,$line)[5];
			if ($nature eq "Coding_exon_1")
			{
				$startPt = $ct;
				last;
			}
			$ct++;
		}
		close FIR;
	}
	
	$startPt++ if ($term1 =~/_A_/);
	$stopPt-- if ($term2 =~/_D_/);
	return ($startPt,$stopPt)	;	
}

1;
