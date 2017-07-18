#!/sw/bin/perl

### Virag Sharma, November 2013.
### This module contains a set of functions associated with the Splice Site Module in the MutatedSS.pl pipeline
### Dependencies -> Uses functions from MH's MyFunctions.pm, MyMaxEnt.pm and VS's useful_functions_VS.pm

package SpliceSiteFunctionsVS;
use strict;
use warnings;
use Exporter;

use lib "$ENV{'genomePath'}/src/LabPerlModules"; 
use MyFunctions;
use MyMaxEnt;

use lib "$ENV{'GeneLossPipeCode'}";
use MyFunctions_VS;
use useful_functions_VS;

our @ISA = ('Exporter');
our @EXPORT = qw(GetJunctionList GetSSList_Ref GetSpliceSite GetSSCds get_1nt);

our $genomePath = $ENV{'genomePath'};
our $dirEST = $ENV{'data_DIR'};

################################################################################################################
####	This function generates a list of search terms (where each term is associated with an exon/intron   ####
####    junction), the input is the number of exons in the gene. For instance, if the number of exons in a  ####
####	gene is 4, the list will have terms: 1_D,2_A,2_D,3_A,3_D and 4_A									####
################################################################################################################

sub GetJunctionList						
{									
	my $noe = shift;  ## The number of exons in the gene
	die "No junctions possible in a gene with only $noe exons\n" if ($noe < 2);
	
	my $startPt = 2;
	my $stopPt  = $noe - 1;
	
	my @array_list = ();
	push(@array_list,"1\_D");
		
	if ($startPt <= $stopPt)
	{	
		while ($startPt <= $stopPt)
		{
			push(@array_list,"$startPt\_A");
			push(@array_list,"$startPt\_D");
			
			$startPt++;
		}
	}
	
	push(@array_list,"$noe\_A"); ## The last element of the array
	return(\@array_list);
}

################################################################################################################
####				Given an accession, get a list of its splice sites for the reference gene				####
################################################################################################################

sub GetSSList_Ref
{
	my($accession,$junctionList,$reference,$data_path) = @_;
	my %ssList = my %ssNucleotides = ();
	
	my @spliceSitesList = @$junctionList;
	my $refFile = "$data_path/$accession/$reference.txt";		
	
	die "The reference file i.e $refFile does not exist\n" if (! -e $refFile);
	
	open(FIR,$refFile) || die "Error opening $refFile for gene '$accession', reference is '$reference'\n";
	while (my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my ($ssPos,$nt) = (split /\t/,$line)[5,13];
		
		$ssNucleotides{$ssPos} = $nt if ($ssPos =~/D_1$/ || $ssPos =~/D_2$/ || $ssPos =~/A_1$/ || $ssPos =~/A_2$/);
	}
	close FIR;
		
	foreach my $junction(@spliceSitesList)
	{
		my $exon_index = $junction;
		$exon_index=~s/[_AD]//g;
		
		my $junc1 = "$junction\_1";
		my $junc2 = "$junction\_2";
		
		my $nt1 = $ssNucleotides{$junc1};
		my $nt2 = $ssNucleotides{$junc2};
		
		my $splice_site = $nt1.$nt2;
		$splice_site    = reverse($splice_site) if ($junction =~/_A/);
		
		$ssList{$junction} = $splice_site;
	}
	return(\%ssList);
}

################################################################################################################
####	This function returns the splice site and the information for the human splice site (useful for     ####
####    creating hyperlink): given the junction (like 1_A, 2_D etc) and the SIQ/human tabular file.	    ####
#### 	Also returns the splice site coordinates of the query species					    ####			
################################################################################################################


sub GetSpliceSite				
{
	my ($junction,$file) = @_;
	
	die "The file $file does not exist. Dying in GetSpliceSite function\n" if (! -e $file);

	my $junc1 = "$junction\_1";  ## The junction is an exon plus the splice site. For example 5_A means 5th exon, acceptor site
	my $junc2 = "$junction\_2";
	
	my %fileToHash = ();
	open(FIR,"$file") || die "Error opening '$file' in GetSpliceSite function\n";
	while (my $line = <FIR>)
	{
		chomp $line;
		my $nature = (split /\t/,$line)[5];
		$fileToHash{$nature} = $line if ($nature ne "-");
	}
	close FIR;
	
	my $nt1 = (split /\t/,$fileToHash{$junc1})[13];
	my $nt2 = (split /\t/,$fileToHash{$junc2})[13];
	
	my $refCds = my $refStr = my $qSpeciesCds = "";
	($refCds,$refStr,$qSpeciesCds) = (split /\t/,$fileToHash{$junc1})[2,4,10];

	if ($junction =~/A/)
	{
		my($dm1,$dm2,$dm3) = (split /\t/,$fileToHash{$junc1})[2,4,10];
		$qSpeciesCds = $dm3;
	}
	
	my $spliceSite = $nt1.$nt2;
	$spliceSite = reverse($spliceSite) if ($junction =~/_A/);
	
	return ($spliceSite,$refCds,$qSpeciesCds,$refStr);
}

#################################################################################################################
####	Get a list of splice site genomic coordinates given the accession, the junction list and the species ####
#### 	The key is the junction, the value is the string where the different features are junction type,	 ####
####	Genomic cds, the strand and the species																 ####
#################################################################################################################

sub GetSSCds
{
	my($accession,$junctionList,$species,$data_path) = @_;
	my %ssCds = ();
	
	my @spliceSitesList = @$junctionList;
	my $speciesFile = "$data_path/$accession/$species.txt";		
	
	die "The species file i.e $speciesFile does not exist\n" if (! -e $speciesFile);
	
	my %chrHash = my %strandHash = my %cdsHash = ();
	open(FIR,$speciesFile) || die "Error opening $speciesFile for gene '$accession', species is '$species'\n";
	while (my $line = <FIR>)
	{
		$line =~s/\s+$//;
		my ($ssPos,$nt,$chr,$strand) = (split /\t/,$line)[5,10,11,12];
		
		if ($ssPos =~/D_1$/ || $ssPos =~/D_2$/ || $ssPos =~/A_1$/ || $ssPos =~/A_2$/)
		{
			$cdsHash{$ssPos} 	= $nt;
			$strandHash{$ssPos} = $strand;
			$chrHash{$ssPos} 	= $chr;
		}
	}
	close FIR;
	
	foreach my $junction(@spliceSitesList)
	{
		my $exon_index = $junction;
		my @tmp = split(/_/,$exon_index);
		my $junction_type = $tmp[1];
		
		my $junc = "";
		$junc = "$junction\_1" if ($junction_type eq "D");
		$junc = "$junction\_2" if ($junction_type eq "A");
	
		my $gCds 	= $cdsHash{$junc};
		my $chr  	= $chrHash{$junc};
		my $strand  = $strandHash{$junc};
	
		my $string = "$junction_type<>$gCds<>$chr<>$strand<>$species";
		$ssCds{$junction}=$string;
	}
	
	return(\%ssCds);
}

###########################################################################################
###	  EXPORTED. Find one nucleotide downstream/upstream of a given genomic coordinate	###	
###   Requires species, acc, junction, reference species and dataDir  					###
###########################################################################################

sub get_1nt
{
	my ($species,$acc,$junction,$Reference,$data_path)=@_;
	
	my @tmp_junction = split(/_/,$junction);
	my $exon_index	 = $tmp_junction[0];
	my $type		 = $tmp_junction[1];
	
	my ($chromosome,$strand) = get_chr_strand($acc,$species,$exon_index,$data_path); 		   ## Get the chromosome and the strand (from SIQ) that aligns to the human exon.
	my $cds = get_anchoring_cds($acc,$species,$exon_index,$type,$Reference,$data_path); ## Get anchoring cds i.e. the genomic coordinate from where to start while extracting the intronic sequence.

	my $seq = "";
	
	if ( ($chromosome ne "NA") and ($cds ne "") )
	{
	
		my $start_cds = $cds;
		my $stop_cds  = "";

		if ($type eq "D")
		{
			if ($strand eq "+") {	$start_cds = $cds + 2; $stop_cds = $start_cds + 1;	}	  ## Fixed
			else {	$start_cds = $cds -3; $stop_cds = $start_cds + 1;	}           		  ## Fixed
		}
	
		else
		{
			if ($strand eq "+") {	$start_cds = $cds - 2; $stop_cds = $start_cds + 1;	}	  ## Fixed
			else {	$start_cds = $cds + 1; $stop_cds = $start_cds + 1;	}					  ## Fixed
		}

		$seq = fetch_seq($start_cds,$stop_cds,$strand,$species,$chromosome);	
	}
	else
	{
		$seq = "-";
	}
	
	return ($seq);
}

###############################################
###	  NOT EXPORTED. Get offset position.	###
###############################################

sub get_offset
{
	my $in = shift;
	my %hash_pairs =("1_2" => 2, "1_3" => 1, "2_1" => 1, "2_3" => 2, "3_1" => 2, "3_2" => 1);
	
	my $offset = $hash_pairs{$in};
	return $offset;
}

1;
