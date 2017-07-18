#!/sw/bin/perl


## Virag Sharma, 2017. MPI-CBG and MPI-PKS, Dresden, Germany

## This script takes as an input, a file in my GenePrediction format
## It outputs a bed file and nature file (described later)
## The bed and the nature file are the first inputs for my EPTM_script

use strict;
use warnings;

#### Minus strand	BAD	NM_032989	chr11	-	64037299	64052176	64037680	64051840	64037680-64037809,64039084-64039275,64051653-64051840
#### Plus strand	ENPP4	NM_014936	chr6	+	46097700	46114436	46107320	46111377	46107320-46108146,46108788-46108959,46111012-46111377

use lib "$ENV{'genomePath'}/src/LabPerlModules";
use MyFunctions;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

sub usage
{
	die "usage: perl $0 -in Input_file -out Output_Directory (if nothing is specified, then the files are written to the present working directory) -length (Length of the flanks around the splice sites, minimum should be 2 so that you can extract the splice site)\n";
}

my $file = my $intronic_flag = my $out_dir = "";
our $lengthRegion = "";  ## set this as a global variable instead of passing it to all the subroutines that might need this value.

GetOptions ("in=s" =>\$file, "out:s" =>\$out_dir, "length=i" =>\$lengthRegion,"intron_flag:s" =>\$intronic_flag) || usage();
usage if ($file eq "");

if ($out_dir eq "")									### There has to be an input file in myGenePrediction format
{											### If no output directory is specified, the script writes to the current directory.
	$out_dir = `pwd|tr -d "\n"`;				
}
else
{
	`mkdir -p $out_dir`;
}

open(FI,"$file") || die "Error opening input file\n";
while (my $l=<FI>)
{
	$l =~s/\s+$//;
	
	my ($acc,$chr,$strand,$exons) = (split /\t/,$l)[1,2,3,8];
	
	## Do some sanity checks for the input file
	die "The strand value $strand is not defined\n" if ( ($strand ne "+") && ($strand ne "-") );	## Die if the strand is neither + nor -.
	#print "WARNING: The chromosome value '$chr' does not start with 'chr'\n" if ($chr !~/^chr/);				## All chromosome values start with "chr". For example, chr7_gl000195_random, chr1, chrUn_gl000211
	
	my $natureFile = GenerateNatureFiles($strand,$exons,$lengthRegion);
	my $bedFile	   = GenerateBedFiles($exons,$lengthRegion,$chr);
	
	### Add -20/-10 or whatever your desired length is to the first element in the bed file and add the same number to the last element in the bed file. 
	## This is the same for both the strands.
	
	my $tmpBed = create_temp("f");
	open(FIB,$bedFile);
	my @bedFileArray = <FIB>;
	close FIB;
	
	my $firstElementBedFile = $bedFileArray[0];
	my $lastElementBedFile  = $bedFileArray[$#bedFileArray];
	
	$firstElementBedFile =~s/\s+$//;
	$lastElementBedFile =~s/\s+$//;
	
	my @tmpFirst = split(/\t/,$firstElementBedFile);
	my $start = $tmpFirst[1];
	my $newElementF = $tmpFirst[1] - $lengthRegion;
	my $firstElement = "$tmpFirst[0]\t$newElementF\t$tmpFirst[2]";
	
	my @tmpLast = split(/\t/,$lastElementBedFile);
	my $stop = $tmpLast[2];
	my $newElementL = $tmpLast[2] + $lengthRegion;
	my $lastElement = "$tmpLast[0]\t$tmpLast[1]\t$newElementL";
	
	if ($firstElementBedFile eq  $lastElementBedFile)		## This would hold true for single-exon genes.
	{
		open(FOB,">$tmpBed");
		print FOB "$tmpFirst[0]\t$newElementF\t$newElementL\n";
		close FOB;
	}
	
	else
	{
		$bedFileArray[0] 	      = $firstElement;
		$bedFileArray[$#bedFileArray] = $lastElement;
	
		open(FOB,">$tmpBed");
		foreach my $lineBed(@bedFileArray)
		{
			$lineBed =~s/\s+$//;
			print FOB "$lineBed\n";
		}
		close FOB;
	}
	
	`mv $tmpBed $bedFile`;
	
	## Now fix the same thing in the nature files
	
	my $fileUp 	 = create_temp("f");
	my $fileDown = create_temp("f");
	
	$start = $start - 1;
	
	open(FOU,">$fileUp");
	my $i = 1;
	while ($i <= $lengthRegion)
	{
		print FOU "-\t$start\n";
		$start--;
		$i++;
	}
	close FOU;
	
	open(FOD,">$fileDown");
	my $j = 1;
	while ($j <= $lengthRegion)
	{
		print FOD "-\t$stop\n";
		$stop++;
		$j++;
	}
	close FOD;
	
	my $combFile = create_temp("f");
	system("cat $fileUp $natureFile $fileDown > $combFile");
	system("mv $combFile $natureFile");
	system("rm $fileUp $fileDown");
	
	system("mv $natureFile $out_dir/$acc.nature");
	system("mv $bedFile $out_dir/$acc.bed");
	
}

close FI;

sub GenerateNatureFiles		## Nature file is a file where each genomic coordinate is labelled. For coding regions, this format has three fields, the coding_exon Number, the genomic coordinate and the codon position (C1|C2|C3).
{
	my ($strand,$exonsListScalar,$lengthRegion) = @_;
	
	my @exonsList = split(",",$exonsListScalar);	
	
	my $noe = scalar(@exonsList);
	my @exons = ();
	@exons = @exonsList if ($strand eq "+");
	@exons = reverse(@exonsList) if ($strand eq "-");	## Reverse the order of the exons if the strand is minus.

	my @lineList = ();
	my $l = my $sum = 0;

	my $exonIndex = 1;
	foreach my $exon(@exons)
	{
		my @tmp = split(/-/,$exon);			## I do not check if the gene has one exons or more. I still get the acceptor/donor positions for genes which have only one exon.
								## However, they do not get printed at the end, because I do not print those positions which come from 1st exon acceptor or last exon donor.
								
		if ($strand eq "-")			## For genes on -ve strand.
		{					
			my $exonStart = $tmp[1];
			my $exonStop  = $tmp[0];
			
			my $CdsStart = $tmp[1] -1;		#### For CDS positions
			my $CdsStop	 = $tmp[0];

			while ($CdsStart >= $CdsStop)
			{
				$sum++;				## Compute this sum, this tells us about the codon positions.
				my $pos = getCodonPos($sum);
				my $line = "Coding_exon\_$exonIndex\t$CdsStart\t$pos";
				
				$lineList[$l] = $line;
				$l++;
				$CdsStart--;
			}
			
			my $start = my $stop = my $i = "";
			
			$start = $exonStop - 1;		### For donor positions
			$stop  = $start - ($lengthRegion - 1);
			$i = 1;
			
			while ($start >= $stop)
			{
				my $line = "$exonIndex\_D\_$i\t$start";		## Suppose it is the first exon, the line would read like "1_D_1 Genomic coordinate"
				$lineList[$l] = $line;
				$l++;
				
				$i++;
				$start--;
			}
			
			$start = $exonStart;		### For acceptor positions
			$stop  = $start + ($lengthRegion - 1);
			$i = 1;
			
			while ($start <= $stop)
			{
				my $line = "$exonIndex\_A\_$i\t$start";		## This should read like this, "2_A_1  Genomic coordinate"
				$lineList[$l] = $line;
				$l++;
				
				$i++;
				$start++;	
			}
		}
		
		else			## For genes on +ve strand.
		{
			my $exonStart = $tmp[0];
			my $exonStop  = $tmp[1];	
			
			my $CdsStart = $tmp[0];		#### For CDS positions
			my $CdsStop	 = $tmp[1] - 1;

			while ($CdsStart <= $CdsStop)
			{
				$sum++;				## Sum tells us about the codon position.
				my $pos = getCodonPos($sum);
				my $line = "Coding_exon\_$exonIndex\t$CdsStart\t$pos";
				
				$lineList[$l] = $line;
				$l++;
				$CdsStart++;
			}
			
			my $start = my $stop = my $i = "";
			$start = $exonStop;		### For donor positions
			$stop  = $start + ($lengthRegion - 1);
			$i = 1;
			
			while ($start <= $stop)
			{
				my $line = "$exonIndex\_D\_$i\t$start";		## For donor positions.
				$lineList[$l] = $line;
				$l++;
				
				$i++;
				$start++;
			}
			
			$start = $exonStart - $lengthRegion;		### For acceptor positions
			$stop  = $exonStart - 1;
			$i = $lengthRegion;
			
			while ($start <= $stop)
			{
				my $line = "$exonIndex\_A\_$i\t$start";
				$lineList[$l] = $line;
				$l++;
				
				$i--;
				$start++;	
			}
		}
		
		$exonIndex++;
	}

	my $natureFile = create_temp("f");
	open(FON,">$natureFile");
	
	foreach my $line(@lineList)
	{											## print everything that is in the @lineList array.
		print FON "$line\n" if ( ($line !~/^$noe\_D/) && ($line !~/^1_A/) );		## However, do not print the last exon donor positions or the firxt exon acceptor positions (they do not exist).	
	}
	close FON;

	return $natureFile;
}


sub GenerateBedFiles		## Given a gene in myGenePreduction format, this function writes a bedFile that should be submitted to mafExtract.
{
	my ($exonsListScalar,$lengthRegion,$chr) = @_;
	
	my @cdsExons = split(",",$exonsListScalar);
	my $noe = scalar(@cdsExons);
	
	my $bedFile = create_temp("f");
	
	open(FOB,">$bedFile");
	my $exonCt = 1;

	if ($noe > 1)
	{
		foreach my $exon(@cdsExons)		### Bed file related stuff --> Same for both the strands:
		{
			my $start = my $stop = "";	
			my @tmpExon  = split("-",$exon);
		
			if ($exonCt == 1)
			{
				$start = $tmpExon[0];
				$stop  = $tmpExon[1] + $lengthRegion;		## No acceptor region for the first exon.
			}
			elsif ($exonCt == $noe)
			{
				$start = $tmpExon[0] - $lengthRegion;		## No donor region for the last exon
				$stop  = $tmpExon[1];
			}
			else
			{							## Otherwise, donor and acceptor regions for intermediate exons.
				$start = $tmpExon[0] - $lengthRegion;
				$stop = $tmpExon[1] + $lengthRegion;
			}
		
			print FOB "$chr\t$start\t$stop\n"; 
		
			$exonCt++;
		}
	}

	else									## Of course, if the gene has only one exon, then get rid of donor /acceptor stuff.
	{
		foreach my $exon(@cdsExons)             ### Bed file related stuff --> Same for both the strands:
        {
            my @tmpExon  = split("-",$exon);
			print FOB "$chr\t$tmpExon[0]\t$tmpExon[1]\n";
		}
	}

	close FOB;	
	return $bedFile;
}


sub getCodonPos
{
	my $in = shift;
	$in = ($in %3);
	
	$in = 3 if ($in == 0);
	my $codonPos = "C".$in;
	return $codonPos;
}
