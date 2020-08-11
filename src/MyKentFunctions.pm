#!/usr/bin/perl

package MyKentFunctions;
use strict;
use Exporter;
use Globals;
use MyFunctions;
use Scalar::Util::Numeric qw(isint);
	
our @ISA = ('Exporter');
our @EXPORT = qw(readChromSizeHash getSeqTwoBitFile getSeqTwoBit getSeqTwoBitQuality getSeqTwoBitQualityBatch getSLineData getELineData getILineData);


#############################################################
# read the chrom.sizes (twoBitInfo output) into a hash
# call: my %chromSize = readChromSizeHash($chromSizeFile);
#############################################################
sub readChromSizeHash {
	my $chromSizeFile = shift; 

	my %hash;
	open(FI,"$chromSizeFile") || die "ERROR: cannot open file $chromSizeFile\n";
	while (my $line=<FI>) {
		chomp($line);
		my @f = split(/\t/, $line);
		$hash{$f[0]} = $f[1];
	}
	close(FI);
	return %hash;
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

#==============  NOTE: This is function we keep as relic to assure compatibility with old code. Use getSeqTwoBitFile instead !! ============== 
#########################################################################################################
# get a sequence given chr:start-end, species and the path to gbdb
#
# + strand: 
#	$seq = getSeqTwoBit($species, $chr, $start-$countUp, $start, $strand) if ($strand eq "+");
# - strand: 	
#	$seq = getSeqTwoBit($species, $chr, $srcSize-$start, $srcSize-$start+$countUp, $strand) if ($strand eq "-");		 
#########################################################################################################
sub getSeqTwoBit {
	my ($species, $chr, $start, $end, $strand, $GBDBpath) = @_;

	# determine the 2bit file; depends on species
	my $TwoBitFile = "$GBDBpath/$species/$species.2bit";
	if (!-e $TwoBitFile) {
		$TwoBitFile = "$GBDBpath/$species/nib/$species.2bit";
		if (!-e $TwoBitFile) {
			print STDERR "ERROR: no 2bit file for $species: $TwoBitFile\n";
			exit -1;
			return "";
		}
	}	
	my $call = "twoBitToFa $TwoBitFile stdout -seq=$chr -start=$start -end=$end";
	$call .= " | faRc stdin stdout -keepCase" if ($strand eq "-");
	print "$call\n" if ($verbose);
	my @res = `$call`;
	die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	print STDERR "ERROR: cannot get seq by calling: $call\n" if ($#res < 1);
	
	# skip first fasta header and remove the \n between the lines
	shift @res; $" = ""; chomp @res; my $seq = "@res";
	return $seq;
}


#########################################################################################################
# get a sequence given chr:start-end, the 2bit file, and a 0/1 flag for unmasking, it runs twoBitToFa and extracts the sequence
#
# unmask is a parameter that does
#  - if set to 1, output all characters as upper case (that means quality-masked or repeats (both lower case) will be upper case)
#  - if set to 2, output lower case characters as N's. Use this to hard-mask quality-masked (lower case character in $species.quality.2bit) bases
#  - if set to neither 1 or 2, return the sequence as it comes from the 2bit file
#
#
# Example: with unmasking on (1)
# + strand: 
#	$seq = getSeqTwoBitFile($chr, $start-$countUp, $start, $strand, "$genomePath/gdbd-HL/$species/$species$TwoBitSuffix", 1) if ($strand eq "+");
# - strand: 	
#	$seq = getSeqTwoBitFile($chr, $srcSize-$start, $srcSize-$start+$countUp, $strand, "$genomePath/gdbd-HL/$species/$species$TwoBitSuffix", 1) if ($strand eq "-");		 
#########################################################################################################
sub getSeqTwoBitFile {
	my ($chr, $start, $end, $strand, $TwoBitFile, $unmask) = @_;

	# run twoBitToFa with -noMask if unmask set to 1
	my $unmaskParameter = "";
	$unmaskParameter = "-noMask" if ($unmask == 1);

	if (! -e $TwoBitFile) {
		die "ERROR in getSeqTwoBitFile: $TwoBitFile does not exist\n";
	}
	
	my $call = "set -o pipefail; echo -e \"$chr\\t$start\\t$end\\tname\\t1000\\t$strand\" | twoBitToFa $TwoBitFile stdout -bed=stdin $unmaskParameter";
	print "$call\n" if ($verbose);
	my @res = `$call`;
	die "ERROR in getSeqTwoBitFile: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	print STDERR "ERROR: cannot get seq by calling: $call\n" if ($#res < 1);
	
	# skip first fasta header and remove the \n between the lines
	shift @res; $" = ""; chomp @res; my $seq = "@res";

	# if unmask set to 2, replace lower case letters by N
	if ($unmask == 2) {
		$seq =~ s/[acgtn]/N/g;
	}

	return $seq;
}

	
#########################################################################################################
# just like getSeqTwoBit but it uses the .quality.2bit files and replaces all lower case letters by N
# .quality.2bit files are produced from the genome 2bit files by masking (lower case) all letters below a certain quality level
#########################################################################################################
sub getSeqTwoBitQuality {
	my ($species, $chr, $start, $end, $strand, $GBDBpath) = @_;

	# determine the quality.2bit file; depends on species
	my $TwoBitFile = "$GBDBpath/$species/$species.quality.2bit";
	if (!-e $TwoBitFile) {
		$TwoBitFile = "$GBDBpath/$species/nib/$species.quality.2bit";
		if (!-e $TwoBitFile) {
			print STDERR "ERROR: no quality.2bit file for $species: $TwoBitFile\n";
			exit -1;
			return "";
		}
	}	
	my $call = "twoBitToFa $TwoBitFile stdout -seq=$chr -start=$start -end=$end";
	$call .= " | faRc stdin stdout -keepCase" if ($strand eq "-");
	print "$call\n" if ($verbose);
	my @res = `$call`;
	die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	print STDERR "ERROR: cannot get seq by calling: $call\n" if ($#res < 1);
	
	# skip first fasta header and remove the \n between the lines
	shift @res; $" = ""; chomp @res; my $seq = "@res";
	# hard mask
	$seq =~ tr/acgt/NNNN/;		

	return $seq;
}
	
	
	
#########################################################################################################
# just like getSeqTwoBitQuality but you give it a pointer to a list of loci (chr:start-end) and strands for one species
# batch calls are much faster with twoBitToFa
# returns a list of sequences where lower case letters (low quality) are replaced by N's and where the seq is reverse complemented if the strand is -
# every entry in the strand array must directly match the coords array
#########################################################################################################
sub getSeqTwoBitQualityBatch {
	my ($species, $coords, $strands, $GBDBpath) = @_;

	# determine the quality.2bit file; depends on species
	my $TwoBitFile = "$GBDBpath/$species/$species.quality.2bit";
	if (!-e $TwoBitFile) {
		$TwoBitFile = "$GBDBpath/$species/nib/$species.quality.2bit";
		if (!-e $TwoBitFile) {
			print STDERR "ERROR: no quality.2bit file for $species: $TwoBitFile\n";
			exit -1;
			return "";
		}
	}	
	my $call = "echo -e \"$coords\" | twoBitToFa $TwoBitFile stdout -seqList=/dev/stdin";
	print "$call\n" if ($verbose);
	my @res = readpipe($call);
	print "res @res\n" if ($verbose);
	print STDERR "ERROR: cannot get seq by calling: $call\n" if ($#res < 1);
	
	
	# skip fasta headers and remove the \n between the lines
	my @seqs;
	my $seq = "";
	foreach my $l (@res) {
		if ($l =~ /^>/) {		
			if ($seq ne "") {
				# hard mask
				$seq =~ tr/acgtn/NNNNN/;		
				my $strand = shift @{$strands};
				$seq = revComp($seq) if ($strand eq "-");
				print "SSS $strand  push $seq\n" if ($verbose);
				push @seqs, $seq if ($seq ne "");
			}
			$seq = "";
		}else{
			chomp $l;
			$seq .= $l;
		}
	}
	# hard mask the last sequence (there is no fasta header following)
	my $strand = shift @{$strands};
	$seq = revComp($seq) if ($strand eq "-");
	$seq =~ tr/acgtn/NNNNN/;		
	push @seqs, $seq if ($seq ne "");

	print "result: @seqs\n" if ($verbose);

	return @seqs;
}
	
	

