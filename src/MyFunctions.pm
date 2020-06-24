#!/usr/bin/perl

package MyFunctions;
use strict;
use Exporter;
use Globals;
use bignum;

our @ISA = ('Exporter');
our @EXPORT = qw(convertToThousandCommaNum getN50 create_temp mylog10 arrayToHash getMotifFreqFile countMotifFreqArray getKaKs_BioPerl_NeiGojobori getIUPAC revComp getSum getMax getMin getAve getSD getSE getMedian getQuantile getSummary_ round min max printHash containsNonACGT swap getSubString makeSeqLogo getKaKs_YN00 getGCContent countMotif);

#########################################################################################################
# Implements a trick Michael found on a the web to add , after every 3 digits
# e.g. 9876543210.892  -->   9,876,543,210.892
#########################################################################################################
sub convertToThousandCommaNum {
	my $n = shift;
	while ($n =~ s/^(-?\d+)(\d{3})/$1,$2/) {};
	return $n;
}

#########################################################################################################
# computes the N50 and L50 for a given set of contig/scaffold lengths
# The second parameter (between 0 and 100) specifies the N** --> set to 90 for N90/L90
# example: my ($N90, $L90) = getN50(90, @contigLenghts);
#########################################################################################################
sub getN50 {
	my ($N, @vector) = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getSum\n";
		return -10000000000;
	}
	# actually to avoid that users give 0.5 instead of 50 for N50, we error if N is less than one
	die "ERROR in getN50: first parameter N was set to $N (must be between 1 and 99)\n" if ($N < 1 || $N > 99);
	my $sum = getSum(@vector);
	my $targetSum = $sum * ($N/100);		# N50 is the contig length when we have reached a sum that is that big
#  print "$N    $sum   $targetSum\n";

	# go through the sorted list and stop if $targetSum is reached
	my @sort = sort {$b <=> $a} @vector;
	my $N50=0;
	my $L50=0;
	my $curSum = 0;
	for (my $i=0; $i<=$#sort; $i++){
     $curSum +=$sort[$i];
     $L50++;
	  $N50 = $sort[$i];
      if($curSum >= $targetSum){
#         print "N$N length reached is $curSum and N$N value is: $N50 and L$N is $L50\n";
         last;
     }
	}
	return ($N50, $L50);
}

#########################################################################################################
# creates a temporary file or directory, if passed with "f" argument -> create a temp file, if passed with "d" argument -> create a temp directory.
#########################################################################################################
sub create_temp 	
{
	my $arg_1=my $temp_dir="";
	
	$arg_1=shift;
	
	if ($arg_1 eq "d")
	{
		$temp_dir=`mktemp -d`;
	}
	
	elsif ($arg_1 eq "f")
	{
		$temp_dir=`mktemp`;
	}
	
	else
	{
		die "Argument not defined for create_temp function\n";
	}
	
	$temp_dir=~s/\s+$//;
	return ($temp_dir);
}


#########################################################################################################
# given an array, the function inits a hash with each element (hash value set to 1)
# returns pointer to hash
#########################################################################################################
sub arrayToHash {
	my %H;
	foreach my $element (@_) {
		$H{$element} = 1;
	};
	return \%H;
}



#########################################################################################################
# returns log10 
#########################################################################################################
sub mylog10 {
  my $n = shift;
  return NaN if ($n == 0);
  return log($n)/log(10);
}


#########################################################################################################
# returns the IUPAC code for a nt or dinucleotide
# e.g. AC --> M
#########################################################################################################
sub getIUPAC {
	my ($input) = @_;
	my %H;
	$H{"A"} = "A";
	$H{"AC"} = "M";
	$H{"AG"} = "R";
	$H{"AT"} = "W";
	$H{"C"} = "C";
	$H{"CA"} = "M";
	$H{"CG"} = "S";
	$H{"CT"} = "Y";
	$H{"G"} = "G";
	$H{"GA"} = "R";
	$H{"GC"} = "S";
	$H{"GT"} = "K";
	$H{"T"} = "T";
	$H{"TA"} = "W";
	$H{"TC"} = "Y";
	$H{"TG"} = "K";
	
	if (! exists $H{$input}) {
		die "no IUPAC code for: $input\n";
	}
	return $H{$input};
}



#########################################################################################################
# return number of occurrences for the given motif
#########################################################################################################
sub countMotif {
	my ($seq, $motif) = @_;

	my $count = 0;
	my $pos = index($seq, $motif, 0);			
	while ($pos != -1) {
		$count++;
		$pos = index($seq, $motif, $pos+1);		
	}
	return $count;
}

#########################################################################################################
# return number of occurrences for the given motif and its frequency given a bunch of sequences as an array
# $seqs is pointer to an array
#########################################################################################################
sub countMotifFreqArray {
	my ($seqs, $motif) = @_;

	my $count = 0;
	my $numWords = 0;
	
	foreach my $seq (@{$seqs}) {
		my $pos = index($seq, $motif, 0);			
		while ($pos != -1) {
			$count++;
			$pos = index($seq, $motif, $pos+1);		
		}
		my $numWordsInSeq = length($seq) - length($motif) + 1;
#		print  length($seq), "  len gives $numWordsInSeq words: $seq\n";
		$numWords += $numWordsInSeq;
	}
	return ($count,$numWords,$count/$numWords);
}

#########################################################################################################
# input is a filename with oneline fasta
# return the freq of the given motif
# e.g. my ($count1, $num1, $freq1) = getMotifFreqFile($ARGV[0]);
#########################################################################################################
sub getMotifFreqFile {
	my ($file, $motif) = @_;

	my @seqs;
	my $line1 = "";
	open(file1, $file);
	while ($line1 = <file1>) {
		chomp($line1);
		next if ($line1 =~ /^>/);
		push @seqs, $line1;
	}
	return countMotifFreqArray(\@seqs, $motif);
}


 


#########################################################################################################
# returns GC content (0..100%) ignorning all non ACGT chars
#########################################################################################################
sub getGCContent {
	my ($seq) = @_;
	
	$seq = uc $seq;
	my $AT = ($seq =~ s/A/A/g);
	$AT += ($seq =~ s/T/T/g);
	my $GC = ($seq =~ s/G/G/g);
	$GC += ($seq =~ s/C/C/g);
#	print "AT $AT  $GC\n";
	return ($GC / ($AT+$GC)) * 100;
}

=pod
########################################################################################################
# creates a sequence logo
# parameters are a pointer to an array with the sequences, a filename and a title for the logo
#
# e.g. makeSeqLogo(\@conf, "PatternsConf$NAGxNAG", $NAGxNAG);
#########################################################################################################
sub makeSeqLogo {
	my ($data, $file, $title) = @_;

	# write the data to a file
	open file1, ">$file";
	foreach my $l (@{$data}) {
		print file1 ">a\n$l\n";
	}
	close file1;
	
	# get the path
	my @path = readpipe("pwd");
	chomp($path[0]);
	
	# call seqlogo
	system "(cd $SOFTWARE_DIR/weblogo/; seqlogo -f $path[0]/$file -o $path[0]/$file -F gif -k 1 -t \"$title\" -y \"Bits\" -a -c -n -Y)";
}
=cut
#########################################################################################################
# return reverse complement seq
#########################################################################################################
sub revComp {
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}

#########################################################################################################
# check if seq contains non ACGT
#########################################################################################################
sub containsNonACGT {
	my $seq = shift;

	$seq = uc $seq;
	$seq =~ tr/uU/TT/;
	$seq =~ s/A//g;
	$seq =~ s/C//g;
	$seq =~ s/G//g;
	$seq =~ s/T//g;
	if (length($seq) != 0) {
#		print STDERR "contains non ACGT: $seq\n";
		return 1;
	}
	return 0;
}



#########################################################################################################
# return the substring and checks for negative indices
#########################################################################################################
sub getSubString {
	my ($seq, $start, $length) = @_;
	
	if ($start >= 0) {
		return substr($seq,$start,$length);
	}else{
		my $newLength = ($length-$start >= 0 ? $length+$start : 0);
		return substr($seq,0,$newLength);
	}
}
	
#########################################################################################################
# compute max of vector
# returns max and the element number having the max
#
#	my ($x_max,$dummy1) = getMax(@tmp);
#########################################################################################################
sub getMax {
	my @vector = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getMax\n";
		return -10000000000;
	}

	my $max = $vector[0];
	my $element = 0;
	for (my $i=0; $i<=$#vector; $i++) {
		$element = $i if ($vector[$i] > $max);
		$max = $vector[$i] if ($vector[$i] > $max);
	}
	return ($max,$element);
}


#########################################################################################################
# compute min of vector
# returns min and the element number having the min
#
#	my ($x_min,$dummy) = getMin(@tmp);
#########################################################################################################
sub getMin {
	my @vector = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getMin\n";
		return 10000000000;
	}

	my $min = $vector[0];
	my $element = 0;
	for (my $i=0; $i<=$#vector; $i++) {
		$element = $i if ($vector[$i] < $min);
		$min = $vector[$i] if ($vector[$i] < $min);
	}
	return ($min,$element);
}


#########################################################################################################
# compute sum of vector
#########################################################################################################
sub getSum {
	my @vector = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getSum\n";
		return -10000000000;
	}

	my $sum = 0;
	my $count = 0;	
	for (my $i=0; $i<=$#vector; $i++) {
		$sum += $vector[$i];
		$count++;
	}
	return $sum;
}


#########################################################################################################
# compute average of vector
#########################################################################################################
sub getAve {
	my @vector = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getAve\n";
		return -10000000000;
	}

	my $sum = 0;
	my $count = 0;	
	for (my $i=0; $i<=$#vector; $i++) {
		$sum += $vector[$i];
		$count++;
	}
	$sum /= $count;
	return $sum;
}

#########################################################################################################
# compute standard deviation of vector
#########################################################################################################
sub getSD {
	my @vector = @_;
	if ($#vector < 1) {			# a vector with one value has SD = 0
		print STDERR "ERROR: empty vector for getSD\n";
		return -1;
	}
	
	# average
	my $ave = 0;
	my $count = 0;	
	for (my $i=0; $i<=$#vector; $i++) {
		$ave += $vector[$i];
		$count++;
	}
	$ave /= $count++;
	
	# variance
	my $sum = 0;
	$count = 0;	
	for (my $i=0; $i<=$#vector; $i++) {
		$sum += (($vector[$i] - $ave) * ($vector[$i] - $ave));
		$count++;
	}
	$sum /= ($count - 1);
	return sqrt($sum);			# return sd
}

#########################################################################################################
# compute standard error of vector
#########################################################################################################
sub getSE {
	my @vector = @_;
	if ($#vector < 1) {			# a vector with one value has SD = 0
		print STDERR "ERROR: empty vector for getSE\n";
		return -1;
	}
	
	my $sd = getSD(@vector);
	my $se = $sd / sqrt(($#vector+1));
	return $se;		
}


#########################################################################################################
# compute median of vector
#########################################################################################################
sub getMedian {
	my @vector = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getMedian\n";
		return -1000000000;
	}
	my @sorted = sort{$a<=>$b}(@vector);
	my $num = $#sorted + 1;
	if ($num % 2 != 0) {		# odd
		return $sorted[int($num/2)];
	}else{
#	print "$num  ", int($num/2)-1,  "    ", int($num/2), "  ", $sorted[int($num/2)], "    ",$sorted[int($num/2)-1] ,"\n";
#	print "@sorted\n";
		return ($sorted[int($num/2)-1] + $sorted[int($num/2)]) / 2;
	}
}


#########################################################################################################
# compute quantil of vector
# just return the rounded value $p*$num
# NOTE: median and 0.5 quantile give different values if there are a even number of elements
# 
# given quantile as first parameter (0 - 1) and the vector 
#########################################################################################################
sub getQuantile {
	my ($p, @vector) = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getQuantile\n";
		return -1000000000;
	}
	my @sorted = sort{$a<=>$b}(@vector);
	my $num = $#sorted+1;

	return $sorted[round($p*$num)-1];
}




#########################################################################################################
# rounding
#########################################################################################################
sub round {
    my ($number) = shift;
	 my $toAdd = ($number < 0 ? -0.5 : 0.5);		# in case of -2.2 adding +0.5 results in -1.7 == -1 but adding -0.5 results in -2.7 == -2
    return int($number + $toAdd);
}

#########################################################################################################
# min
#########################################################################################################
sub min {
    my ($a,$b) = @_;
	 return ($a < $b ? $a : $b);
}

#########################################################################################################
# max
#########################################################################################################
sub max {
    my ($a,$b) = @_;
	 return ($a > $b ? $a : $b);
}

#########################################################################################################
# print hash
# optional parameter delimiter (default \n)
#########################################################################################################
sub printHash {
	my @a = @_;
	my $hash = shift;
	my $delimiter = "\n";
	if ($#a > 0) {
		$delimiter = $a[1];
	}

	foreach my $key (keys %{$hash}) {
		print "$key  =>  $hash->{$key};$delimiter";
	}

	if ($delimiter ne "\n") {print "\n";};
}

#########################################################################################################
# simple swap  ($x, $y) = swap($x,$y);
#########################################################################################################
sub swap {
	my ($a, $b) = @_;
	return ($b, $a);
}


#########################################################################################################
# return min, first-second-third quartile, max, average, sd
#########################################################################################################
sub getSummary_ {
	my (@vector) = @_;
	if ($#vector < 0) {
		print STDERR "ERROR: empty vector for getSummary\n";
		return -1000000000;
	}
	my @sorted = sort{$a<=>$b}(@vector);

	my $ave = getAve(@sorted);
	my $sd = getSD(@sorted);
	my $first= getQuantile(0.25,@sorted);
	my $median = getMedian(@sorted);
	my $third = getQuantile(0.75,@sorted);
	my ($min,$el1) = getMin(@sorted);
	my ($max,$el2) = getMax(@sorted);

	return ($min, $first, $median, $third, $max, $ave, $sd);
}



=pod
############################################################################
# get the Ka Ks values using my own version of YN00 from PAML
# expect codon alignment sequences like CAA GCA GAG GCT CCC GTG CAG GAA GAG AAG CTG
# return pointer to results hash
# my $r = getKaKs_YN00("GTTCTTGACAACCTGGACAGCCATCTGAAAAAATCCGAATACTTTCGCTTCCTCTGGTTCCCGCATAGTGAGAATGTCAGCATCATCTACCAGGACCACACCAACAAG","GTTCTCGACAACCTGGACACCCATCTGAAGAAGTCGGAATACTTTCGCTTCCTCTGGTTTCCGCACAGTGAGAACGTCAGTGTCATCTACCAGGACCACACCAACAAG");
# print "$r->{omega} $r->{dN} $r->{dS}\n";
############################################################################
sub getKaKs_YN00 {
	my ($seq1, $seq2) = @_;
	
	my @dir = readpipe("pwd");
	my $WorkingPATH = $dir[0]; chomp($WorkingPATH);
	
	# get the length of nucleotide sequence
	my $s = $seq1;
	$s =~ s/ //g;			  # remove all spaces
	my $len = length($s);
	
#	print "$seq1\n$seq2\n$s\n$len\n\n";
	
	# print input file in phylip format
	open YNInput, ">$WorkingPATH/tmp.yn00.input";
	print YNInput "2 $len\nsrc       $seq1\ndest      $seq2\n";
	close YNInput;
	
	# use my own command-line-parameter version of YN00
	my @res = readpipe "$PAMLPath/bin/myyn00 -seqfile $WorkingPATH/tmp.yn00.input";

	# parse the output
	my $results;
	my $dataline;
	for (my $i=0; $i<=$#res; $i++) {
		my $line = $res[$i];
		
		if ($line =~ /^seq\.\s+seq\./) {
			$dataline = $res[$i+2];
			if ($dataline !~ /^\s+\d+\s+\d+\s+[\d\.]+\s+[\d\.]+\s+[\d\.\-]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.\-]+\s*\+\-\s*[\d\.]+\s+[\d\.\-]+\s*\+\-\s*[\d\.]+/)  	{
				print STDERR "ERROR: cannot parse YN dataline: $dataline\n";
		      $results = { 'S' => 0, 'N' => 0,'t' => 0, 'kappa' => 0,  'omega' => 0,  'KaKs' => 0,  'dN' => 0, 'dN_SE' => 0, 'dS' => 0, 'dS_SE' => 0, 'S_changed' => 0, 'N_changed' => 0, 'valid' => 0};
				return $results;
			}
			
			#  seq1 seq1  S        N       t     kappa   omega    dN       dN_SE    dS       dS_SE
			#   5    2   190.4   664.6   1.2984  1.5804  0.1416 0.1842 +- 0.0184  1.3008 +- 0.1995
			my @f = split(/\s+/, $dataline);
#			print "$dataline\n";
#			for (my $i=0; $i<=$#f; $i++) {
#				print "$i:     $f[$i]\n";
#			}

			# output from myyn00.c is 			
         #   fprintf(fout, " %3d  %3d ", is+1,js+1);		 the two sequences
         #   fprintf(fout,"%7.1f %7.1f %8.4f %7.4f %7.4f %6.4f +- %6.4f %7.4f +- %6.4f\n", S,N,t,com.kappa,com.omega,dN,SEdN,dS,SEdS);
			if ($#f == 13) {
	   	   $results = { 
				 	'S' => $f[3], 
				 	'N' => $f[4],	
					't' => $f[5],			# time, measured by the expected number of nt substitutions per codon
					'kappa' => $f[6],		# transition / transversion mutation rate ratio
					'omega' => $f[7],		# dN / dS
					'KaKs' => $f[7],		#	== omega
					'dN' => $f[8],			# rate of non-synonymous mutations
					'dN_SE' => $f[10],	# std error
					'dS' => $f[11],		# rate of synonymous mutations	
					'dS_SE' => $f[13]};	# std error
				$results->{'S_changed'} = $results->{'S'} * $results->{'dS'};		# number of synonymous mutations
				$results->{'N_changed'} = $results->{'N'} * $results->{'dN'};		# number of non-synonymous mutations
		      $results->{'valid'} = 1;
		      $results->{'wholeLine'} = $dataline;
					
			}else{
				die "ERROR: cannot parse YN dataline: $dataline\n";
			}
			last;
		}
	}

	# remove tmp files
#	system "rm $WorkingPATH/tmp.yn00.input $WorkingPATH/2YN.dN $WorkingPATH/2YN.dS $WorkingPATH/2YN.t $WorkingPATH/rst $WorkingPATH/rst1 $WorkingPATH/rub";
	
	# set valid = -2 if dS = 0 (i.e. there are no synonymous changes and Ka/Ks = 99)
	if ($results->{dS} == 0) {
		$results->{valid} = -2;
	}
	
	# there can be ERRORs by yn00 (for example when the ClustalW ali is wrong) that produce dN or dS > 1 --> exclude those cases
	if ($results->{dS} > 50 || $results->{dS} < 0 || $results->{dN} > 50 || $results->{dN} < 0) {
		print  STDERR "Warning: yn00 error: $dataline";
      $results->{valid} = -1;
	}
	
	return $results;
}
=cut


############################################################################
# code fragment that shows how to compute Ka/Ks using bioperl with a pairwise all-to-all method
############################################################################
sub getKaKs_BioPerl_NeiGojobori {
	use lib "$ENV{'genomePath'}/src/BioPerl";
	use Bio::AlignIO;
	use Bio::Align::DNAStatistics;

	my $stats = new Bio::Align::DNAStatistics;
	my $in = new Bio::AlignIO(-format => 'fasta',
	#                            -file   => '../Bioperl/bioperl-1.4/t/data/nei_gojobori_test.aln');
   	                         -file   => 'x');
	my $alnobj = $in->next_aln;
	my $results2 = $stats->calc_all_KaKs_pairs($alnobj);
	for my $an (@$results2){
		print "comparing ". $an->{'Seq1'}." and ". $an->{'Seq2'}. " D_n $an->{D_n} D_s $an->{D_s} D_n/D_s ", $an->{D_n}/$an->{D_s},"\n";
		for (sort keys %$an ){
			next if /Seq/;
			printf("%-9s %.4f \n",$_ , $an->{$_});
		}
		print "\n\n";
	}
}
















