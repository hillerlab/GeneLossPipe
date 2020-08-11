#!/sw/bin/perl

# Michael Hiller, 2009
# combine several BDB files that start with a common prefix and write just 1 BDB file

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BerkeleyDB;

$| = 1;		# == fflush(stdout)
# options
my $verbose = 0;
sub usage {
	die "usage: $0 DirAndORBDBprefix outFile\nDirAndORBDBprefix can be dir/prefix or just the prefix (e.g. BDB/STITCH or ./STITCH or STITCH)\n";
};
GetOptions ("v|verbose"  => \$verbose) || usage();
usage() if ($#ARGV < 1);
my $BDBprefix = $ARGV[0];				# prefix for the BDB file to be read
my $outFile = $ARGV[1];					# output file
usage() if ($BDBprefix eq "" || $outFile eq "");

die "ERROR: $outFile exists. Delete it first.\n" if (-e $outFile);

# just count how many keys were written
my $numKeys = 0;

# do it
# first lock the output file
print "waiting for lockfile lockFile.$outFile ....";
system "set -o pipefail; lockfile -1 lockFile.$outFile" || die "ERROR: cannot get lock file: lockFile.$outFile";
print " got it\n";
# create the output BDB
my %outhash;
tie %outhash, "BerkeleyDB::Btree",	-Filename => "$outFile", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";

# read all input BDB files
my $basename = `set -o pipefail; basename $BDBprefix`;
die "ERROR: basename $BDBprefix\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp($basename);
my $dirname = `set -o pipefail; dirname $BDBprefix`;
die "ERROR: dirname $BDBprefix\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp($dirname);
print "reading from directory $dirname all files starting with $basename ...\n";

opendir(DIR, "$dirname") or die "can't open dir $dirname: $!\n";
while (defined(my $file = readdir(DIR))) {
     if ($file =~ /^$basename/) {
         print "read $dirname/$file\n";
			doit("$dirname/$file");
     }
}
closedir(DIR);

# close output file and release the lock
untie %outhash;
system "set -o pipefail; rm -f lockFile.$outFile" || die "ERROR: cannot delete lockFile.$outFile";

print "All done. Wrote a total of $numKeys key-value pairs to $outFile\n";



########################################################################
# read the entire file and write every key-value pair to the output file
sub doit {
	my ($file) = @_;

	# input	
	my %h;
	tie %h, "BerkeleyDB::Btree",	-Filename => "$file", -Flags => DB_RDONLY || die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";
	foreach my $key (keys %h) { 
		$numKeys ++;
		print "$key --> $h{$key}\n" if ($verbose);
		$outhash{$key} = $h{$key};
	}
	untie %h;
}	
	
