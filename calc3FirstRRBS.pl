#!/usr/bin/perl -w
use 5.010;
use strict;
use warnings;
use File::Basename;

my $mypath = $ARGV[0];
my @res;
die "\n  [!] Error: Please specify a valid path to FASTQ. \n\n" if (!defined $mypath);

if( -f $mypath ) {
	
	@res = readNcalc3First($mypath);
	print "sample\ttotreads\trrbsreads\tratio\n";
	print basename($mypath) ."\t". $res[0] ."\t". $res[1] ."\t". $res[2] . "\n";

} elsif( -d $mypath ) {

	my $dirfile;
	opendir(DH, $mypath);
	my @files = readdir(DH);
	closedir(DH);

	print "sample\ttotreads\trrbsreads\tratio\n";
	foreach my $file (@files) {
		next if($file !~ /\.(fastq|fq)(\.gz|)$/);
		$dirfile = join("/", $mypath, $file);
		@res = readNcalc3First($dirfile);
		print $file ."\t". $res[0] ."\t". $res[1] ."\t". $res[2] . "\n";

	}
}

# subroutine definition ---
sub readNcalc3First
{	
	my $fastq = $_[0];

	if ($fastq =~ /.gz$/) {
		open(my $fh, "gzip -cd $fastq|") or die "\n[!] Error: Could not open file '$fastq' ($!)\n\n";
		calc3First($fh);
	} else {
		open(my $fh, '<:encoding(UTF-8)', $fastq) or die "\n[!] Error: Could not open file '$fastq' ($!)\n\n";
		calc3First($fh);
	}
}

sub calc3First
{
	my $fh = $_[0];
	my $readcount = 0;
	my $rrbscount = 0;

	while(!eof($fh)) {
		my $id = <$fh>;
		my $seq = <$fh>;
		my $plus = <$fh>;
		my $qual = <$fh>;
		chomp($seq);
		die "\n[!] Error: invalid FASTQ file\n\n" if ($id !~ m/^\@/ || $plus !~ m/^\+/);
		$readcount++;
		$rrbscount++ if ($seq =~ m/^[CT]GG/);
	}

	my $ratio = $rrbscount/$readcount;
	
	return($readcount, $rrbscount, $ratio);
}



