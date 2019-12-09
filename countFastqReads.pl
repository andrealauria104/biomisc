#!/usr/bin/perl -w
use 5.010;
use strict;
use warnings;
use File::Basename;

my $mypath = $ARGV[0];
my $count;
die "\n  [!] Error: Please specify a valid path to FASTQ. \n\n" if (!defined $mypath);

if( -f $mypath ) {
	
	$count = countreads($mypath);
	print basename($mypath) ."\t". $count . "\n";

} elsif( -d $mypath ) {

	my $dirfile;
	opendir(DH, $mypath);
	my @files = readdir(DH);
	closedir(DH);
	
	foreach my $file (@files) {
		next if($file !~ /(\.fastq|\.fq)$/);
		$dirfile = join("/", $mypath, $file);
		$count = countreads($dirfile);
		print $file ."\t". $count . "\n";

	}
}
# subroutine definition ---
sub countreads
{
	my $fastq = $_[0];
	my $line = 0;
	my $readcount = 0;

	if ($fastq =~ /.gz$/) {
		open(my $fh, "gzip -cd $fastq|") or die "\n[!] Error: Could not open file '$fastq' ($!)\n\n";
		my $id = <$fh>;
		my $seq = <$fh>;
		my $plus = <$fh>;
		my $qual = <$fh>;
		die "\n[!] Error: invalid FASTQ file\n\n" if ($id !~ m/^\@/ || $plus !~ m/^\+/);
		$line += 4;
		while(<$fh>) {
			$line++;
		}
		close($fh);
	} else {
		open(my $fh, '<:encoding(UTF-8)', $fastq) or die "\n[!] Error: Could not open file '$fastq' ($!)\n\n";
		my $id = <$fh>;
		my $seq = <$fh>;
		my $plus = <$fh>;
		my $qual = <$fh>;
		die "\n[!] Error: invalid FASTQ file\n\n" if ($id !~ m/^\@/ || $plus !~ m/^\+/);
		$line += 4;
		while(<$fh>) {
			$line++;
		}
		close($fh);
	} 

	$readcount = $line/4;
}