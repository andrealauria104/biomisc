#!/usr/bin/perl -w
use 5.010;
use strict;
use warnings;
use File::Spec::Functions 'catfile';

my $dirpath = $ARGV[0];
die "\n  [!] Error: Please specify a valid path for demultiplexing report. \n\n" if (!defined $dirpath || !-d $dirpath . "/Reports");
 
my $file_patterns = catfile($dirpath, 'Reports/*/*/default/*/all/lane.html');
my %results;

# Read bcl2fastq report ---
while ( my $file = glob($file_patterns) ) {
	
	my $name = $file;
	$name =~ s/(^(.*)default\/|\/all(.*)$)//g;
	print "[+] Processing sample: $name\n";

	open(my $fh, '<:encoding(UTF-8)', $file) or die "\n[!] Error: Could not open file '$file' ($!)\n\n";
	my @rows = <$fh>;
	close($fh);

	my @summary;
	foreach my $row (@rows) {
		chomp($row);
		if($row =~ m/(<td>\d+,(.*?)<\/td>?|<td>\d+<\/td>?)/g) {
			$row =~ s/(<td>|<\/td>)//g;
			$row =~ s/(,)/./g;
			push(@summary, $row);
		}
	}
	$results{$name} = \@summary;
}


# save summary table to file ---
my $outfile = catfile($dirpath, 'demultiplex_summary.txt');
print " -- output file: " . $outfile . "\n";

# header --
my @tbhead = ("sample", "TOT reads", "PF", "Yeld(Mb)", "lane1", "lane2", "lane3", "lane4");
open (my $OUT, "> $outfile") || die "[!] Error: Could not open file '$outfile' ($!)\n";
print $OUT join("\t",@tbhead) . "\n";
close($OUT);

# data --
foreach my $key (sort {uc($a) cmp uc($b)} keys %results)  
{ 	
	my $tbrow; 
	open (my $OUT, ">> $outfile") || die "[!] Error: Could not open file '$outfile' ($!)\n";
	$tbrow = join("\t"
		, $results{$key}[0]
		, $results{$key}[1]
		, $results{$key}[2]
		, $results{$key}[4]
		, $results{$key}[7]
		, $results{$key}[10]
		, $results{$key}[13]);
	print $OUT $key . "\t" . $tbrow . "\n";
	close($OUT);
}
