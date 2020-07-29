#!/usr/bin/env perl
=head1 NAME

features-from-gtf.pl - GTF feature extraction tool

=cut
=head1 DESCRIPTION

Retrieves genomic coordinates for TSS, exons 
(collapsing overlapping isoforms) and introns 
from GTF annotation file. 
=cut

=head1 SYNOPSIS

features-from-gtf.pl [Options] <path_to_gtf>

=head1 VERSION

0.01
=cut

=head1 OPTIONS

=over 5

=item B<-h|--help>

Show help information and exit.

=item B<-o|--outdir>

Output directory (default: .)

=item B<-f|--force>

Overwrites existing files.

=item B<-m|--manual>

Read manual.

=item B<-v|--version>

Show the version number and exit

=back

=cut
=head1 AUTHOR

Andrea Lauria
=cut

# 0. Resources ----
use strict;
use warnings;
use File::Basename;
use List::Util qw[min max];
use Array::Utils qw(:all);
use Getopt::Long;
use Pod::Usage;

my $version = '0.1';
my $out_fhdir = "./";
my ($force, $gtf, $in_fh, $out_fh, $checked_bed);

parseArguments();

die "\n  [!] Error: Please provide a valid path for input GTF file. \n\n" if (!defined $gtf || !-f $gtf);

print "============================================== \n";
print "==== Extract features from GTF annotation ==== \n";
print "============================================== \n";
print "\n -- annotation = " . $gtf . "\n";

# 1. Extract TSS ---
my $tss_bed = $out_fhdir . basename($gtf);
$tss_bed =~ s/\.gtf(\.gz)?$/\.tss.bed/g;

print "\n[+] Extract TSS, output file (BED format) = " . $tss_bed . "\n";
$checked_bed = check_bed($tss_bed, $force);

if($checked_bed) {
    open($in_fh, "awk \'\$3 == \"gene\"\' " . $gtf . " |"); # open "gene" features
    while(<$in_fh>) {
        my $tss;
        my @line = split /\t/;
        next if ($line[2] ne "gene");

        $_ =~ m/gene_id "(.+?)";/;    
        my $gene_id = $1;
        my $chr = $line[0];
        my $strand = $line[6];

        # Retrieve TSS
        if ($strand eq "+") {
            $tss = $line[3];
        } elsif ($strand eq "-") {
            $tss = $line[4];
        }
        # Write in BED format
        open (my $out_fh, ">> $tss_bed") || die "[!] Error: Could not open file '$tss_bed' ($!)\n";
        print $out_fh join("\t", $chr, $tss-1, $tss, $gene_id, ".", $strand) . "\n";
        close($out_fh);
        
    }
    close($in_fh);
}

# # 2. Extract exons ----
my $exon_bed = $out_fhdir . basename($gtf);
$exon_bed =~ s/\.gtf(\.gz)?$/\.exons.bed/g;

print "\n[+] Extract exons, output file (BED format) = " . $exon_bed . "\n";
$checked_bed = check_bed($exon_bed, $force);

if($checked_bed) {
    my %gene_exons;
    open($in_fh, "awk \'\$3 == \"exon\"\' " . $gtf . " |"); # open "exon" features
    while(<$in_fh>) {
        chomp;
        my @line = split /\t/;
        next if ($line[2] ne "exon");

        # Retrieve exons
        $_ =~ m/gene_id "(.+?)";/;    
        my $gene_id = $1;
        $_ =~ m/exon_id "(.+?)";/;    
        my $exon_id = $1;
        my $exon_name = join("_","gene_id", $gene_id, "exon_id", $exon_id); 

        my $chr = $line[0];
        my $start = $line[3];
        my $end = $line[4];
        my $strand = $line[6];

        push(@{$gene_exons{$gene_id}}, [$chr, $start-1, $end, $strand]);
    }
    close($in_fh);
    foreach my $id (keys %gene_exons) {
        my @exons = sort {$a->[1] <=> $b->[1]} @{$gene_exons{$id}};
        # Merge overlapping exons
        for (my $i = 0; $i < $#exons; $i++) {

            my @interval_1 = @{$exons[$i]}[1,2];
            my @interval_2 = @{$exons[$i+1]}[1,2];

            if (intersect(@interval_1, @interval_2)) {
                my $start = min($interval_1[0], $interval_2[0]);
                my $end = max($interval_1[1], $interval_2[1]);

                @{$exons[$i+1]}[1] = $start;
                @{$exons[$i+1]}[2] = $end;

                splice(@exons, $i, 1);
                $i--;
            }
        }
        my $exon_n = 1;
        for (@exons) {
            my $exon_name = $id . "_exon_" . $exon_n;
            # Write in BED format
            open ($out_fh, ">> $exon_bed") || die "[!] Error: Could not open file '$exon_bed' ($!)\n";
            print $out_fh join("\t", @{$_}[0..2], $exon_name, ".", @{$_}[3]) . "\n";
            close($out_fh);
            $exon_n++;
        }
    }

    system("sort -k1,1 -k2,2n -V $exon_bed > tmp_exon.bed");
    system("mv tmp_exon.bed $exon_bed");
}

# 3. Extract introns ----
my $intron_bed = $out_fhdir . basename($gtf);
$intron_bed =~ s/\.gtf(\.gz)?$/\.introns.bed/g;

print "\n[+] Extract introns, output file (BED format) = " . $intron_bed . "\n";
$checked_bed = check_bed($intron_bed, $force);

if($checked_bed) {
    system("mergeBed -c 4 -o distinct -i $exon_bed > tmp_exon.bed");
    my %gene_exons_i;
    open($in_fh, "< tmp_exon.bed"); # open "exon" features
    while(<$in_fh>) {
        chomp;
        my @line = split /\t/;
        my $gene_id = $line[3];
        $gene_id =~ s/_exon.*//;    
        my $chr = $line[0];
        my $start = $line[1];
        my $end = $line[2];

        push(@{$gene_exons_i{$gene_id}}, [$chr, $start, $end]);
        
    }
    close($in_fh);

    foreach my $id (keys %gene_exons_i) {
        my @exons_i = sort {$a->[1] <=> $b->[1]} @{$gene_exons_i{$id}};
        
        for (my $i = 0; $i < $#exons_i; $i++) {
            my $intron_start = $exons_i[$i][2];
            my $intron_end = $exons_i[$i+1][1]; 
            my $intron_name = $id . "_intron_" . $i;

            open ($out_fh, ">> $intron_bed") || die "[!] Error: Could not open file '$exon_bed' ($!)\n";
            print $out_fh join("\t", $exons_i[$i][0], $intron_start, $intron_end, $intron_name,".","+") . "\n";
            close($out_fh);
            }
        }

    system("sort -k1,1 -k2,2n -V $intron_bed > tmp_intron.bed");
    system("bedtools subtract -a tmp_intron.bed -b $exon_bed > $intron_bed");
    system("rm tmp_exon.bed tmp_intron.bed");
}

# Subroutines ----
sub check_bed 
{
	my $bed = $_[0];
    my $forced = $_[1];
	if (-f $bed && $forced) {
        print "\n  [!] Overwriting existing file\n";
        system("rm $bed");
        return 1;
	} elsif (!-f $bed) {
        return 1;
    } else {
        print "\n  [!] Output file exists. Re-run with -f|--force option for replacing it.\n";
        return 0;
    }
}

sub parseArguments
{
    my $HELP    = 0;   # Show help overview.
    my $MANUAL  = 0;   # Show manual
    my $VERSION = 0;   # Show version number and exit.
    
    GetOptions(
        "h|help"   => \$HELP,
        "m|manual"   => \$MANUAL,
        "v|version"  => \$VERSION,
        "o|outdir=s" => \$out_fhdir,
        "f|force"  => \$force
        );
    
    $gtf = $ARGV[0];
    pod2usage(1) if $HELP;
    pod2usage(-verbose => 2 ) if $MANUAL;

    if ( $VERSION )
    {
        print "Version: $version \n";
        exit;
    }
}
