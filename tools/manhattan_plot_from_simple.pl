#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use FileHandle;

die "Usage: manhattan_plot_from_simple.pl <simple_file> <pdf_output> <gff3_file> <title>" unless defined $ARGV[3];

my ($input, $output, $gff3) = @ARGV;

my $gene_coords = get_coords($gff3);

my $qqman_input = "$input.qqman";
open OUT, ">$qqman_input";
print OUT "SNP\tCHR\tBP\tP\n"; # Header line

my %hash;
my ($geneid, $c, $coord, $p);
open FH, $ARGV[0];
while (<FH>) {
    chomp;
    my @i = split ' ', $_; # Split on any number of whitespace characters

    next if $i[0] eq "Gene" or $_ =~ /#/; # Skip the header
    $geneid = $i[0];
    if (defined $gene_coords->{$geneid}) {
        ($c, $coord) = @{$gene_coords->{$geneid}};
    } else {next;}

    $c =~ s/chr//;

    if ($c eq "X") { $c = 23; }
    elsif ($c eq "Y") {$c = 24; }
    $p = $i[4];
    print OUT "$geneid\t$c\t$coord\t$p\n";
}
close FH;
close OUT;

my $title = $ARGV[3];
my $cmd = "Rscript $FindBin::Bin/manhattan.R $qqman_input $output $title";
print "$cmd";
system($cmd);
#-------------------------
# Subs -------------------
#-------------------------
sub get_coors {
    my $file = shift;

    my $fh = FileHandle->new();

    $fh->open($file);
    my %gene_coords;
    while (my $line = <$fh>) {
        chomp $line;
        my @i = split /\t/, $line;

        next unless @i==9 and $i[2] eq "gene";

        $i[8] .= ";";
        my ($id) = $i[8] =~ /ID=(\S+?);/; # Non-space character followed by anything?

        $gene_coords{$id} = [$i[0], $i[3]];
    }
    $fh->close;

    return \%gene_coords;
}