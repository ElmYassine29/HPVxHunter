#!/usr/bin/perl
use strict;
use warnings;

die "Usage: $0 <fasta_file> <id_list_file>\n" unless @ARGV == 2;
my ($fasta, $id_list) = @ARGV;

open my $id_fh, '<', $id_list or die "Can't open $id_list: $!";
my %wanted;
while (<$id_fh>) {
    chomp;
    $wanted{$_} = 1 if $_ ne '';
}
close $id_fh;

open my $in,  '<', $fasta or die "Can't open $fasta: $!";

my ($header, $seq) = ('', '');
while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        # Nouveau header → traiter la précédente séquence si valide
        if ($header && exists $wanted{$header}) {
            my $length = length($seq);
            if ($length >= 200 && $seq !~ /[^ACGTacgt]/) {
                print ">$header|$header\n$seq\n";
            }
        }
        # Réinitialiser pour la prochaine
        $header = $1;
        $seq = '';
    } else {
        $seq .= $line;
    }
}
# Dernière séquence à traiter
if ($header && exists $wanted{$header}) {
    my $length = length($seq);
    if ($length >= 200 && $seq !~ /[^ACGTacgt]/) {
        print ">$header|$header\n$seq\n";
    }
}
close $in;
