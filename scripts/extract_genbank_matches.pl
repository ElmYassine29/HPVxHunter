#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my ($tsv_file, $fasta_file) = @ARGV;
die "Usage: $0 <input.tsv> <genbank.fasta>\n" unless @ARGV == 2;

# Lire les IDs à partir du fichier TSV (colonne 2)
my %wanted_ids;
open(my $tsv, '<', $tsv_file);
while (<$tsv>) {
    chomp;
    my @cols = split(/\t/);
    next if $. == 1;  # Skip header
    $wanted_ids{$cols[1]} = 1 if $cols[1];  # On stocke le match ID
}
close $tsv;

# Lire le fichier fasta et extraire les séquences qui matchent
open(my $in_fasta, '<', $fasta_file);


my $write = 0;
while (my $line = <$in_fasta>) {
    if ($line =~ /^>(\S+)/) {
        my $id = $1;
        $write = exists $wanted_ids{$id};
    }

    print $line if $write;
}

close $in_fasta;
