#!/usr/bin/perl
use strict;
use warnings;
use autodie;

# Input is the results file from vsearch
my ($file, $output_dir,$DB,$fasta) = @ARGV;

my @prob_newtypes;


my %fasta_lengths = %{ get_all_lengths($fasta) };
my %db_lengths    = %{ get_all_lengths($DB) };

# Header of output (to STDOUT)
print join("\t", "Sequence", "Match", "Type", "Lineage", "Sub-lineage", "Similarity", "Length Difference", "Results"), "\n";

# Open the input file
open(my $fh, '<', $file);

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;  # skip empty lines

    my @columns = split("\t", $line);

    if (@columns < 3) {
        warn "Ligne mal formatée (moins de 3 colonnes): $line\n";
        next;
    }

    my $seq_id = $columns[0];
    my $match = $columns[1] // '';
    my $perc  = $columns[2] // 0;

    process_cluster($seq_id, $match, $perc);

    if ($perc < 95) {
        push @prob_newtypes, $seq_id;
    }
}

if (@prob_newtypes) {
    open(my $out, '>', "$output_dir/probable_newtypes.txt");
    print $out "$_\n" for @prob_newtypes;
    close $out;
}

sub process_cluster {
    my ($seq, $hit, $perc) = @_;
    my ($type, $lineage, $sub_lineage, $res) = classify_similarity($perc);


   

    my $seq_value = $fasta_lengths{$seq} // 0;

    my $hit_value = 0;
    for my $header (keys %db_lengths) {
        if ($header =~ /\Q$hit\E/) {  # \Q...\E échappe les caractères spéciaux dans $hit
            $hit_value = $db_lengths{$header};
            last;  # On arrête après la première correspondance
        }
    }

    my $diff = abs($seq_value - $hit_value);

    print join("\t", $seq, $hit, $type, $lineage, $sub_lineage, $perc, $diff, $res), "\n";
}

sub classify_similarity {
    my ($perc) = @_;
    

    my ($type, $lineage, $sub_lineage, $res) = ('ND', 'ND', 'ND', "--deep_analysis 1 recommended - no strong match ($perc%)");

    if ($perc > 99.5) {
        ($type, $lineage, $sub_lineage, $res) = ('Same type', 'Same lineage', 'Same sublineage', 'Perfect or near-perfect match');
    }
    elsif ($perc > 98.5) {
        ($type, $lineage, $sub_lineage, $res) = ('Same type', 'Different lineage', 'Different sublineage', 'Close match - minor divergence');
    }
    elsif ($perc > 98.0) {
        ($type, $lineage, $sub_lineage, $res) = ('Same type', 'Different lineage', 'Different sublineage', 'Minor divergent - review recommended');
    }
    elsif ($perc > 95.0) {
        ($type, $lineage, $sub_lineage, $res) = ('Same type', 'Different lineage', 'Different sublineage', 'Moderately divergent - review recommended');
    }
    elsif ($perc > 90.0) {
        ($type, $lineage, $sub_lineage, $res) = ('Same type', 'Different lineage', 'Different sublineage', 'Highly divergent - review recommended');
    }

    return ($type, $lineage, $sub_lineage, $res);
}




sub get_all_lengths {
    my ($fasta_file) = @_;
    my %lengths;

    open my $fh, '<', $fasta_file or die "Cannot open $fasta_file: $!";

    my $header = '';
    my $seq = '';

    while (<$fh>) {
        chomp;
        next if /^$/;

        if (/^>(.*)/) {
            if ($header) {
                $lengths{$header} = length($seq);
            }
            $header = $1;
            $seq = '';
        } else {
            s/\s+//g;
            $seq .= $_;
        }
    }

    # Dernière séquence
    if ($header) {
        $lengths{$header} = length($seq);
    }

    close $fh;
    return \%lengths;
}