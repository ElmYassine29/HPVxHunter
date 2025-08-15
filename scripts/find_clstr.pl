#!/usr/bin/perl
use strict;
use warnings;
use autodie;


# Input is the results file from v search
my ($file, $output_dir,$DB,$fasta) = @ARGV;

my @prob_newtypes;


my %fasta_lengths = %{ get_all_lengths($fasta) };
my %db_lengths    = %{ get_all_lengths($DB) };

# Open the file from vsearch...
open(my $fh, '<', $file);
print join("\t", "Sequence", "Type", "Lineage", "Sub-lineage", "Similarity", "Length Difference", "Results"), "\n";


while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;  # skip empty lines

    my @columns = split("\t", $line);
    
    # Vérifie qu'on a au moins 3 colonnes
    if (@columns < 3) {
        warn "Ligne mal formatée (moins de 3 colonnes): $line\n";
        next;
    }

    my $seq_id = $columns[0];
    my $match = $columns[1] // '';  # définit une chaîne vide si undef
    my $perc = $columns[2] // 0;    # définit 0 si undef

    process_cluster($seq_id, $match, $perc);

    if($perc < 95){
        push @prob_newtypes, $seq_id;

    }


}

if (@prob_newtypes) {
    open(my $out, '>', "$output_dir/probable_newtypes.txt");
    for my $seq_id (@prob_newtypes) {
        print $out "$seq_id\n";
    }
    close $out;
}

sub process_cluster {
    my ($seq, $hit, $perc) = @_;
    
    # Initialisation des variables pour éviter "uninitialized" warnings
    my ($type, $lineage, $sub_lineage) = ('ND', 'ND', 'ND');
    
    # PREMIER CAS: si REF
    if (defined $hit && $hit =~ /(HPV\d+)REF/) {
        $type = $1; # exemple HPV68
        ($lineage, $sub_lineage) = ('A', 'A-1'); # ref toujours le A-1
    }
    # si c'est un variant
    elsif (defined $hit && $hit =~ /^HPV(\d+)-([A-Z])(\d*)/i) {
        $type = "HPV$1";
        $lineage = $2;
        $sub_lineage = $3 ? "$2-$3" : "$2-1";
    }
    
    my ($ls, $sls, $res) = classify_similarity($perc, $lineage, $sub_lineage,$type);
    
   if ($perc < 90.0) {
    $type = 'ND';
    }
    
    my $seq_value = 0;

if (exists $fasta_lengths{$seq}) {
    $seq_value = $fasta_lengths{$seq};
} else {
    for my $header (keys %fasta_lengths) {
        if ($header =~ /\b\Q$seq\E\b/) {
            $seq_value = $fasta_lengths{$header};
            last;
        }
    }

    
}

    my $hit_value = 0;
    

    if (exists $db_lengths{$hit}) {
    $hit_value = $db_lengths{$hit};
} else {
    # Sinon, chercher un header contenant le hit comme sous-chaîne
    for my $header (keys %db_lengths) {
        if ($header =~ /\b\Q$hit\E\b/) {  # \b pour délimitation de mot
            $hit_value = $db_lengths{$header};
            last;
        }
    }

    
}

    my $diff = abs($seq_value - $hit_value);
    
    print join("\t", $seq, $type, $ls, $sls, $perc, $diff, $res), "\n";
}

sub classify_similarity {
    my ($perc, $ref_lineage, $ref_sub_lineage, $type_ref) = @_;
    
    # Initialisation des valeurs par défaut
    $ref_lineage //= 'ND';
    $ref_sub_lineage //= 'ND';
    
    return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, the closest type is $type_ref") if $perc <= 85;
    
    if ($perc > 99.5) {
        return ($ref_lineage, $ref_sub_lineage, 'Same sublineage');
    }
    elsif ($perc > 99.0) {
        return ($ref_lineage, 'ND', 'Probably new sublineage');
    }
    elsif ($perc > 98.0) {
        return ('ND', 'ND', 'Probably new lineage - Minor divergent');
    }
    elsif ($perc > 95.0) {
        return ('ND', 'ND', 'Probably new lineage - Moderately divergent');
    }
    elsif ($perc > 90.0) {
        return ('ND', 'ND', 'Probably new lineage - Highly divergent');
    }
    else {
        return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, the closest type is $type_ref");
    }
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