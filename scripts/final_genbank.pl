#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use Excel::Writer::XLSX;

# Get arguments
my ($tsv_one, $tsv_two, $out_dir,$format,$session) = @ARGV;
die "Usage: $0 <main.tsv> <pave.tsv> <out_dir>\n" unless $tsv_one && $tsv_two && $out_dir;

my %pave;
my @hr_hpvs = qw(35 31 16 33 58 52 18 45 59 39 68 56 66 51); # high-risk HPV types
my %hr_map = map { ("HPV$_" => 1) } @hr_hpvs;

# Step 1: Parse VSEARCH result to build %pave hash
open(my $fh, '<', $tsv_two);
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @columns = split("\t", $line);
    next unless @columns >= 3;

    my ($seq_id, $match, $perc) = @columns[0,1,2];
    $match //= '';
    $perc //= 0;

    if($format eq 'compact'){
    my $res_pave = process_cluster($seq_id, $match, $perc);
    $pave{$seq_id} = $res_pave;
    }
    else{
    my $res_pave = process_cluster_detailed($seq_id, $match, $perc); # all columns in one argument joined by "-"
    $pave{$seq_id} = $res_pave;

    }
}
close $fh;

# Step 2: Prepare Excel output
my $xlsx_file = "$out_dir/${session}_results.xlsx";
my $workbook  = Excel::Writer::XLSX->new($xlsx_file);
my $worksheet = $workbook->add_worksheet();

# Formats
my $header_fmt = $workbook->add_format(bg_color => 'yellow', bold => 1);
my $hr_fmt     = $workbook->add_format(bg_color => '#FFC7CE');  # light red
my $normal_fmt = $workbook->add_format();

# Step 3: Process TSV and write to Excel
open(my $tsv1, '<', $tsv_one);
my $row = 0;

while (<$tsv1>) {
    chomp;
    my @cols = split(/\t/);
    
    if ($cols[0] eq "Sequence") {
        
        if ($format eq "detailed") {
    push @cols, (
        "Sequence blasted in PAVE",
        "Type",
        "Lineage",
        "Sub-lineage",
        "Percentage",
        "Length Difference",
        "PAVE result"
        );
        } else {
        push @cols, "PAVE match";  # Format compact
            }

        for my $col (0..$#cols) {
            $worksheet->write($row, $col, $cols[$col], $header_fmt);
        }
    } else {

        my $paver;

        if ($format eq "compact") {
        my $match_id = $cols[1];
         $paver    = $pave{$match_id} // "NA";
        push @cols, $paver;
        }
        else{
            my $match_id = $cols[1];
             $paver = $pave{$match_id} // "NA";
            if ($paver eq "NA") {
    push @cols, ("NA") x 6;  # 6 fois NA pour les 6 colonnes
} else {
    my @parts = split /-/, $paver;
    push @cols, @parts;
}

        }

        if ($format eq "compact") {
        for my $col (0..$#cols) {
            if ($col == $#cols && $paver =~ /(HPV\d+)/ && exists $hr_map{$1}) {
                $worksheet->write($row, $col, $cols[$col], $hr_fmt);
            } else {
                $worksheet->write($row, $col, $cols[$col], $normal_fmt);
            }
        }
        }
        else{
            for my $col (0..$#cols) {
    # Si c'est la 10e colonne (index 9) et que le motif est trouvé
    if ($col == 9 && $paver =~ /(HPV\d+)/ && exists $hr_map{$1}) {
        # Applique le format rouge à la 10e colonne
        $worksheet->write($row, $col, $cols[$col], $hr_fmt);
    } else {
        # Applique le format normal pour les autres
        $worksheet->write($row, $col, $cols[$col], $normal_fmt);
    }
}
        }
    }
    $row++;
}
close $tsv1;

# Subroutines

#if compact
sub process_cluster {
    my ($seq, $hit, $perc) = @_;
    my ($type, $lineage, $sub_lineage) = ('ND', 'ND', 'ND');

    if ($hit =~ /(HPV\d+)REF/) {
        $type = $1;
        ($lineage, $sub_lineage) = ('A', 'A1');
    }
    elsif ($hit =~ /^HPV(\d+)-([A-Z])(\d*)/i) {
        $type = "HPV$1";
        $lineage = $2;
        $sub_lineage = defined($3) && $3 ne '' ? "$2$3" : "$21";
    }

    my ($ls, $sls, $res) = classify_similarity_detailed($perc, $lineage, $sub_lineage, $type);
    return "$type-$ls-$sls";
}

sub classify_similarity {
    my ($perc, $ref_lineage, $ref_sub_lineage, $type_ref) = @_;
    $ref_lineage //= 'ND';
    $ref_sub_lineage //= 'ND';

    return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, closest type: $type_ref") if $perc <= 85;
    return ($ref_lineage, $ref_sub_lineage, 'Same sublineage') if $perc > 99.5;
    return ($ref_lineage, 'ND', 'Probably new sublineage') if $perc > 98.5;
    return ($ref_lineage, 'ND', 'Probably new lineage') if $perc > 98.0;
    return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, closest type: $type_ref");
}



#if detailed 

sub process_cluster_detailed {
    my ($seq, $hit, $perc) = @_;
    
    # Initialisation des variables pour éviter "uninitialized" warnings
    my ($type, $lineage, $sub_lineage) = ('ND', 'ND', 'ND');
    
    # PREMIER CAS: si REF
    if (defined $hit && $hit =~ /(HPV\d+)REF/) {
        $type = $1; # exemple HPV68
        ($lineage, $sub_lineage) = ('A', 'A1'); # ref toujours le A-1
    }
    # si c'est un variant
    elsif (defined $hit && $hit =~ /^HPV(\d+)-([A-Z])(\d*)/i) {
        $type = "HPV$1";
        $lineage = $2;
        $sub_lineage = defined($3) && $3 ne '' ? "$2$3" : "$21";
    }
    
    my ($ls, $sls, $res) = classify_similarity($perc, $lineage, $sub_lineage,$type);
    my $diff = sprintf("%.2f", 100 - $perc);
    return("$seq-$type-$ls-$sls-$perc-$diff-$res");
}

sub classify_similarity_detailed {
    my ($perc, $ref_lineage, $ref_sub_lineage, $type_ref) = @_;
    
    # Initialisation des valeurs par défaut
    $ref_lineage //= 'ND';
    $ref_sub_lineage //= 'ND';
    
    return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, the closest type is $type_ref") if $perc <= 85;
    
    if ($perc > 99.5) {
        return ($ref_lineage, $ref_sub_lineage, 'Same sublineage');
    }
    elsif ($perc > 98.5) {
        return ($ref_lineage, 'ND', 'Probably new sublineage');
    }
    elsif ($perc > 98.0) {
        return ($ref_lineage, 'ND', 'Probably new lineage');
    }
    else {
        return ('ND', 'ND', "No significant match ($perc), try --deep_analysis, the closest type is $type_ref");
    }
}
