#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use Excel::Writer::XLSX;

# Arguments d'entrée
my ($file, $out_dir, $session) = @ARGV;

# Liste des types HPV à haut risque
my @hr_hpvs = qw(35 31 16 33 58 52 18 45 59 39 68 56 66 51);
my %hr_map = map { ("HPV$_" => 1) } @hr_hpvs;

# Création du fichier Excel
my $xlsx_file = "$out_dir/${session}_results.xlsx";
my $workbook  = Excel::Writer::XLSX->new($xlsx_file);
my $worksheet = $workbook->add_worksheet();

# Définir les formats
my $header_fmt = $workbook->add_format(bg_color => 'yellow', bold => 1);
my $hr_fmt     = $workbook->add_format(bg_color => '#FFC7CE');  # rouge clair
my $normal_fmt = $workbook->add_format();

# Écriture de l'en-tête dans le fichier Excel
my @headers = ("Sequence", "Type", "Lineage", "Sub-lineage", "Similarity", "Length Difference", "Results");
my $row = 0;
for my $col (0..$#headers) {
    $worksheet->write($row, $col, $headers[$col], $header_fmt);
}

# Traitement du fichier d'entrée et écriture dans le fichier Excel
$row++;
open(my $fh, '<', $file);
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;  # Skip empty lines
    next if $line =~ /Sequence/;
    
    my @cols = split("\t", $line);
    
    # Vérifie qu'on a au moins 3 colonnes
    if (@cols < 3) {
        warn "Ligne mal formatée (moins de 3 colonnes): $line\n";
        next;
    }

    my $seq_id = $cols[0];
    my $hit = $cols[1] // '';
    my $perc = $cols[2] // 0;
    
    # Identifie le type HPV
    my $type = 'ND';
    if (defined $hit && $hit =~ /(HPV\d+)/) {
        $type = $1;
    }

    # Vérifie si c'est un HPV à haut risque
    my $is_hr = exists $hr_map{$type} ? 1 : 0;

    # Applique les formats conditionnels
    for my $col (0..$#cols) {
        if ($is_hr) {
            $worksheet->write($row, $col, $cols[$col], $hr_fmt);
        } else {
            $worksheet->write($row, $col, $cols[$col], $normal_fmt);
        }
    }
    $row++;
}

# Fermeture du fichier
close $fh;

# Finalisation du fichier Excel
$workbook->close();

print "Fichier Excel créé : $xlsx_file\n";
