#!/usr/bin/perl
use strict;
use warnings;
use autodie;

# Input is the results file from v search
my ($file) = shift or die "Usage: $0 <vsearch_results_file>\n";


# Open the file from vsearch...
open(my $fh, '<', $file);
print join("\t", "Potential New Type Sequence", "Final Assignment"), "\n";

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

    if($perc<90){

        print"$seq_id\tis a new HPV type\n";
    }

	else{ print"$seq_id\tis NOT a new HPV type\n";}



}
