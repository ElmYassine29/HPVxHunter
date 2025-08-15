use strict;
use warnings;
use File::Find;

my $directory = shift;

my $ORF = "L1";


sub process_file {
    my $file = $_;
    return unless $file =~ /\.csv$/;  # Vérifie que c'est un fichier .csv

    open my $fh, '<', $file or die "Impossible d'ouvrir $file: $!";
    (my $filename = $file) =~ s/\.csv$//;  # Supprime l'extension .csv

    while (my $line = <$fh>) {
        chomp $line;

        my @cols = split(/,/, $line);

        # Vérifie qu'il y a au moins 4 colonnes et que la 2e est "L1"
        if (scalar(@cols) >= 4 && $cols[1] eq "$ORF") {
            #print $out ">$filename\n$cols[3]\n";
            print  ">$filename\n$cols[3]\n";  # FASTA format
            last;  # Une seule séquence "L1" par fichier
        }
    }

    close $fh;
}

find(\&process_file, $directory);

#print "All L1 ORFs are done\n";
