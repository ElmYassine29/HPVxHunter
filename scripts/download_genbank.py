import argparse
from Bio import Entrez
import time
import datetime

# Argument parser to handle output path
parser = argparse.ArgumentParser(description="Download complete HPV genomes from GenBank.")
parser.add_argument("--output", required=True, help="Path to output FASTA file.")
args = parser.parse_args()

# Your email for NCBI Entrez
Entrez.email = "ton_email@example.com"  # Replace with your real email

# GenBank search term optimized for complete HPV genomes
search_term = '(("Human papillomavirus"[Organism] OR "Human papillomavirus"[All Fields] OR HPV[All Fields]) ' \
              'AND 7000:10000[Sequence Length]) AND viruses[filter]'

# First, get the number of sequences matching the query
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=0)
record = Entrez.read(handle)
total_sequences = int(record["Count"])
handle.close()

print(f"[INFO] Total number of sequences found: {total_sequences}")

# Fetch all IDs
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=total_sequences)
record = Entrez.read(handle)
id_list = record["IdList"]
handle.close()

print(f"[INFO] Number of IDs retrieved: {len(id_list)}")

# Download sequences in batches
batch_size = 500
fasta_data = ""

for start in range(0, len(id_list), batch_size):
    end = min(start + batch_size, len(id_list))
    batch_ids = id_list[start:end]
    
    try:
        print(f"[INFO] Downloading sequences {start + 1} to {end}...")
        handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="fasta", retmode="text")
        fasta_data += handle.read()
        handle.close()
    except Exception as e:
        print(f"[ERROR] Problem with batch {start}-{end}: {str(e)}")
    
    time.sleep(1)  # Respect NCBI usage limits

# Save FASTA file
with open(args.output, "w") as f:
    f.write(fasta_data)

nb_sequences = fasta_data.count(">")

# Save metadata file
meta_path = args.output + ".meta"
with open(meta_path, "w") as meta_file:
    meta_file.write(f"count={nb_sequences}\n")
    meta_file.write(f"date={datetime.datetime.now().isoformat()}\n")

print(f"[INFO] {nb_sequences} sequences saved to {args.output}")
print(f"[INFO] Metadata saved to {meta_path}")
