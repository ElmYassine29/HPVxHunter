from Bio import Entrez
Entrez.email = "your_email@example.com"

search_term = '(("Human papillomavirus"[Organism] OR "Human papillomavirus"[All Fields] OR HPV[All Fields]) ' \
              'AND 7000:10000[Sequence Length]) AND viruses[filter]'

handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=0)
record = Entrez.read(handle)
print(int(record["Count"]))
