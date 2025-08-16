# HPVxHunter

**HPVxHunter** is a command-line pipeline for the **classification and genomic analysis** of Human Papillomavirus (HPV) sequences.  
It identifies HPV types, lineages, and sublineages, and can flag potential novel HPV types following established classification rules.  
The tool supports curated, public, and user-provided reference databases, with **automatic GenBank database updating**.

---

##  Main Features

- **Multiple database modes**  
  - **PaVE** (Papillomavirus Episteme) – curated HPV reference database with complete lineage and sublineage annotation.  
  - **GenBank** – automatically downloaded and updated via NCBI Entrez; broader coverage but may lack complete annotation.  
  - **Custom database** – any user-provided FASTA file (e.g., new or unpublished HPV sequences).  

- **Automatic GenBank DB update**  
  - Checks NCBI for new sequences and updates the local GenBank database before each run.  

- **Adaptive pipeline logic**  
  - **Direct classification** with PaVE (type, lineage, sublineage, high-risk flagging).  
  - **Two-step classification** for GenBank/custom DB:  
    1. Identify best match in chosen DB.  
    2. Map the match against PaVE for improved annotation.  

- **High-quality outputs**  
  - Excel reports with classification results and length-based false-positive filtering.  
  - Automatic highlighting of **high-risk HPV types**.  
  - Graphical summaries:  
    - Top 20 most frequent HPV types  
    - Lineage distribution for top 12 HPV types  
    - Sublineage distribution per lineage for top 12 HPV types  

- **Novel type detection** *(optional)*  
  - Triggered with `--deep_analysis 1` when using the PaVE database.  
  - Extracts the **L1 ORF** with **PuMA** and compares it against PaVE's L1 references.  
  - Flags sequences with >10% L1 divergence as potential new HPV types.  

- **Clean & reproducible execution**  
  - Intermediate files automatically removed unless `--debug 1` is set.  
  - Organized outputs in a dedicated `results/<run_name>` directory.  

---

## ⚙ How It Works

1. **Input parsing & configuration**  
   - User specifies reference DB (`--db pave | genbank | <custom_path>`), sample FASTA, and optional parameters (`--deep_analysis`, `--format`, `--plot`, etc.).  

2. **Database preparation**  
   - **PaVE**: uses included curated reference.  
   - **GenBank**: automatically checks NCBI for new sequences; downloads and replaces local DB if updates are found.  
   - **Custom**: uses provided FASTA file directly.  

3. **Global alignment**  
   - Uses **VSEARCH** for pairwise similarity search against the selected DB.  

4. **Primary analysis**  
   - **PaVE mode**: Classification + annotation, generation of Excel file.  
   - **GenBank/Custom**:  
     - Identify best hits.  
     - Map those hits to PaVE to infer type, lineage, and sublineage.  

5. **Optional deep analysis** (`--deep_analysis 1`, PaVE only)  
   - Detects divergent sequences.  
   - Extracts L1 ORF with PuMA.  
   - Compares L1 sequences against PaVE L1 database.  
   - Generates a dedicated Excel report for potential novel HPV types.  

6. **Visualization**  
   - Generates barplots and lineage/sublineage distributions.  

7. **Final output**  
   - All results saved in `results/<run_name>` including:  
     - Excel reports  
     - Graphical plots  
     - (Optional) deep analysis results  

---

##  Installation

We recommend installing HPVxHunter in a **conda environment**:

```bash
git clone https://github.com/ElmYassine29/HPVxHunter.git
cd HPVxHunter
conda env create -f env/environment.yml
conda activate HPVxHunter

sudo apt update
sudo apt install gnumeric
```


## Usage 

Make sure your input sequences and any custom database are in FASTA format.

Basic example using PaVE database and test fasta file:

```bash
./HPVxHunter.sh --db pave --sample test/test.fasta --run_name example_run

```

For full option details, run:
```bash
./HPVxHunter.sh --help

```

If you use HPVxHunter in your research, please cite the associated paper (in preparation).
