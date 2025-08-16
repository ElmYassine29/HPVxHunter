#!/bin/bash
set -euo pipefail

#  ======= HELP SECTION =======
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo "Usage: $0 [options]

Available options:

  --db <pave|genbank|path_to_custom_db>  : Reference database to use for analysis.
                                          By default, the PAVE HPV reference DB (data/DB_full.fasta) will be used.
                                          You can specify 'genbank' to use the GenBank database,
                                          or provide the path to a custom database file in FASTA format.

  --sample <fasta_file>                  : Input sample sequences in FASTA format.
                                          Default: test/test.fasta

  --run_name <name>                      : Name of the run/session. Used for naming outputs.
                                          If not specified, a timestamp-based name will be used.

  --debug 0|1                            : If set to 1, keeps intermediate files for debugging.
                                          Default: 0.

  --deep_analysis 0|1                    : If set to 1, runs extra steps for divergent or novel HPV types.
                                          Default: 0. (only available when using the PAVE database as reference)

  --format <compact|detailed>            : Specifies the format of the output Excel file if genbank or custom db analysis.
                                          'compact' will generate a simpler output, while 'detailed' will include more extensive
                                          information. Default: 'compact'.

  --plot <pdf|png>                       : Specifies the format of the plot (pdf or png).
                                          Default is 'png'.
                                          

  -h, --help                             : Show this help message and exit.

"
    exit 0
fi




# ======= DEFAULT VALUES =======
SCRIPT_DIR="$(dirname "$0")"

DB="$SCRIPT_DIR/data/DB_full.fasta"                      # default: PAVE reference
SAMPLE_FILE="$SCRIPT_DIR/test/test.fasta"
SESSION=""
DEBUG=0
DEEP_ANALYSIS=0
FORMAT="compact"
PLOT='png'

# ======= PARSE OPTIONS =======
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)
      if [[ -z "$2" ]]; then
        echo "Pave database used by default"
        DB="$SCRIPT_DIR/data/DB_full.fasta"
        shift 1
      elif [[ "${2,,}" == "pave" ]]; then
        DB="$SCRIPT_DIR/data/DB_full.fasta"
        shift 2
      elif [[ "${2,,}" == "genbank" ]]; then
        DB="genbank"
        shift 2
      else
        if [[ -f "$2" ]]; then
          DB="$2"
          shift 2
        else
          echo "Error: Database file '$2' does not exist."
          exit 1
        fi
      fi
      ;;
      
    --sample)
      if [[ -f "$2" ]]; then
        SAMPLE_FILE="$2"
        shift 2
      else
        echo "Error: Sample file '$2' does not exist."
        exit 1
      fi
      ;;
      
    --run_name)
      SESSION=$(echo "$2" | tr -cd '[:alnum:]_-')
      RANDOM_SUFFIX=$(head /dev/urandom | tr -dc a-z0-9 | head -c6)
      SESSION="${SESSION}_${RANDOM_SUFFIX}"
      shift 2
      ;;
      
    --debug)
      if [[ "$2" == "0" || "$2" == "1" ]]; then
        DEBUG="$2"
        shift 2
      else
        echo "Invalid value for --debug: $2"
        echo "Allowed values are 0 or 1."
        exit 1
      fi
      ;;
      
    --deep_analysis)
      if [[ "$2" == "0" || "$2" == "1" ]]; then
        DEEP_ANALYSIS="$2"
        shift 2
      else
        echo "Invalid value for --deep_analysis: $2"
        echo "Allowed values are 0 or 1."
        exit 1
      fi
      ;;
      
    --format)
      if [[ "${2,,}" == "compact" || "${2,,}" == "detailed" ]]; then
        FORMAT="${2,,}"  # Accepts "compact" or "detailed"
        shift 2
      else
        echo "Invalid format option: $2"
        echo "Valid formats: compact, detailed."
        exit 1
      fi
      ;;
      
    --plot)
      if [[ "${2,,}" == "png" || "${2,,}" == "pdf" ]]; then
        PLOT="${2,,}"  # Accepts "png" or "pdf"
        shift 2
      else
        echo "Invalid format option: $2"
        echo "Valid formats: png, pdf."
        exit 1
      fi
      ;;
      
    *)
      echo "Unknown option: $1"
      echo "Use --help to see available options."
      exit 1
      ;;
  esac
done





# ======= CHECK IF DB IS 'genbank' AND HANDLE DOWNLOAD =======

if [[ "$DB" == "genbank" ]]; then
    echo "[INFO] Option --db genbank détectée."

    META="$SCRIPT_DIR/data/genbank_db.meta"
    GENBANK_DB="$SCRIPT_DIR/data/genbank_db.fasta"

    count_remote=$(python3 "$SCRIPT_DIR/scripts/genbank_check_count.py" | tr -d '[:space:]')
    echo "[DEBUG] Remote GenBank count (from Entrez): $count_remote"

    if [[ -f "$META" ]]; then
        count_local=$(grep "^count=" "$META" | cut -d= -f2)
        echo "[DEBUG] count local est : $count_local"
    else
        count_local=0
    fi

    if [[ "$count_remote" -eq "$count_local" ]]; then
        echo "[INFO] GenBank DB is up to date. Skipping download."
    else
        echo "[INFO] New sequences available. Downloading..."
        python3 "$SCRIPT_DIR/scripts/download_genbank.py" --output "$GENBANK_DB"
        echo "[INFO] New database downloaded and saved to: $GENBANK_DB"

        # Update metadata
        echo "count=$count_remote" > "$META"
        echo "date=$(date -Iseconds)" >> "$META"
    fi

    # Use downloaded GenBank DB for analysis
    DB="$GENBANK_DB"
fi






# ======= GENERATE SESSION NAME IF NEEDED =======
SESSION="${SESSION:-run_$(date +%Y%m%d_%H%M%S)_$(head /dev/urandom | tr -dc a-z0-9 | head -c6)}"

# ======= OUTPUT SETUP =======

mkdir -p $SCRIPT_DIR/results/"$SESSION"

DIR="$SCRIPT_DIR/results/"$SESSION""

mkdir -p $DIR/intermediate_results
mkdir -p $DIR/plots
mkdir -p $DIR/final_files

OUT_CSV="results_${SESSION}.tsv"

# ======= ALIGNMENT =======
vsearch --usearch_global "$SAMPLE_FILE" --db "$DB" --id 0  --maxaccepts 1 --threads 8 --blast6out  $DIR/intermediate_results/sample.besthit.tsv

# ======= MAIN ANALYSIS =======

if [[ "$DB" == "$SCRIPT_DIR/data/DB_full.fasta" ]]; then
    echo "[INFO] Main Analysis Based on Database PAVE"
    # Ecrire les réslutats dans le fichier tsv + ecrire les probables nouveaux types dans le fichier $SCRIPT_DDIR/analyze_results/probable_newtypes.txt
    perl "$SCRIPT_DIR/scripts/find_clstr.pl" $DIR/intermediate_results/sample.besthit.tsv "$DIR/intermediate_results" "$DB" "$SAMPLE_FILE" > "$DIR/intermediate_results/$OUT_CSV"

    #Make the Excel file 
    perl "$SCRIPT_DIR/scripts/pave_excel.pl" "$DIR/intermediate_results/$OUT_CSV" "$DIR/final_files" "$SESSION"
    echo "[INFO] XLSX file ready: $DIR/final_files/"


elif [[ "$DB" == "$SCRIPT_DIR/data/genbank_db.fasta" ]]; then
    echo "[INFO] Main Analysis Based on Database GenBank"
    # Analyze on genbank
    perl "$SCRIPT_DIR/scripts/find_genbank_2.pl" "$DIR/intermediate_results/sample.besthit.tsv" "$DIR/intermediate_results" "$DB" "$SAMPLE_FILE"  > "$DIR/intermediate_results/$OUT_CSV"

    #Extract sequences of IDs matched with our samples and make a blast 
    perl "$SCRIPT_DIR/scripts/extract_genbank_matches.pl" "$DIR/intermediate_results/$OUT_CSV" "$DB" >"$DIR/intermediate_results/genbank_matches.fasta"
    vsearch --usearch_global "$DIR/intermediate_results/genbank_matches.fasta" --db "$SCRIPT_DIR/data/DB_full.fasta" --id 0 --maxaccepts 1 --threads 8 --blast6out $DIR/intermediate_results/besthit_genbank_match.tsv

    #Join both analyzes and Make xlsx file
    echo "[INFO] Starting XLSX generation..."
    perl "$SCRIPT_DIR/scripts/final_genbank.pl" "$DIR/intermediate_results/$OUT_CSV" "$DIR/intermediate_results/besthit_genbank_match.tsv" "$DIR/final_files" "$FORMAT" "$SESSION"
    echo "[INFO] XLSX file ready: $DIR/final_files/"


else
    echo "[INFO] Main Analysis Based on Custom Database"
    # Script pour une DB personnalisée
    perl "$SCRIPT_DIR/scripts/find_genbank_2.pl" "$DIR/intermediate_results/sample.besthit.tsv" "$DIR/intermediate_results" "$DB" "$SAMPLE_FILE"  > "$DIR/intermediate_results/$OUT_CSV"
    
    #Extract sequences of IDs matched with our samples and make a blast 
    perl "$SCRIPT_DIR/scripts/extract_genbank_matches.pl" "$DIR/intermediate_results/$OUT_CSV" "$DB" >"$DIR/intermediate_results/genbank_matches.fasta"
    vsearch --usearch_global "$DIR/intermediate_results/genbank_matches.fasta" --db "$SCRIPT_DIR/data/DB_full.fasta" --id 0 --maxaccepts 1 --threads 8 --blast6out $DIR/intermediate_results/besthit_genbank_match.tsv

    #Join both analyzes and Make xlsx file
    echo "[INFO] Starting XLSX generation..."
    perl "$SCRIPT_DIR/scripts/final_genbank.pl" "$DIR/intermediate_results/$OUT_CSV" "$DIR/intermediate_results/besthit_genbank_match.tsv" "$DIR/final_files" "$FORMAT" "$SESSION"
    echo "[INFO] XLSX file ready: $DIR/final_files/"
fi




# ======= OPTIONAL DEEP ANALYSIS =======
if [ "$DEEP_ANALYSIS" -eq 1 ] && [ "$DB" = "$SCRIPT_DIR/data/DB_full.fasta" ]; then
    if [ -f "$DIR/intermediate_results/probable_newtypes.txt" ]; then
        echo "[INFO] Deep analysis of divergent types..."
        perl $SCRIPT_DIR/scripts/catch_new.pl "$SAMPLE_FILE" "$DIR/intermediate_results/probable_newtypes.txt" > $DIR/intermediate_results/prob_newtypes.fasta

        mkdir -p $DIR/intermediate_results/puma_results
        python3 $SCRIPT_DIR/scripts/puma/scripts/run_puma.py -i "$DIR/intermediate_results/prob_newtypes.fasta" -d "$SCRIPT_DIR/scripts/puma/data_dir/" -o "$DIR/intermediate_results/puma_results/"

        perl $SCRIPT_DIR/scripts/cut_L1.pl "$DIR/intermediate_results/puma_results/" > $DIR/intermediate_results/${SESSION}_samples_L1.fasta

        vsearch --usearch_global $DIR/intermediate_results/${SESSION}_samples_L1.fasta \
                --db $SCRIPT_DIR/data/DB_full_L1.fasta \
                --id 0 --maxaccepts 1 --threads 8 \
                --blast6out $DIR/intermediate_results/${SESSION}_sample_L1.besthit.tsv

        perl $SCRIPT_DIR/scripts/find_newtype.pl $DIR/intermediate_results/${SESSION}_sample_L1.besthit.tsv > $DIR/intermediate_results/${SESSION}_result_deep_analysis.tsv
        ssconvert $DIR/intermediate_results/${SESSION}_result_deep_analysis.tsv  $DIR/final_files/${SESSION}_result_deep_analysis.xlsx

    else
        echo "[INFO] There are no potential new HPV types."
    fi
else
    echo "[INFO] Deep analysis is only available when using the PAVE reference database."
fi



# ======= FINAL PLOT =======

if [[ "$DB" == "$SCRIPT_DIR/data/DB_full.fasta" ]]; then
    if Rscript "$SCRIPT_DIR/scripts/plot_hpvs.R" "$DIR/final_files/${SESSION}_results.xlsx" "$SESSION" "$DIR/plots/" "$PLOT"; then
        echo "[INFO] Plots generated: Distribution of Types, Lineages and Sub-Lineages"
    else
        echo "[WARNING] Plot generation failed."
    fi
else
    if Rscript "$SCRIPT_DIR/scripts/plot_hpvs_genbank.R" "$DIR/final_files/${SESSION}_results.xlsx" "$SESSION" "$DIR/plots/" "$PLOT" "$FORMAT" ; then
        echo "[INFO] Plots generated: Distribution of Types, Lineages and Sub-Lineages"
    else
        echo "[WARNING] Plot generation failed."
    fi
fi



# ======= CLEANUP / DEBUG =======
if [ "$DEBUG" -eq 0 ]; then
    rm -rf $DIR/intermediate_results/*
     
else
    echo "[DEBUG] Kept intermediate files in : $DIR/intermediate_results/"
fi

