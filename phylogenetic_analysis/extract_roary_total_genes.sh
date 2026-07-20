#!/bin/bash

BASE_DIR="/mmfs1/projects/barney.geddes/tania.gupta/TG_comparitive_genomics/whole_genome/roary_array_final_assembly_185_st_TG"
OUT="roary_total_gene_clusters_by_blastp.csv"

echo "blastp_identity,total_gene_clusters" > "$OUT"

for DIR in "$BASE_DIR"/*; do
    [[ -d "$DIR" ]] || continue

    BASENAME=$(basename "$DIR")

    # Only use folders that are purely numeric, e.g. 60, 61, 62...
    if [[ "$BASENAME" =~ ^[0-9]+$ ]]; then
        FILE="$DIR/summary_statistics.txt"

        if [[ -f "$FILE" ]]; then
            TOTAL=$(awk -F'\t' '$1 == "Total genes" {print $3}' "$FILE")

            if [[ -n "$TOTAL" ]]; then
                echo "${BASENAME},${TOTAL}" >> "$OUT"
            else
                echo "WARNING: Could not find Total genes in $FILE" >&2
            fi
        else
            echo "WARNING: Missing $FILE" >&2
        fi
    fi
done

sort -t, -k1,1n "$OUT" -o "$OUT"

echo "Done. Output saved to:"
echo "$OUT"
cat "$OUT"



