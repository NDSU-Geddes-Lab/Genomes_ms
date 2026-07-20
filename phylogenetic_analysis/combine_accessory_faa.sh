#!/bin/bash



PROKKA_DIR="/mmfs1/projects/barney.geddes/tania.gupta/TG_comparitive_genomics/accessory/prokka_out_no_distinct_clade"



OUTDIR="/mmfs1/projects/barney.geddes/tania.gupta/TG_comparitive_genomics/accessory/prokka_out_no_distinct_clade/accessory_functional_categories_clean"



mkdir -p "$OUTDIR"



COMBINED="$OUTDIR/all_accessory_proteins.faa"

LOG="$OUTDIR/combined_faa_log.txt"

MISSING="$OUTDIR/missing_faa.txt"



> "$COMBINED"

> "$LOG"

> "$MISSING"



echo "Combining FAA files from:"

echo "$PROKKA_DIR"

echo



for folder in "$PROKKA_DIR"/TG*_accessory*; do



    if [ -d "$folder" ]; then



        name=$(basename "$folder")



        faa=$(find "$folder" -maxdepth 1 -type f -name "*.faa" | head -n 1)



        if [ -f "$faa" ]; then

            echo "Adding $name"

            echo "$name	$faa" >> "$LOG"



            awk -v prefix="$name" '

            /^>/ {

                sub(/^>/, ">" prefix "|")

                print

                next

            }

            {print}

            ' "$faa" >> "$COMBINED"



        else

            echo "No FAA found in $name"

            echo "$name" >> "$MISSING"

        fi



    fi



done



echo

echo "Done."

echo "Combined file:"

echo "$COMBINED"

echo

echo "Log file:"

echo "$LOG"

echo

echo "Missing FAA file list:"

echo "$MISSING"
