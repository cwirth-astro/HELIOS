#!/usr/bin/env bash

# Directory containing the archives
ARCHIVE_DIR="./opacity"   # change if needed

for file in "$ARCHIVE_DIR"/*.tar*; do
    [[ -e "$file" ]] || continue

    fname=$(basename "$file")

    # Extract the molecule portion before the first underscore
    molecule_full="${fname%%_*}"

    # Remove isotopes:
    #   - strip leading digits before element symbols
    #   - keep element letters (1–2 letters) + count digits
    # Then remove dashes
    #
    # Example: 31P-1H3 → P-H3 → PH3
    molecule_clean=$(echo "$molecule_full" \
        | sed -E 's/[0-9]+([A-Za-z]{1,2}[0-9]*)/\1/g' \
        | sed 's/-//g')

    target_dir="$ARCHIVE_DIR/$molecule_clean"

    echo "Creating directory: $target_dir"
    mkdir -p "$target_dir"

    echo "Extracting $fname into $target_dir/"

    # Detect gzip vs plain tar
    if file "$file" | grep -qi "gzip compressed"; then
        tar -xzvf "$file" -C "$target_dir"
    else
        tar -xvf "$file" -C "$target_dir"
    fi

done
