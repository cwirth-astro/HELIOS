#!/usr/bin/env bash

# Directory containing the archives
ARCHIVE_DIR="./opacity"   # change this if needed

# Loop through all tar and tar.gz files in the archive directory
for file in "$ARCHIVE_DIR"/*.tar*; do
    # skip if nothing matches
    [[ -e "$file" ]] || continue

    # Just the filename
    fname=$(basename "$file")

    # Everything before the first underscore
    prefix="${fname%%_*}"

    # Full path to the extraction directory
    target_dir="$ARCHIVE_DIR/$prefix"

    echo "Creating directory: $target_dir"
    mkdir -p "$target_dir"

    echo "Extracting $fname into $target_dir/"

    # Determine gzip vs plain tar
    if file "$file" | grep -q "gzip compressed"; then
        tar -xzvf "$file" -C "$target_dir"
    else
        tar -xvf "$file" -C "$target_dir"
    fi

done