#!/bin/bash

# Input files
species_file="species_only.txt"
sequence_file="sequences.csv"
output_file="output_species_with_sequences.csv"

# Initialize the output file with headers
echo "Species,Sequence" > "$output_file"

# Read the sequences from the sequences.csv file into an associative array
declare -A sequences
while IFS=, read -r species sequence
do
    # Skip the header line
    if [[ "$species" != "Bacteria Name" ]]; then
        # Remove extra spaces or line breaks from both species and sequence
        species=$(echo "$species" | tr -d '\r\n' | sed 's/^[ \t]*//;s/[ \t]*$//')
        sequence=$(echo "$sequence" | tr -d '\r\n' | sed 's/^[ \t]*//;s/[ \t]*$//')
        sequences["$species"]="$sequence"
    fi
done < "$sequence_file"

# Process the species file
while IFS= read -r species
do
    # Remove extra spaces or line breaks from species name in species.txt
    species=$(echo "$species" | tr -d '\r\n' | sed 's/^[ \t]*//;s/[ \t]*$//')

    # Check if the species exists in the sequences array
    if [[ -n "${sequences["$species"]}" ]]; then
        # If found, prepend the sequence to the species name
        echo "${species},${sequences["$species"]}" >> "$output_file"
    else
        # If not found, just write the species name and leave it unchanged
        echo "${species},${species}" >> "$output_file"
    fi
done < "$species_file"

echo "Process complete. The output file is saved as $output_file."

