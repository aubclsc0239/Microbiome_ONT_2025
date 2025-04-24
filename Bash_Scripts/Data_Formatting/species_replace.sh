#!/bin/bash

# This script is written to create a ref_seq file after Epi2me run by search for sequences of found species from NCBI database
# This script takes the NCBI 18S/16S database as a reference file (sequences.csv)
# Note that the database in not uploaded becuase of the size. You can download the database follwoing Epi2me run. 
# The species list is made from Epi2me output
# This script output a fasta.csv file with species and sequence to be used for alignment, then phylogenetic tree

# Input files
species_file="species_only.txt"
sequence_file="sequences.csv"
output_file="fasta.csv"

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

