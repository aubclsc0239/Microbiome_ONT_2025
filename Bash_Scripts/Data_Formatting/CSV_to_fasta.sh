#!/bin/bash

# this script was written to convert a formatted refSeq in CSV to a Fasta file for next step (alignment)
# this script take a refSeq (CSV) and output a fasta file which can be used for alignment

# Specify the input CSV file and output FASTA file
input_csv="fasta.csv"  # Replace with the path to your CSV file
output_fasta="fasta_file.fasta"  # Replace with the desired path for the FASTA file

# Check if input file exists
if [[ ! -f "$input_csv" ]]; then
    echo "Input file does not exist: $input_csv"
    exit 1
fi

# Clear the output FASTA file if it already exists
> "$output_fasta"

# Process the CSV file line by line
while IFS=',' read -r name sequence || [[ -n "$name" ]]; do
    # Write the header to the FASTA file
    echo ">$name" >> "$output_fasta"
    
    # Break the sequence into lines of 80 characters and append to the FASTA file
    echo "$sequence" | fold -w 80 >> "$output_fasta"
done < "$input_csv"

echo "FASTA file created successfully: $output_fasta"
