#!/bin/bash

# Specify the input folder containing fastq.gz files
input_folder="../../221107_Diabetes/dia_hla2/"

# Specify the input file containing patterns
patterns_file="patterns.txt"

# Specify the number of CPU cores to use for parallel processing
num_cores=18  # Change this to the desired number of cores

# Check if the patterns file exists
if [ ! -e "${patterns_file}" ]; then
  echo "Patterns file not found: ${patterns_file}"
  exit 1
fi

# Combine all patterns from the input file into a single pattern
patterns=$(cat "${patterns_file}" | sed 's/^[^ ]* //; s/|$//')

# Function to process a single file
process_file() {
  input_file="$1"
  output_file="${input_file%.fastq.gz}.txt"
  
  zcat "${input_file}" | awk 'NR%4==2' | grep -E "${patterns}" > "${output_file}"
  echo "Processed ${input_file} and saved the result in ${output_file}"
}

export -f process_file

# List all .fastq.gz files in the input folder
files=("${input_folder}"/*.fastq.gz)

# Use parallel to process the files in parallel
parallel -j "${num_cores}" process_file ::: "${files[@]}"

