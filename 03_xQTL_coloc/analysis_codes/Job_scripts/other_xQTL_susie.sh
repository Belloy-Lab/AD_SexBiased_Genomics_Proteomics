#!/bin/bash
# Path to the data file
DATA_FILE="../reference_files/other_QTL_reference.csv"

# Read the data file and submit jobs
tail -n +2 "$DATA_FILE" | while IFS=',' read -r study tissue filepath qtl_type || [ -n "$study" ]; do
    # Trim leading/trailing spaces and remove newline or carriage return characters from each variable
    study=$(echo "$study" | tr -d '\r' | xargs)
    tissue=$(echo "$tissue" | tr -d '\r' | xargs)
    filepath=$(echo "$filepath" | tr -d '\r' | xargs)
    qtl_type=$(echo "$qtl_type" | tr -d '\r' | xargs)


    # Print the values for verification
    echo "$study" "$tissue" "$filepath" "$qtl_type" 

    # Construct unique job name and output filenames
    job_name="${study}_${tissue}_other_${qtl_type}_susie"
    output_file="logs/${job_name}.out"
    error_file="logs/${job_name}.err"

    # Submit job
    bsub \
    -J "$job_name" \
    -Ne \
    -oo "$output_file" \
    -eo "$error_file" \
    -q subscription \
    -R "rusage[mem=150GB] span[hosts=1]" \
    -a 'docker(satheshsiva27/multiomics-toolkit:0.1)' \
    Rscript ../other_xQTL_master_susie.R \
		"$study" "$tissue" "$filepath" "$qtl_type" 
done