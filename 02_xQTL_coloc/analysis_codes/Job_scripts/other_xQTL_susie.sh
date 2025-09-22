#!/bin/bash
# Path to the data file
DATA_FILE="/storage2/fs1/belloy2/Active/04_Code/$USER/xQTL/reference_files/other_QTL_reference.csv"

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
    output_file="/storage2/fs1/belloy2/Active/05_Projects/$USER/xQTL/logs/${job_name}.out"
    error_file="/storage2/fs1/belloy2/Active/05_Projects/$USER/xQTL/logs/${job_name}.err"

    # Submit job
    bsub \
    -g /satheshkumar/compute-belloy \
    -J "$job_name" \
    -Ne \
    -oo "$output_file" \
    -eo "$error_file" \
    -q general \
    -R "rusage[mem=150GB] span[hosts=1] select[hname!='compute1-exec-40.ris.wustl.edu' && hname!='compute1-exec-160.ris.wustl.edu' && hname!='compute1-exec-106.ris.wustl.edu' && hname!='compute1-exec-175.ris.wustl.edu']" \
    -sp 95 \
    -G compute-belloy-t1 \
    -sla compute-belloy-t1 \
    -a 'docker(happygoat1946/sctakin:4.6.2)' \
    Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/xQTL/other_xQTL_master_susie.R \
		"$study" "$tissue" "$filepath" "$qtl_type" 
done