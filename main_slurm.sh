#!/bin/bash

#SBATCH -p pool1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o slurm_out/%x_%a.out
#SBATCH -e slurm_out/%x_%a.err


# Check if the number of arguments is correct
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 project_name"
    exit 1
fi

# Access command-line arguments
project_name="$1"
idx=$SLURM_ARRAY_TASK_ID
data_dir="data/$project_name"

i=0

find $data_dir/*.maf -type f | while IFS= read -r file; do
    filesize=$(stat -f%z "$file")
    fname=$(basename "$file")
    barcode="${fname%.*}"
    if [ "$i" -eq "$idx" ] && [ "$filesize" -gt 100000 ] && [ "$filesize" -lt 2000000 ]; then
        echo "Processing $barcode"
        bash infer_sample.sh $project_name $barcode
    fi
    i=$((i+1))
done