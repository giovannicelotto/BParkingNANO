#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=resubmit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G                       # 2G for Data needed   
#SBATCH --partition=short              # Specify your cluster partition
#SBATCH --time=01:00:00  
#SBATCH --output=/t3home/gcelotto/slurm/output/resubmit.out  # Output file for stdout
#SBATCH --error=/t3home/gcelotto/slurm/output/resubmit.out    # Output file for stderr
#SBATCH --dependency=singleton

# Base directory to start the search from
BASE_DIR="$1"

# Check if a directory is provided
if [ -z "$BASE_DIR" ]; then
  echo "Not an argument provided for base directory"
  echo "Using ./ as argument"
  BASE_DIR="./"
fi

# Find all directories recursively and execute `crab resubmit` on them
find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r folder; do
  echo "Executing: crab resubmit $folder"
  crab resubmit "$folder"
done
