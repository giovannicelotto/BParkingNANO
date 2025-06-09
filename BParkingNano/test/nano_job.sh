#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                       # 2G for Data needed   
#SBATCH --partition=short              # Specify your cluster partition
#SBATCH --time=1:00:00  
#SBATCH --output=/t3home/gcelotto/slurm/output/miniToNano.out  # Output file for stdout
#SBATCH --error=/t3home/gcelotto/slurm/output/miniToNano.out    # Output file for stderr
#SBATCH --dependency=singleton

# Load the CMS environment
cd /work/gcelotto/CMSSW_12_4_8/src
cmsenv

# Run the CMS analysis script
cmsRun /work/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/run_nano_cfg.py inputFiles="$1" outNumber="$2" maxEvents=-1 massHypo="$3"

#file_path="/scratch/ZJetsToQQ_100to200_${2}.root"
file_path="/scratch/GluGluSpin0_M${3}_Run2_mc_124X_${2}.root"

# Check if the file size is larger than 2 MB (2*1024*1024 bytes)
#if [ $(stat -c%s "$file_path") -gt 1024 ]; then
    # Perform the file transfer using xrdcp
xrdcp -f -N "$file_path" root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/MCfiducial_corrections2025Mar10/GluGluSpin0_M${3}/
#else
#    echo "File size is less than 2 MB, skipping transfer."
#fi
