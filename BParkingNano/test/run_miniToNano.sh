#!/bin/bash

# Check if the argument mHypo is provided
if [ -z "$1" ]; then
  echo "Please provide mHypo as an argument."
  exit 1
fi
mHypo="$1"
# Assign the DIRECTORY based on mHypo
case "$mHypo" in
    50)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/ratramon/GluGluSpin0ToBBbar_W_1p0_M_50/RunIISummer20UL18_MINI_fullStat/240815_194329/0000"
        ;;
    70)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/ratramon/GluGluSpin0ToBBbar_W_1p0_M_70/RunIISummer20UL18_MINI/240815_194311/0000"
        ;;
    100)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/ratramon/GluGluSpin0ToBBbar_W_1p0_M_100/RunIISummer20UL18_MINI/240815_194320/0000"
        ;;
    200)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/ratramon/GluGluSpin0ToBBbar_W_1p0_M_200/RunIISummer20UL18_MINI/240815_194329/0000"
        ;;
    300)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/ratramon/GluGluSpin0ToBBbar_W_1p0_M_300/RunIISummer20UL18_MINI/240815_194338/0000"
        ;;
    *)
        echo "Invalid mHypo value. Please choose from 50, 70, 100, 200, or 300."
        exit 1
        ;;
esac

# Print the selected directory
echo "Selected DIRECTORY: $DIRECTORY"

#DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/miniAODZJets/ZJetsToQQ_HT100to200/RunIISummer20UL18/240523_144055/0000"
# Define the CMS base directory
cmsbase="/t3home/gcelotto/CMSSW_12_4_8/src"

# Define the path to the run_nano.py script
SCRIPT_PATH="$cmsbase/PhysicsTools/BParkingNano/test/run_nano_cfg.py"
cd $cmsbase || { echo "Failed to change directory to $cmsbase"; exit 1; }
cmsenv
echo "CMS environment base: $CMSSW_BASE"

# Verify the run_nano.py script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script $SCRIPT_PATH does not exist."
    exit 1
fi


#count=0
for FILE in "$DIRECTORY"/*; do
    if [ -f "$FILE" ]; then
        
        #echo "$FILE"
        BASENAME=$(basename "$FILE")
        NUMBER=$(echo "$BASENAME" | awk -F'MINI_|\\.root' '{print $2}')
        echo "Extracted number: $NUMBER"
        
        if [ -f "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_M"$mHypo"/GluGluSpin0_M"$mHypo"_"$NUMBER".root" ]; then
            echo "File exists: $FILE"
        else
            #if [ $count -lt 10 ]; then
            #    echo "Processing file: $FILE"
            #    # Add your processing logic here
            #    ((count++))
            #else
            #    break
            #fi
            MOD_RESULT=$((NUMBER % 25))
            sbatch -J "ggH_M"$mHypo"_$MOD_RESULT" /t3home/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/nano_job.sh "root://t3dcachedb.psi.ch:1094//$FILE" "$NUMBER" "$mHypo"
        fi
        
        
        #cmsRun /t3home/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/run_nano_cfg.py inputFiles= outNumber=$NUMBER

    fi
done


