#!/bin/bash
DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/miniAODZJets/ZJetsToQQ_HT100to200/RunIISummer20UL18/allFiles"

cmsbase="/work/gcelotto/CMSSW_12_4_8/src"

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
        NUMBER=$(echo "$BASENAME" | awk -F'_|\\.root' '{print $2}')
        echo "Extracted number: $NUMBER"
        echo "File $FILE"
        
        if [ -f "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ZJetsToQQ_noTrig2025Jun03/ZJetsToQQ_HT-100to200/ZJetsToQQ_HT100to200_Run2_mc_124X_"$NUMBER".root" ]; then
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
            echo "OutNumber : $NUMBER"
            sbatch -J "ZJets_$MOD_RESULT" /work/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/nano_job_ZJets100To200.sh "root://t3dcachedb.psi.ch:1094//$FILE" "$NUMBER" 
        fi
        
        
        #cmsRun /t3home/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/run_nano_cfg.py inputFiles= outNumber=$NUMBER

    fi
done


