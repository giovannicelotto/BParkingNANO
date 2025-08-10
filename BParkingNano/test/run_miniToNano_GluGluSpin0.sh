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
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_50_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250626_110446/0000"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M50"
        ;;
    70)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_70_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250626_110510/0000"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M70"
        ;;
    100)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_100_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250626_110534/0000"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M100"
        ;;
    125)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD_PrivateMC/combined"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M125"
        ;;
    200)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_200_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250626_110601/0000"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M200"
        ;;
    300)
        DIRECTORY="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_300_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250626_110626/0000"
        DIRECTORY_NANO="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M300"
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
cmsbase="/work/gcelotto/CMSSW_12_4_8/src"

# Define the path to the run_nano.py script
SCRIPT_PATH="$cmsbase/PhysicsTools/BParkingNano/test/run_nano_cfg.py"
cd $cmsbase || { echo "Failed to change directory to $cmsbase"; exit 1; }
#cmsenv
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
        #NUMBER=$(echo "$BASENAME" | awk -F"GluGluSpin0_M${mHypo}_mini_|\\.root" '{print $2}')
        NUMBER=$(echo "$BASENAME" | awk -F"MINLO_Private_|\\.root" '{print $2}')
        echo "Extracted number: $NUMBER"
        echo "File $FILE"
        
        if [ -f $DIRECTORY_NANO"/GluGluSpin0_M"$mHypo"_Run2_mc_124X_"$NUMBER".root" ]; then
            echo "File exists: $FILE"
        else
            #if [ $count -lt 10 ]; then
            #    echo "Processing file: $FILE"
            #    # Add your processing logic here
            #    ((count++))
            #else
            #    break
            #fi
            MOD_RESULT=$((NUMBER % 100))
            sbatch -J "ggH_M"$mHypo"_$MOD_RESULT" /work/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/nano_job_new.sh "root://t3dcachedb.psi.ch:1094//$FILE" "$NUMBER" "$mHypo"
            #exit 1

        fi
        
        
        #cmsRun /t3home/gcelotto/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/run_nano_cfg.py inputFiles= outNumber=$NUMBER

    fi
done