#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                       # 2G for Data needed   
#SBATCH --partition=short              # Specify your cluster partition
#SBATCH --time=1:00:00  
#SBATCH --output=/t3home/gcelotto/slurm/output/miniToNano_70.out  # Output file for stdout
#SBATCH --error=/t3home/gcelotto/slurm/output/miniToNano_70.out    # Output file for stderr
#SBATCH --dependency=singleton
# ---- Parse Input Arguments ----
inputFiles="$1"
outNumber="$2"
massHypo="$3"

# ---- Set path to CMSSW release ----
CMSSW_DIR="/work/gcelotto/CMSSW_12_4_8"

# ---- Run inside container ----
cmssw-el7 --bind /scratch --bind /work/gcelotto/CMSSW_12_4_8 --bind /pnfs/psi.ch/cms/trivcat/store/user/gcelotto << EOF
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd $CMSSW_DIR/src
eval \$(scramv1 runtime -sh)

echo "Mass is $massHypo"
echo "Number is $outNumber"
echo "Input files: $inputFiles"

cmsRun PhysicsTools/BParkingNano/test/run_nano_cfg.py \
    inputFiles=$inputFiles \
    outNumber=$outNumber \
    maxEvents=-1 \
    massHypo=$massHypo \
    outputName="GluGluSpin0_M"

file_path="/scratch/GluGluSpin0_M${massHypo}_Run2_mc_124X_${outNumber}.root"

if [ -f "\$file_path" ]; then
    echo "Transferring \$file_path"
    xrdcp -f -N "\$file_path" root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/GluGluSpin0_private/M${massHypo}/
else
    echo "ERROR: File not found: \$file_path"
fi
EOF

#CMSSW_DIR="/work/gcelotto/CMSSW_12_4_8"
#cmssw-el7 --bind $CMSSW_DIR << EOF
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#
#cd $CMSSW_DIR/src
#eval \$(scramv1 runtime -sh)

#inputFiles="file:/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/pMINLO/GluGluSpin0ToBBbar_W_1p0_M_50_MuEnriched_TuneCP5_13TeV_pythia8_cfi/RunIISummer20UL18MiniAOD_PrivateMC/250625_121624/0000/mini_1.root"\
#outNumber=1 \
#maxEvents=-1 \
#massHypo=50 \
#outputName="GluGluSpin0_M"

