from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)
options.register('lhcRun', 2,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "LHC Run 2 or 3 (default)"
)
options.register('outNumber', -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "number of the outFile"
)
options.register('outputName', "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "output Name"
)
options.register('massHypo', -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "mass of the Spin 0 particle"
)
options.setDefault('maxEvents', 10000)
options.setDefault('tag', '124X')
options.parseArguments()
print(options)

globaltag = None
if   options.lhcRun == 3: globaltag = '124X_mcRun3_2022_realistic_v11' if options.isMC else '124X_dataRun3_Prompt_v4'
elif options.lhcRun == 2: globaltag = '106X_upgrade2018_realistic_v16' if options.isMC else '106X_dataRun2_v35'
if options._beenSet['globalTag']: globaltag = options.globalTag
ext1 = {2:'Run2', 3:'Run3'}
ext2 = {False:'data', True:'mc'}


processName = "DataFiltered" if options.outputName=="" else options.outputName

print(options.outputName=="")

if options.outNumber!=-1:
    print("You are here")
    print("Your outNumber is ", str(options.outNumber))

    #print('/scratch/'+'_'.join(['ZJetsToQQ_HT-100to200',
    #                                            ext1[options.lhcRun],
    #                                            ext2[options.isMC],
    #                                            options.tag,
    #                                            options.outNumber])+'.root')
    if str(options.massHypo)!="-1":
        outputFileNANO = cms.untracked.string('/scratch/'+'_'.join([processName+str(options.massHypo),
                                                ext1[options.lhcRun],
                                                ext2[options.isMC],
                                                options.tag,
                                                str(options.outNumber)])+'.root')
    else:
        #outputFileNANO = cms.untracked.string('/scratch/'+'_'.join([processName,
        #                                        ext1[options.lhcRun],
        #                                        ext2[options.isMC],
        #                                        options.tag,
        #                                        options.outNumber])+'.root')
        outputFileNANO = cms.untracked.string('_'.join([processName,
                                                ext1[options.lhcRun],
                                                ext2[options.isMC],
                                                options.tag,
                                                str(options.outNumber)])+'.root')
        print(outputFileNANO)

else:
    outputFileNANO = cms.untracked.string('_'.join([processName,
                                                ext1[options.lhcRun],
                                                ext2[options.isMC],
                                                options.tag])+'.root')

outputFileFEVT = cms.untracked.string('_'.join(['BParkingFullEvt',
                                                ext1[options.lhcRun],
                                                ext2[options.isMC],
                                                options.tag])+'.root')
if not options.inputFiles:
    if options.lhcRun == 2:
        options.inputFiles = [
'/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/00F7273C-6F52-7D4E-8175-E86320D6068A.root'
#'/store/mc/RunIISummer20UL18MiniAODv2/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2520000/002F46B9-4287-194A-BA9B-469CFB34D146.root'
#'/store/mc/RunIISummer20UL18MiniAODv2/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/270000/08E9309F-52F4-C74A-B740-06BD4A21A61E.root',
#'file:/t3home/gcelotto/moreMC/mini.root'
] if options.isMC else [
'/store/data/Run2018A/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2430000/004BEEAD-CCCD-4A4F-9217-91A5A28EA0C8.root'
]
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

# Process
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.BParkingNano.modifiers_cff import *
process = None
if   options.lhcRun == 3: process = cms.Process('BParkNANO',eras.Run3,BToKEE_DiEle)
elif options.lhcRun == 2: process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring(
        'keep *',
        'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
    'drop *',  # Drop everything by default
    "keep nanoaodFlatTable_*Table*_*_*",  # Keep event-level FlatTables
    "keep nanoaodUniqueString_nanoMetadata_*_*",  # Keep basic metadata
    "keep nanoaodMergeableCounterTable_*_*_*",
    #"keep TTree_Runs_*_*"

)

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.BParkingNano.nanoBPark_cff import *

# Definition of the NanoSequence

process = nanoAOD_customizeMuonTriggerBPark(process)
process = nanoAOD_customizeElectronFilteredBPark(process)

#process = nanoAOD_customizeTrackFilteredBPark(process)     #Tracks are removed




# Path and EndPath definitions
process.nanoAOD_Jets_step = cms.Path(process.nanoSequence)

# customisation of the process.
# nanoAOD_customizeMC modifies the Path. It needs to be defined after nanoAOD_Jets_step
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)


process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)





process.schedule = cms.Schedule(
    process.nanoAOD_Jets_step,
    process.endjob_step,
    process.NANOAODoutput_step
    )
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Save only events that pass skimming
process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring(
        'nanoAOD_Jets_step',
    )
)

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))        #Disable ROOT MultiThreading (can conflict with Crab)
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)                              # Used to trick CRAB into recognizing the output file even if it has a dynamic name or is created programmatically

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

print("\n========= Schedule =========")
print(process.schedule)
print("============================\n")
for name, path in process.paths.iteritems():
    print("Name : %s\nPath : %s"%(name, path))
