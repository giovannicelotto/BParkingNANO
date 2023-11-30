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
options.register('reportEvery', 10,
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

options.setDefault('maxEvents', 20000)
options.setDefault('tag', '124X')
options.parseArguments()
print(options)

globaltag = None
if   options.lhcRun == 3: globaltag = '124X_mcRun3_2022_realistic_v11' if options.isMC else '124X_dataRun3_Prompt_v4'
elif options.lhcRun == 2: globaltag = '102X_upgrade2018_realistic_v15' if options.isMC else '102X_dataRun2_v11'
if options._beenSet['globalTag']: globaltag = options.globalTag

ext1 = {2:'Run2', 3:'Run3'}
ext2 = {False:'data', True:'mc'}
outputFileNANO = cms.untracked.string('_'.join(['Hbb_noTrig',
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
'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/01799CFD-BA2E-844B-87CB-273185CF1A4A.root',
'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/90000/CFBFCA29-1649-2345-970E-731824064446.root',
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/2EA135EF-43C9-624B-A50B-B72B2E2D5393.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/A386A139-01EC-3C4C-AB48-F2FC8814FDF2.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/23EBA0E0-6729-7C41-9888-F8E07B20C677.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/5A00A84C-ED50-B94B-A25D-EE2A6EAD3D4C.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/7010DF71-3814-3047-8786-3A3A894D0E02.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/F338472B-D37A-8040-B482-111D21820AD7.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/DE59C01E-2184-334D-A77A-6D5B18697B02.root', 
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/F7C64AA3-BA68-0F4D-B845-D59C363DF441.root', 
#'/store/mc/RunIISummer20UL18MiniAODv2/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/270000/36E1C10A-F287-A74E-B874-6C359FACD755.root'
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/2C416EBB-5DEE-E811-B2B4-D48564592B02.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/14A607B2-5DEE-E811-B6B1-509A4C74D08F.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/3A86E8AC-5DEE-E811-8B2F-002590907826.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/988407DC-5DEE-E811-A28C-A4BF0101DB93.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/C074AF04-2BED-E811-B7D8-0CC47AFCC6A6.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/1A79E2C6-5DEE-E811-8692-0025905C3D6C.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/F4947FA5-25EC-E811-BE29-90B11C0DCA4B.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/D2946C2A-35EC-E811-94D2-002590E3A224.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/5A82EFA7-3CEC-E811-876D-0CC47AD98B8E.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/E0B52B60-48EC-E811-8B72-0CC47AD98D6E.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/90F5498C-59EC-E811-B2D1-1C6A7A26BCDB.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/54CA8ECF-4BED-E811-AD93-002590E39D52.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/1876BAD5-5DEE-E811-BE0C-0CC47AD9908C.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/28FEBABF-3CEC-E811-90CF-002590D9D8AE.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/E8F803A0-80EC-E811-98D3-0CC47AB0B704.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/7C54C8B8-B7ED-E811-804F-002590FD5A72.root'
] if options.isMC else [
        '/store/data/Run2018A/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2430000/004BEEAD-CCCD-4A4F-9217-91A5A28EA0C8.root'
        ]
    elif options.lhcRun == 3:
        options.inputFiles = [
            'root://cms-xrd-global.cern.ch//store/user/jodedra/BuTOjpsiKEE20221103FIFTYMminiaod/BuTOjpsiKEE20221103FIFTYM/SUMMER22_MINIAOD/221106_001759/0000/step1_inMINIAODSIM_1.root',
        ] if options.isMC else [
            'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleElectronLowMass0/MINIAOD/PromptReco-v1/000/356/170/00000/45c0f2ed-eb5b-4292-abc8-3117424d9432.root'
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
      'drop *',
      "keep nanoaodFlatTable_*Table*_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.BParkingNano.nanoBPark_cff import *
if options.lhcRun == 2:
    process = nanoAOD_customizeMuonTriggerBPark(process)
    process = nanoAOD_customizeElectronFilteredBPark(process)
    process = nanoAOD_customizeTrackFilteredBPark(process)
    #process = nanoAOD_customizeBToKLL(process)
    #process = nanoAOD_customizeBToKstarEE(process)
    #process = nanoAOD_customizeBToKstarMuMu(process)
elif options.lhcRun == 3:
    from PhysicsTools.BParkingNano.electronsTrigger_cff import *
    process = nanoAOD_customizeDiEle(process)
    process = nanoAOD_customizeElectronFilteredBPark(process)
    process = nanoAOD_customizeTriggerBitsBPark(process)
    process = nanoAOD_customizeTrackFilteredBPark(process)
    process = nanoAOD_customizeBToKLL(process)

# Path and EndPath definitions
if options.lhcRun == 2:
    process.nanoAOD_Jets_step = cms.Path(process.nanoSequence)
    #process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.nanoBKMuMuSequence + CountBToKmumu )
    #process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.nanoBKeeSequence   + CountBToKee   )
    #process.nanoAOD_KstarMuMu_step = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.KstarToKPiSequence + process.nanoBKstarMuMuSequence + CountBToKstarMuMu )
    #process.nanoAOD_KstarEE_step  = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.KstarToKPiSequence + process.nanoBKstarEESequence + CountBToKstarEE  )
elif options.lhcRun == 3:
    process.nanoAOD_DiEle_step = cms.Path(process.nanoSequence
                                          #+process.nanoDiEleSequence
                                          #+process.nanoTracksSequence
                                          #+process.nanoBKeeSequence
                                          #+CountBToKee
                                          )

# customisation of the process.
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
if options.lhcRun == 3:

    process.schedule = cms.Schedule(process.nanoAOD_DiEle_step,
                                    process.endjob_step,
                                    process.NANOAODoutput_step)
    if options.wantFullRECO:
        process.schedule = cms.Schedule(process.nanoAOD_DiEle_step,
                                        process.endjob_step,
                                        process.FEVTDEBUGHLToutput_step,
                                        process.NANOAODoutput_step)
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)
    process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('nanoAOD_DiEle_step')
    )

elif options.lhcRun == 2:

    process.schedule = cms.Schedule(
        process.nanoAOD_Jets_step,
        #process.nanoAOD_Kee_step,
        #process.nanoAOD_KstarMuMu_step,
        #process.nanoAOD_KstarEE_step,
        process.endjob_step,
        process.NANOAODoutput_step
    )
    if options.wantFullRECO:
        process.schedule = cms.Schedule(
            process.nanoAOD_KMuMu_step,
            #process.nanoAOD_Kee_step,
            #process.nanoAOD_KstarMuMu_step,
            #process.nanoAOD_KstarEE_step,
            process.endjob_step,
            process.FEVTDEBUGHLToutput_step,
            process.NANOAODoutput_step
        )
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)
    process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'nanoAOD_Jets_step',
            #'nanoAOD_Kee_step',
            #'nanoAOD_KstarMuMu_step',
            #'nanoAOD_KstarEE_step',
        )
    )

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
