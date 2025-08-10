from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
#from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.BParkingNano.trgbits_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.boostedTaus_cff import *
#from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jets_cff import *
from PhysicsTools.NanoAOD.dijet_cff import dijetBuilder, dijetTable, dijetSequence, dijetPtFilter

#from PhysicsTools.NanoAOD.jetsAK4_Puppi_cff import *
#from PhysicsTools.NanoAOD.jetMC_cff import *
from PhysicsTools.NanoAOD.muons_cff import *

##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.tracksBPark_cff import *

## B collections
#from PhysicsTools.BParkingNano.BToKLL_cff import *
#from PhysicsTools.BParkingNano.BToKstarLL_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)
linkedObjectsNew = cms.EDProducer("PATObjectCrossLinker",
   jets=cms.InputTag("finalJets"),
   muons=cms.InputTag("muonTrgSelector:SelectedMuons"),
   electrons=cms.InputTag("slimmedElectrons"),
   lowPtElectrons=cms.InputTag(""),
   taus=cms.InputTag("slimmedTaus"),
   boostedTaus=cms.InputTag(""),
   photons=cms.InputTag("slimmedPhotons"),
   vertices=cms.InputTag("slimmedSecondaryVertices")
)

#linkedObjectsCR = cms.EDProducer("PATObjectCrossLinker",
#   jets=cms.InputTag("finalJets"),
#   muons=cms.InputTag("finalMuons"),
#   electrons=cms.InputTag("slimmedElectrons"),
#   lowPtElectrons=cms.InputTag(""),
#   taus=cms.InputTag("slimmedTaus"),
#   boostedTaus=cms.InputTag(""),
#   photons=cms.InputTag("slimmedPhotons"),
#   vertices=cms.InputTag("slimmedSecondaryVertices")
#)

nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequenceCommon = cms.Sequence(nanoMetadata +
                                 muonBParkSequence + 
                                 jetSequence +                                          # dont know how to handle the muonsubptraw missing
                                 #cms.Sequence(jetTask) +
                                    vertexTables+
                                 linkedObjectsNew +
                                 jetTables+
                                 dijetSequence
                            #cms.Sequence(vertexTask) + 
                             )
nanoSequence = cms.Sequence(nanoSequenceCommon + nanoSequenceOnlyFullSim)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genWeightsTableSequence+
                                genParticleBParkSequence + nanoSequenceCommon + jetMC +
                              globalTablesMC  + genParticleBParkTables + lheInfoTable
                              )




def nanoAOD_customizeMuonTriggerBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence #+ muonBParkSequence 
                                        + muonBParkTables
                                        )# +
                                        #slimmedMuonsUpdated+isoForMu+ptRatioRelForMu+slimmedMuonsWithUserData+ finalMuons + linkedObjectsCR+muonTable)#+ muonTriggerMatchedTables)   ###comment in this extra table in case you want to create the TriggerMuon collection again.
    return process

def nanoAOD_customizeTrackFilteredBPark(process):
    process.nanoSequence = cms.Sequence(process.nanoSequence+ tracksBParkSequence + tracksBParkTables)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoSequence     = cms.Sequence(process.nanoSequence+ electronsBParkSequence + electronBParkTables)
    #process.nanoBKstarEESequence = cms.Sequence( electronsBParkSequence + electronBParkTables)
    return process

def nanoAOD_customizeTriggerBitsBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeBToKLL(process):
    #process.nanoBKeeSequence   = cms.Sequence( process.nanoBKeeSequence + BToKEESequence    + BToKeeTable   )
    #process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
    return process

#three possibilities for K*LL
def nanoAOD_customizeBToKstarLL(process):
    #process.nanoBKstarLLSequence   = cms.Sequence( KstarToKPiSequence + BToKstarLLSequence + KstarToKPiTable + BToKstarLLTables )
    return process

def nanoAOD_customizeBToKstarEE(process):
    #process.nanoBKstarEESequence   = cms.Sequence( process.nanoBKstarEESequence + BToKstarEESequence + BToKstarEETable + KstarToKPiTable )
    return process

def nanoAOD_customizeBToKstarMuMu(process):
    #process.nanoBKstarMuMuSequence = cms.Sequence( BToKstarMuMuSequence + BToKstarMuMuTable + KstarToKPiTable )
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    print("Called nanoADO_customizeMC")
    for name, path in process.paths.iteritems():
        print("Name : %s\nPath : %s"%(name, path))
        # replace all the non-match embedded inputs with the matched ones (are the same muons but with matching to genMuons)
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        #massSearchReplaceAnyInputTag(path, 'electronTrgSelector:SelectedElectrons', 'selectedElectronsMCMatchEmbedded') # Is this needed if the trigger is emulated ???
        massSearchReplaceAnyInputTag(path, 'electronsForAnalysis:SelectedElectrons', 'selectedElectronsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # modify the path to include mc-specific info
        #print("This is nanoSequenceMC ", nanoSequenceMC)
        path.remove(nanoSequenceCommon)
        path.insert(0, nanoSequenceMC)
        # muonBParkMC = cms.Sequence(muonBParkSequence + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable)
        # muonBparkMC includes also bpark sequence
        path.replace(process.muonBParkSequence, process.muonBParkMC)
        path.replace(process.electronsBParkSequence, process.electronBParkMC)
        path.replace(process.tracksBParkSequence, process.tracksBParkMC)

        print("\n\nNew Path : %s"%path)
