# linkedObjects -> linkedOBjectsNew
# Comment the following lines:
#   nThreads = cms.uint32(1),
#   singleThreadPool = cms.string("no_threads"),
import os 
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_nanoAOD_94X2016_cff import run2_nanoAOD_94X2016

from CondCore.CondDB.CondDB_cfi import CondDB
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.NanoAOD.common_cff import Var



from PhysicsTools.NanoAOD.nanoDQM_tools_cff import *

from  PhysicsTools.NanoAOD.common_cff import *
from RecoJets.JetProducers.ak4PFJetsBetaStar_cfi import *


##################### User floats producers, selectors ##########################
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

chsForSATkJets = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string('charge()!=0 && pvAssociationQuality()>=5 && vertexRef().key()==0'))
softActivityJets = ak4PFJets.clone(src = 'chsForSATkJets', doAreaFastjet = False, jetPtMin=1) 
softActivityJets10 = cms.EDFilter("CandPtrSelector", src = cms.InputTag("softActivityJets"), cut = cms.string('pt>10'))
softActivityJets5 = cms.EDFilter("CandPtrSelector", src = cms.InputTag("softActivityJets"), cut = cms.string('pt>5'))
softActivityJets2 = cms.EDFilter("CandPtrSelector", src = cms.InputTag("softActivityJets"), cut = cms.string('pt>2'))

from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
# Note: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
#      (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CMSSW_7_6_4_and_above )
jetCorrFactorsNano = patJetCorrFactors.clone(src='slimmedJets',
    levels = cms.vstring('L1FastJet',
        'L2Relative',
        'L3Absolute',
	'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
#jetCorrFactorsAK8 = patJetCorrFactors.clone(src='slimmedJetsAK8',
#    levels = cms.vstring('L1FastJet',
#        'L2Relative',
#        'L3Absolute',
#	'L2L3Residual'),
#    payload = cms.string('AK8PFPuppi'),
#    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#)
#run2_miniAOD_80XLegacy.toModify(jetCorrFactorsAK8, payload = cms.string('AK8PFchs')) # ak8PFJetsCHS in 2016 80X miniAOD

from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *

updatedJets = updatedPatJets.clone(
	addBTagInfo=False,
	jetSource='slimmedJets',
	jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsNano") ),
)

#updatedJetsAK8 = updatedPatJets.clone(
#	addBTagInfo=False,
#	jetSource='slimmedJetsAK8',
#	jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsAK8") ),
#)
jecSources = [
    "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
    "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta",
    "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
    "RelativeBal", "RelativeSample", "RelativeFSR", "PileUpDataMC",
    "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"
]


looseJetId = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('WINTER16'),
			    quality = cms.string('LOOSE'),
			  ),
                          src = cms.InputTag("updatedJets")
)
tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('WINTER17'),
			    quality = cms.string('TIGHT'),
			  ),
                          src = cms.InputTag("updatedJets")
)
tightJetIdLepVeto = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('WINTER17'),
			    quality = cms.string('TIGHTLEPVETO'),
			  ),
                          src = cms.InputTag("updatedJets")
)
for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
    modifier.toModify( tightJetId.filterParams, version = "WINTER16" )
    modifier.toModify( tightJetIdLepVeto.filterParams, version = "WINTER16" )


#looseJetIdAK8 = cms.EDProducer("PatJetIDValueMapProducer",
#			  filterParams=cms.PSet(
#			    version = cms.string('WINTER16'),
#			    quality = cms.string('LOOSE'),
#			  ),
#                          src = cms.InputTag("updatedJetsAK8")
#)
#tightJetIdAK8 = cms.EDProducer("PatJetIDValueMapProducer",
#			  filterParams=cms.PSet(
#			    version = cms.string('WINTER17'),
#			    quality = cms.string('TIGHT'),
#			  ),
#                          src = cms.InputTag("updatedJetsAK8")
#)
#tightJetIdLepVetoAK8 = cms.EDProducer("PatJetIDValueMapProducer",
#			  filterParams=cms.PSet(
#			    version = cms.string('WINTER17'),
#			    quality = cms.string('TIGHTLEPVETO'),
#			  ),
#                          src = cms.InputTag("updatedJetsAK8")
#)
#for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
#    modifier.toModify( tightJetIdAK8.filterParams, version = "WINTER16" )
#    modifier.toModify( tightJetIdLepVetoAK8.filterParams, version = "WINTER16" )


bJetVars = cms.EDProducer("JetRegressionVarProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("updatedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
    #musrc = cms.InputTag("slimmedMuons"),
    #elesrc = cms.InputTag("slimmedElectrons")
)

jercVars = cms.EDProducer("BetaStarPackedCandidateVarProducer",
    srcJet = cms.InputTag("updatedJets"),    
    srcPF = cms.InputTag("packedPFCandidates"),
    maxDR = cms.double(0.4)
)

#variables for 2018 regression
updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
     src = cms.InputTag("updatedJets"),
     userFloats = cms.PSet(
         leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
         leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
         leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
         leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
         leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
         leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
         leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
         leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
         leptonPt = cms.InputTag("bJetVars:leptonPt"),
         vtxPt = cms.InputTag("bJetVars:vtxPt"),
         vtxMass = cms.InputTag("bJetVars:vtxMass"),
         vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
         vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
         ptD = cms.InputTag("bJetVars:ptD"),
         genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
         qgl = cms.InputTag('qgtagger:qgLikelihood'),
         #jercCHPUF = cms.InputTag("jercVars:chargedHadronPUEnergyFraction"),
         #jercCHF = cms.InputTag("jercVars:chargedHadronCHSEnergyFraction"),
         ),
     userInts = cms.PSet(
        tightId = cms.InputTag("tightJetId"),
        tightIdLepVeto = cms.InputTag("tightJetIdLepVeto"),
        vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
        leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
        
     ),
)
### Evaluate particle-net regression + classification                    
pnetDiscriminatorNames = [];
pnetDiscriminatorLabels = [];

#evaluator from https://gitlab.cern.ch/rgerosa/particlenetstudiesrun2/-/blob/cmssw_13X_new_features/TrainingNtupleMakerAK4/plugins/ParticleNetFeatureEvaluator.cc?ref_type=heads
# included under PhysicsTools/NanoAOD/plugind
pfParticleNetAK4LastJetTagInfos = cms.EDProducer("ParticleNetFeatureEvaluator",
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondary_vertices = cms.InputTag("slimmedSecondaryVertices"),
    #jets = cms.InputTag("slimmedJetsCalibrated" if options.applyJECs else ("slimmedJetsPuppi" if (options.usePuppiJets and not options.reRunPuppi) else ("patJetsPuppi" if (options.usePuppiJets and options.reRunPuppi) else "slimmedJets"))),
    jets = cms.InputTag("updatedJetsWithUserData"),# if options.applyJECs else ("slimmedJetsPuppi" if (options.usePuppiJets and not options.reRunPuppi) else ("patJetsPuppi" if (options.usePuppiJets and options.reRunPuppi) else "slimmedJets"))),
    losttracks = cms.InputTag("lostTracks"),
    jet_radius = cms.double(0.4),
    min_jet_pt = cms.double(10),
    max_jet_eta = cms.double(2.5),
    min_jet_eta = cms.double(0),
    min_pt_for_pfcandidates = cms.double(0.1),
    min_pt_for_track_properties = cms.double(-1),
    min_pt_for_losttrack = cms.double(1.0),
    max_dr_for_losttrack = cms.double(0.2),
    min_pt_for_taus = cms.double(20.),
    max_eta_for_taus = cms.double(2.5),
    dump_feature_tree = cms.bool(False),
    use_puppiP4 = cms.bool(False),
    puppi_weights = cms.InputTag(""),
)
#negative tags for calib against fake - just run for now, not tabled
pfParticleNetAK4LastNegativeJetTagInfos = pfParticleNetAK4LastJetTagInfos.clone(
    flip_ip_sign = cms.bool(True)
)

from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
pfParticleNetAK4LastJetTags = boostedJetONNXJetTagsProducer.clone();
pfParticleNetAK4LastJetTags.src = cms.InputTag("pfParticleNetAK4LastJetTagInfos");
pfParticleNetAK4LastJetTags.flav_names = cms.vstring('probmu','probele','probtaup1h0p','probtaup1h1p','probtaup1h2p','probtaup3h0p','probtaup3h1p','probtaum1h0p','probtaum1h1p','probtaum1h2p','probtaum3h0p','probtaum3h1p','probb','probc','probuds','probg','ptcorr','ptreshigh','ptreslow','ptnu');

#loading model and preproc for CHS jets 
pfParticleNetAK4LastJetTags.preprocess_json = cms.string('PhysicsTools/NanoAOD/data/ParticleNetStudiesRun2/ParticleNetAK4/PNETUL/ClassRegTwoTargets/preprocess.json');
pfParticleNetAK4LastJetTags.model_path = cms.FileInPath('PhysicsTools/NanoAOD/data/ParticleNetStudiesRun2/ParticleNetAK4/PNETUL/ClassRegTwoTargets/particle-net.onnx');
pfParticleNetAK4LastJetTags.debugMode = cms.untracked.bool(False)

pfParticleNetAK4LastNegativeJetTags = pfParticleNetAK4LastJetTags.clone();
pfParticleNetAK4LastNegativeJetTags.src = cms.InputTag("pfParticleNetAK4LastNegativeJetTagInfos");

pnetDiscriminatorNames.extend([
    "pfParticleNetAK4LastJetTags:probmu",
    "pfParticleNetAK4LastJetTags:probele",
    "pfParticleNetAK4LastJetTags:probtaup1h0p",
    "pfParticleNetAK4LastJetTags:probtaup1h1p",
    "pfParticleNetAK4LastJetTags:probtaup1h2p",
    "pfParticleNetAK4LastJetTags:probtaup3h0p",
    "pfParticleNetAK4LastJetTags:probtaup3h1p",
    "pfParticleNetAK4LastJetTags:probtaum1h0p",
    "pfParticleNetAK4LastJetTags:probtaum1h1p",
    "pfParticleNetAK4LastJetTags:probtaum1h2p",
    "pfParticleNetAK4LastJetTags:probtaum3h0p",
    "pfParticleNetAK4LastJetTags:probtaum3h1p",
    "pfParticleNetAK4LastJetTags:probb",
    "pfParticleNetAK4LastJetTags:probc",
    "pfParticleNetAK4LastJetTags:probuds",
    "pfParticleNetAK4LastJetTags:probg",
    "pfParticleNetAK4LastJetTags:ptcorr",
    "pfParticleNetAK4LastJetTags:ptreslow",
    "pfParticleNetAK4LastJetTags:ptreshigh",
    "pfParticleNetAK4LastJetTags:ptnu",
    "pfParticleNetAK4LastNegativeJetTags:probmu",
    "pfParticleNetAK4LastNegativeJetTags:probele",
    "pfParticleNetAK4LastNegativeJetTags:probtaup1h0p",
    "pfParticleNetAK4LastNegativeJetTags:probtaup1h1p",
    "pfParticleNetAK4LastNegativeJetTags:probtaup1h2p",
    "pfParticleNetAK4LastNegativeJetTags:probtaup3h0p",
    "pfParticleNetAK4LastNegativeJetTags:probtaup3h1p",
    "pfParticleNetAK4LastNegativeJetTags:probtaum1h0p",
    "pfParticleNetAK4LastNegativeJetTags:probtaum1h1p",
    "pfParticleNetAK4LastNegativeJetTags:probtaum1h2p",
    "pfParticleNetAK4LastNegativeJetTags:probtaum3h0p",
    "pfParticleNetAK4LastNegativeJetTags:probtaum3h1p",
    "pfParticleNetAK4LastNegativeJetTags:probb",
    "pfParticleNetAK4LastNegativeJetTags:probc",
    "pfParticleNetAK4LastNegativeJetTags:probuds",
    "pfParticleNetAK4LastNegativeJetTags:probg",
    "pfParticleNetAK4LastNegativeJetTags:ptcorr",
    "pfParticleNetAK4LastNegativeJetTags:ptreslow",
    "pfParticleNetAK4LastNegativeJetTags:ptreshigh",
    "pfParticleNetAK4LastNegativeJetTags:ptnu",
    
])

pnetDiscriminatorLabels = [name.replace("pfParticleNetAK4LastJetTags:","").replace("pfParticleNetAK4LastNegativeJetTags:","neg") for name in pnetDiscriminatorNames]
## ParT training
parTDiscriminatorNames = [];
parTDiscriminatorLabels = [];

#evaluator from https://gitlab.cern.ch/rgerosa/particlenetstudiesrun2/-/blob/cmssw_13X_new_features/TrainingNtupleMakerAK4/plugins/ParTFeatureEvaluator.cc?ref_type=heads
#    
# included under PhysicsTools/NanoAOD/plugind
pfParTAK4LastJetTagInfos = cms.EDProducer("ParTFeatureEvaluator",
      muons = cms.InputTag("slimmedMuons"),
      electrons = cms.InputTag("slimmedElectrons"),
      photons = cms.InputTag("slimmedPhotons"),
      taus = cms.InputTag("slimmedTaus"),
      vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      secondary_vertices = cms.InputTag("slimmedSecondaryVertices"),
      #jets = cms.InputTag("slimmedJetsCalibrated" if options.applyJECs else ("slimmedJetsPuppi" if options.usePuppiJets else "slimmedJets")),
      jets = cms.InputTag("updatedJetsWithUserData"),# if options.applyJECs else ("slimmedJetsPuppi" if (options.usePuppiJets and not options.reRunPuppi) else ("patJetsPuppi" if (options.usePuppiJets and options.reRunPuppi) else "slimmedJets"))),
      losttracks = cms.InputTag("lostTracks"),
      jet_radius = cms.double(0.4),
      min_jet_pt = cms.double(10),
      max_jet_eta = cms.double(2.5),
      min_jet_eta = cms.double(0),
      min_pt_for_pfcandidates = cms.double(0.1),
      min_pt_for_track_properties = cms.double(-1),
      min_pt_for_losttrack = cms.double(1.0),
      max_dr_for_losttrack = cms.double(0.2),
      min_pt_for_taus = cms.double(20.),
      max_eta_for_taus = cms.double(2.5),
)
# only used in calib against fakes 
pfParTAK4LastNegativeJetTagInfos = pfParTAK4LastJetTagInfos.clone();
pfParTAK4LastNegativeJetTagInfos.flip_ip_sign = cms.bool(True);

import json
cmssw_base_dir = os.getenv("CMSSW_BASE");
with open(cmssw_base_dir+"/src/PhysicsTools/NanoAOD/data/ParTAK4/particle-transformer-3d.json") as json_data:
      jd = json.load(json_data)
      output_nodes = jd["output_names"]
      output_nodes = [node.replace("label","prob").replace("_","").replace("targetpt","ptcorr") for node in output_nodes]
  
pfParTAK4LastJetTags = boostedJetONNXJetTagsProducer.clone();
pfParTAK4LastJetTags.src = cms.InputTag("pfParTAK4LastJetTagInfos");
pfParTAK4LastJetTags.flav_names = cms.vstring(output_nodes);
pfParTAK4LastJetTags.preprocess_json = cms.string('PhysicsTools/NanoAOD/data/ParTAK4/particle-transformer-3d.json');
pfParTAK4LastJetTags.model_path = cms.FileInPath('PhysicsTools/NanoAOD/data/ParTAK4/particle-transformer-3d.onnx');

pfParTAK4LastNegativeJetTags = pfParTAK4LastJetTags.clone()
pfParTAK4LastNegativeJetTags.src = cms.InputTag("pfParTAK4LastNegativeJetTagInfos");

for node in output_nodes:
    parTDiscriminatorNames.append("pfParTAK4LastJetTags:"+node);
  
parTDiscriminatorLabels = [name.replace("pfParTAK4LastJetTags:","").replace("pfParTAK4LastNegativeJetTags:","neg") for name in parTDiscriminatorNames]
  

#print(parTDiscriminatorNames)
   

## Update final jet collection                                                                                                                                                    
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
slimmedJetsUpdated = updatedPatJets.clone(
    jetSource = "updatedJetsWithUserData",
    addJetCorrFactors = False,
)

slimmedJetsUpdated.discriminatorSources += pnetDiscriminatorNames        
slimmedJetsUpdated.discriminatorSources += parTDiscriminatorNames        
    

for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
    modifier.toModify( updatedJetsWithUserData.userInts,
            looseId = cms.InputTag("looseJetId"),
    )

#updatedJetsAK8WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
#     src = cms.InputTag("updatedJetsAK8"),
#     userFloats = cms.PSet(),
#     userInts = cms.PSet(
#        tightId = cms.InputTag("tightJetIdAK8"),
#        tightIdLepVeto = cms.InputTag("tightJetIdLepVetoAK8"),
#     ),
#)
#for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
#    modifier.toModify( updatedJetsAK8WithUserData.userInts,
#            looseId = cms.InputTag("looseJetIdAK8"),
#    )


#finalJets = cms.EDFilter("PATJetRefSelector",
#   src = cms.InputTag("updatedJetsWithUserData"),
#   cut = cms.string("pt > 0")#, eta<4.5, area>0.2, area<0.82
#
# add JES uncertainties BEFORE running JER smearing - is this correct?
slimmedJetsJESUpdated = cms.EDProducer(
    "CorrectionLibJECUpdater",
    jets=cms.InputTag("slimmedJetsUpdated"),
    correctionFile=cms.FileInPath("PhysicsTools/NanoAOD/data/jet_jerc.json")
)
# JER smearing
slimmedSmearedJets = cms.EDProducer('SmearedPATJetProducer',
       src = cms.InputTag('slimmedJetsJESUpdated'),
       enabled = cms.bool(True),
       rho = cms.InputTag("fixedGridRhoFastjetAll"),
       algo = cms.string('AK4PFchs'),
       algopt = cms.string('AK4PFchs_pt'),
       #resolutionFile = cms.FileInPath('Autumn18_V7_MC_PtResolution_AK4PFchs.txt'),
       #scaleFactorFile = cms.FileInPath('combined_SFs_uncertSources.txt'),

       genJets = cms.InputTag('slimmedGenJets'),
       dRMax = cms.double(0.2),
       dPtMaxFactor = cms.double(3),

       debug = cms.untracked.bool(False),
   # Systematic variation
   # 0: Nominal
   # -1: -1 sigma (down variation)
   # 1: +1 sigma (up variation)
   variation = cms.int32(0),  # If not specified, default to 0
   uncertaintySource = cms.string(""), # If not specified, default to Total
       )

# example for computing total uncertainties on jet resolution
slimmedSmaredJetsDown=slimmedSmearedJets.clone(variation=cms.int32(-1))
slimmedSmearedJetsUp=slimmedSmearedJets.clone(variation=cms.int32(1))

finalJets = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("slimmedJetsJESUpdated"),
    cut = cms.string("pt > 20 && abs(eta) < 2.5")
)

#finalJetsAK8 = cms.EDFilter("PATJetRefSelector",
#    src = cms.InputTag("updatedJetsAK8WithUserData"),
#    cut = cms.string("pt > 170")
#)





##################### Tables for final output and docs ##########################



jetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("linkedObjectsNew","jets"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Jet"),
    doc  = cms.string("slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (" + finalJets.cut.value()+")"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the jets
    externalVariables = cms.PSet(
        bReg2018 = ExtVar(cms.InputTag("bjetNN2018:corr"),float, doc="pt correction for b-jet energy regression",precision=12),
        #bRegMVA = ExtVar(cms.InputTag("bjetMVA"),float, doc="pt corrected with b-jet regression",precision=14),
        #bRegNN = ExtVar(cms.InputTag("bjetNN:corr"),float, doc="pt correction for b-jet energy regression",precision=12),
        #bRegNN2 = ExtVar(cms.InputTag("bjetNN2:corr"),float, doc="pt correction for b-jet energy regression",precision=12),
 #       bRegRes = ExtVar(cms.InputTag("bjetNN:res"),float, doc="res on pt corrected with b-jet regression",precision=8),
    ),
    variables = cms.PSet(P4Vars,
        area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
        nMuons = Var("?hasOverlaps('muons')?overlaps('muons').size():0", "uint8", doc="number of muons in the jet"),
        muonIdx1 = Var("?overlaps('muons').size()>0?overlaps('muons')[0].key():-1", int, doc="index of first matching muon"),
        muonIdx2 = Var("?overlaps('muons').size()>1?overlaps('muons')[1].key():-1", int, doc="index of second matching muon"),
        electronIdx1 = Var("?overlaps('electrons').size()>0?overlaps('electrons')[0].key():-1", int, doc="index of first matching electron"),
        electronIdx2 = Var("?overlaps('electrons').size()>1?overlaps('electrons')[1].key():-1", int, doc="index of second matching electron"),
        nElectrons = Var("?hasOverlaps('electrons')?overlaps('electrons').size():0", int, doc="number of electrons in the jet"),
        btagCMVA = Var("bDiscriminator('pfCombinedMVAV2BJetTags')",float,doc="CMVA V2 btag discriminator",precision=10),
        btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
        btagDeepFlavB = Var("bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')",float,doc="DeepFlavour b+bb+lepb tag discriminator",precision=10),
        btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
        btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",float,doc="DeepCSV charm btag discriminator",precision=10),
        btagDeepFlavC = Var("bDiscriminator('pfDeepFlavourJetTags:probc')",float,doc="DeepFlavour charm tag discriminator",precision=10),
	
	#Particle NET - weights from
	#https://gitlab.cern.ch/rgerosa/particlenetstudiesrun2/-/tree/cmssw_13X_new_features/TrainingNtupleMakerAK4/data/ParticleNetAK4/PNETUL/ClassRegTwoTargets?ref_type=heads
 
	btagPNetB = Var("?bDiscriminator('pfParticleNetAK4LastJetTags:probb')>0?bDiscriminator('pfParticleNetAK4LastJetTags:probb'):-1",float,precision=10,doc="ParticleNet probb"),     
	PNetRegPtRawCorr = Var("?abs(eta())<2.5?bDiscriminator('pfParticleNetAK4LastJetTags:ptcorr'):bDiscriminator('pfParticleNetAK4LastJetTags:ptcorr')",float,precision=10,doc="ParticleNet universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT"),
        PNetRegPtRawCorrNeutrino = Var("?abs(eta())<2.5?bDiscriminator('pfParticleNetAK4LastJetTags:ptnu'):bDiscriminator('pfParticleNetAK4LastJetTags:ptnu')",float,precision=10,doc="ParticleNet universal flavor-aware pT regression neutrino correction, relative to visible. To apply full regression, multiply raw jet pT by both PNetRegPtRawCorr and PNetRegPtRawCorrNeutrino."),
        PNetRegPtRawRes = Var("?abs(eta())<2.5?0.5*(bDiscriminator('pfParticleNetAK4LastJetTags:ptreshigh')-bDiscriminator('pfParticleNetAK4LastJetTags:ptreslow')):0.5*(bDiscriminator('pfParticleNetAK4LastJetTags:ptreshigh')-bDiscriminator('pfParticleNetAK4LastJetTags:ptreslow'))",float,precision=10,doc="ParticleNet universal flavor-aware jet pT resolution estimator, (q84 - q16)/2"),

	#UparTAK4  - weights from 
	#https://gitlab.cern.ch/rgerosa/particlenetstudiesrun2/-/tree/cmssw_13X_new_features/TrainingNtupleMakerAK4/data/ParTAK4?ref_type=heads

	tagUParTAK4B = Var("?bDiscriminator('pfParTAK4LastJetTags:probb')>0?bDiscriminator('pfParTAK4LastJetTags:probb'):-1",float,precision=10,doc="UnifiedParT b vs. udscg"),
	ParTAK4RegPtRawCorr = Var("?bDiscriminator('pfParTAK4LastJetTags:ptcorr')>0?bDiscriminator('pfParTAK4LastJetTags:ptcorr'):-1",float,precision=10,doc="UnifiedParT universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT"),
        UParTAK4RegPtRawCorrNeutrino = Var("?bDiscriminator('pfParTAK4LastJetTags:ptcorrnu')>0?bDiscriminator('pfParTAK4LastJetTags:ptcorrnu'):-1",float,precision=10,doc="UnifiedParT universal flavor-aware pT regression neutrino correction, relative to visible. To apply full regression, multiply raw jet pT by both UParTAK4RegPtRawCorr and UParTAK4RegPtRawCorrNeutrino."),
        UParTAK4RegPtRawRes = Var("?(bDiscriminator('pfParTAK4LastJetTags:ptcorrq84')+bDiscriminator('pfParTAK4LastJetTags:ptcorrq16'))>0?0.5*(bDiscriminator('pfParTAK4LastJetTags:ptcorrq84')-bDiscriminator('pfParTAK4LastJetTags:ptcorrq16')):-1",float,precision=10,doc="UnifiedParT universal flavor-aware jet pT resolution estimator, (q84 - q16)/2"),


	    #puIdDisc = Var("userFloat('pileupJetId:fullDiscriminant')",float,doc="Pilup ID discriminant",precision=10),
        puId = Var("userInt('pileupJetId:fullId')",int,doc="Pilup ID flags"),
        jetId = Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')",int,doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
        qgl = Var("userFloat('qgl')",float,doc="Quark vs Gluon likelihood discriminator",precision=10),
        nConstituents = Var("numberOfDaughters()",int,doc="Number of particles in the jet"),
        rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),
        chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
        neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
        chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
        neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
        muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
        #jercCHPUF = Var("userFloat('jercCHPUF')", float, doc="Pileup Charged Hadron Energy Fraction with the JERC group definition", precision= 6),
        #jercCHF = Var("userFloat('jercCHF')", float, doc="Charged Hadron Energy Fraction with the JERC group definition", precision= 6),
        
        #add variables computed in updatedJetsWithUserData
        leadTrackPt = Var("userFloat('leadTrackPt')", float, doc="pt of the leading track", precision= 6),  
        leptonPtRel = Var("userFloat('leptonPtRel')", float, doc="lepton pT relative to jet pT", precision= 6),
        leptonPtRatio = Var("userFloat('leptonPtRatio')", float, doc="Ratio of lepton pT to jet pT", precision= 6),
        leptonPtRelInv = Var("userFloat('leptonPtRelInv')", float, doc="Inverse of lepton pT relative to jet pT", precision= 6),
        leptonPtRelv0 = Var("userFloat('leptonPtRelv0')", float, doc="lepton pT relative to jet pT for v0", precision= 6),
        leptonPtRatiov0 = Var("userFloat('leptonPtRatiov0')", float, doc="Ratio of lepton pT to jet pT for v0", precision= 6),
        leptonPtRelInvv0 = Var("userFloat('leptonPtRelInvv0')", float, doc="Inverse of lepton pT relative to jet pT for v0", precision= 6),
        leptonDeltaR = Var("userFloat('leptonDeltaR')", float, doc="Delta R value for lepton", precision= 6),
        leptonPt = Var("userFloat('leptonPt')", float, doc="lepton pT", precision= 6),
        vtxPt = Var("userFloat('vtxPt')", float, doc="Vertex pT", precision= 6),
        vtxMass = Var("userFloat('vtxMass')", float, doc="Vertex mass", precision= 6),
        vtx3dL = Var("userFloat('vtx3dL')", float, doc="3D length of the vertex", precision= 6),
        vtx3deL = Var("userFloat('vtx3deL')", float, doc="3D length error of the vertex", precision= 6),
        ptD = Var("userFloat('ptD')", float, doc="pT D value", precision= 6),
        genPtwNu = Var("userFloat('genPtwNu')", float, doc="Generated pT with neutrinos", precision= 6),
        #add jec systematics 
        sys_JECAbsoluteStat_Up = Var("userFloat('jecSysAbsoluteStatUp')",float,doc="JEC uncertainty AbsoluteStat  up",precision=10),
    	sys_JECAbsoluteStat_Down = Var("userFloat('jecSysAbsoluteStatDown')",float,doc="JEC uncertainty AbsoluteStat down",precision=10),
        sys_JECAbsoluteScale_Up = Var("userFloat('jecSysAbsoluteScaleUp')",float,doc="JEC uncertainty AbsoluteScale  up",precision=10),
        sys_JECAbsoluteScale_Down = Var("userFloat('jecSysAbsoluteScaleDown')",float,doc="JEC uncertainty AbsoluteScale down",precision=10),
        sys_JECAbsoluteMPFBias_Up = Var("userFloat('jecSysAbsoluteMPFBiasUp')",float,doc="JEC uncertainty AbsoluteMPFBias  up",precision=10),
        sys_JECAbsoluteMPFBias_Down = Var("userFloat('jecSysAbsoluteMPFBiasDown')",float,doc="JEC uncertainty AbsoluteMPFBias down",precision=10),
        sys_JECFragmentation_Up = Var("userFloat('jecSysFragmentationUp')",float,doc="JEC uncertainty Fragmentation  up",precision=10),
        sys_JECFragmentation_Down = Var("userFloat('jecSysFragmentationDown')",float,doc="JEC uncertainty Fragmentation down",precision=10),
        sys_JECSinglePionECAL_Up = Var("userFloat('jecSysSinglePionECALUp')",float,doc="JEC uncertainty SinglePionECAL  up",precision=10),
        sys_JECSinglePionECAL_Down = Var("userFloat('jecSysSinglePionECALDown')",float,doc="JEC uncertainty SinglePionECAL down",precision=10),
        sys_JECSinglePionHCAL_Up = Var("userFloat('jecSysSinglePionHCALUp')",float,doc="JEC uncertainty SinglePionHCAL  up",precision=10),
        sys_JECSinglePionHCAL_Down = Var("userFloat('jecSysSinglePionHCALDown')",float,doc="JEC uncertainty SinglePionHCAL down",precision=10),
        sys_JECFlavorQCD_Up = Var("userFloat('jecSysFlavorQCDUp')",float,doc="JEC uncertainty FlavorQCD  up",precision=10),
        sys_JECFlavorQCD_Down = Var("userFloat('jecSysFlavorQCDDown')",float,doc="JEC uncertainty FlavorQCD down",precision=10),
        sys_JECTimePtEta_Up = Var("userFloat('jecSysTimePtEtaUp')",float,doc="JEC uncertainty TimePtEta  up",precision=10),
        sys_JECTimePtEta_Down = Var("userFloat('jecSysTimePtEtaDown')",float,doc="JEC uncertainty TimePtEta down",precision=10),
        sys_JECRelativeJEREC1_Up = Var("userFloat('jecSysRelativeJEREC1Up')",float,doc="JEC uncertainty RelativeJEREC1  up",precision=10),
        sys_JECRelativeJEREC1_Down = Var("userFloat('jecSysRelativeJEREC1Down')",float,doc="JEC uncertainty RelativeJEREC1 down",precision=10),
        sys_JECRelativeJEREC2_Up = Var("userFloat('jecSysRelativeJEREC2Up')",float,doc="JEC uncertainty RelativeJEREC2  up",precision=10),
        sys_JECRelativeJEREC2_Down = Var("userFloat('jecSysRelativeJEREC2Down')",float,doc="JEC uncertainty RelativeJEREC2 down",precision=10),
        sys_JECRelativeJERHF_Up = Var("userFloat('jecSysRelativeJERHFUp')",float,doc="JEC uncertainty RelativeJERHF  up",precision=10),
        sys_JECRelativeJERHF_Down = Var("userFloat('jecSysRelativeJERHFDown')",float,doc="JEC uncertainty RelativeJERHF down",precision=10),
        sys_JECRelativePtBB_Up = Var("userFloat('jecSysRelativePtBBUp')",float,doc="JEC uncertainty RelativePtBB  up",precision=10),
        sys_JECRelativePtBB_Down = Var("userFloat('jecSysRelativePtBBDown')",float,doc="JEC uncertainty RelativePtBB down",precision=10),
        sys_JECRelativePtEC1_Up = Var("userFloat('jecSysRelativePtEC1Up')",float,doc="JEC uncertainty RelativePtEC1  up",precision=10),
        sys_JECRelativePtEC1_Down = Var("userFloat('jecSysRelativePtEC1Down')",float,doc="JEC uncertainty RelativePtEC1 down",precision=10),
        sys_JECRelativePtEC2_Up = Var("userFloat('jecSysRelativePtEC2Up')",float,doc="JEC uncertainty RelativePtEC2  up",precision=10),
        sys_JECRelativePtEC2_Down = Var("userFloat('jecSysRelativePtEC2Down')",float,doc="JEC uncertainty RelativePtEC2 down",precision=10),
        sys_JECRelativePtHF_Up = Var("userFloat('jecSysRelativePtHFUp')",float,doc="JEC uncertainty RelativePtHF  up",precision=10),
        sys_JECRelativePtHF_Down = Var("userFloat('jecSysRelativePtHFDown')",float,doc="JEC uncertainty RelativePtHF down",precision=10),
        sys_JECRelativeBal_Up = Var("userFloat('jecSysRelativeBalUp')",float,doc="JEC uncertainty RelativeBal  up",precision=10),
        sys_JECRelativeBal_Down = Var("userFloat('jecSysRelativeBalDown')",float,doc="JEC uncertainty RelativeBal down",precision=10),
        sys_JECRelativeSample_Up = Var("userFloat('jecSysRelativeSampleUp')",float,doc="JEC uncertainty RelativeSample  up",precision=10),
        sys_JECRelativeSample_Down = Var("userFloat('jecSysRelativeSampleDown')",float,doc="JEC uncertainty RelativeSample down",precision=10),
        sys_JECRelativeStatEC_Up = Var("userFloat('jecSysRelativeStatECUp')",float,doc="JEC uncertainty RelativeStatEC up",precision=10),
        sys_JECRelativeStatEC_Down = Var("userFloat('jecSysRelativeStatECDown')",float,doc="JEC uncertainty RelativeStatEC up",precision=10),
        sys_JECRelativeStatFSR_Up = Var("userFloat('jecSysRelativeStatFSRUp')",float,doc="JEC uncertainty RelativeStatFSR up",precision=10),
        sys_JECRelativeStatFSR_Down = Var("userFloat('jecSysRelativeStatFSRUp')",float,doc="JEC uncertainty RelativeStatFSR up",precision=10),
        sys_JECRelativeStatHF_Up = Var("userFloat('jecSysRelativeStatHFUp')",float,doc="JEC uncertainty RelativeStatHF up",precision=10),
        sys_JECRelativeStatHF_Down = Var("userFloat('jecSysRelativeStatHFDown')",float,doc="JEC uncertainty RelativeStatHF up",precision=10),
        sys_JECRelativeFSR_Up = Var("userFloat('jecSysRelativeFSRUp')",float,doc="JEC uncertainty RelativeFSR  up",precision=10),
        sys_JECRelativeFSR_Down = Var("userFloat('jecSysRelativeFSRDown')",float,doc="JEC uncertainty RelativeFSR down",precision=10),
        sys_JECPileUpDataMC_Up = Var("userFloat('jecSysPileUpDataMCUp')",float,doc="JEC uncertainty PileUpDataMC  up",precision=10),
        sys_JECPileUpDataMC_Down = Var("userFloat('jecSysPileUpDataMCDown')",float,doc="JEC uncertainty PileUpDataMC down",precision=10),
        sys_JECPileUpPtRef_Up = Var("userFloat('jecSysPileUpPtRefUp')",float,doc="JEC uncertainty PileUpPtRef  up",precision=10),
        sys_JECPileUpPtRef_Down = Var("userFloat('jecSysPileUpPtRefDown')",float,doc="JEC uncertainty PileUpPtRef down",precision=10),
        sys_JECPileUpPtBB_Up = Var("userFloat('jecSysPileUpPtBBUp')",float,doc="JEC uncertainty PileUpPtBB  up",precision=10),
        sys_JECPileUpPtBB_Down = Var("userFloat('jecSysPileUpPtBBDown')",float,doc="JEC uncertainty PileUpPtBB down",precision=10),
        sys_JECPileUpPtEC1_Up = Var("userFloat('jecSysPileUpPtEC1Up')",float,doc="JEC uncertainty PileUpPtEC1  up",precision=10),
        sys_JECPileUpPtEC1_Down = Var("userFloat('jecSysPileUpPtEC1Down')",float,doc="JEC uncertainty PileUpPtEC1 down",precision=10),
        sys_JECPileUpPtEC2_Up = Var("userFloat('jecSysPileUpPtEC2Up')",float,doc="JEC uncertainty PileUpPtEC2  up",precision=10),
        sys_JECPileUpPtEC2_Down = Var("userFloat('jecSysPileUpPtEC2Down')",float,doc="JEC uncertainty PileUpPtEC2 down",precision=10),
        sys_JECPileUpPtHF_Up = Var("userFloat('jecSysPileUpPtHFUp')",float,doc="JEC uncertainty PileUpPtHF  up",precision=10),
        sys_JECPileUpPtHF_Down = Var("userFloat('jecSysPileUpPtHFDown')",float,doc="JEC uncertainty PileUpPtHF down",precision=10),
        


        # add the userInts as well
        tightId = Var("userInt('tightId')", int, doc="tightId flag", precision= 6),
        tightIdLepVeto = Var("userInt('tightIdLepVeto')", int, doc="tightIdLepVeto flag", precision= 6),
        vtxNtrk = Var("userInt('vtxNtrk')", int, doc="vtxNtrk flag", precision= 6),
        leptonPdgId = Var("userInt('leptonPdgId')", int, doc="leptonPdgId flag", precision= 6),   
    )
)

#jets are not as precise as muons
jetTable.variables.pt.precision=10



bjetNN2018= cms.EDProducer("BJetEnergyRegressionMVA",
    backend = cms.string("TF"),
    src = cms.InputTag("linkedObjectsNew","jets"),
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    rhosrc = cms.InputTag("fixedGridRhoFastjetAll"),

    #weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2017.pb"),         # 2017
    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/model-37.pb"),                  # 2018
    name = cms.string("JetRegNN"),
    isClassifier = cms.bool(False),
    variablesOrder = cms.vstring(["Jet_pt","Jet_eta","rho","Jet_mt","Jet_leadTrackPt","Jet_leptonPtRel","Jet_leptonDeltaR","Jet_neHEF","Jet_neEmEF","Jet_vtxPt","Jet_vtxMass","Jet_vtx3dL","Jet_vtxNtrk","Jet_vtx3deL","Jet_numDaughters_pt03","Jet_energyRing_dR0_em_Jet_rawEnergy","Jet_energyRing_dR1_em_Jet_rawEnergy","Jet_energyRing_dR2_em_Jet_rawEnergy","Jet_energyRing_dR3_em_Jet_rawEnergy","Jet_energyRing_dR4_em_Jet_rawEnergy","Jet_energyRing_dR0_neut_Jet_rawEnergy","Jet_energyRing_dR1_neut_Jet_rawEnergy","Jet_energyRing_dR2_neut_Jet_rawEnergy","Jet_energyRing_dR3_neut_Jet_rawEnergy","Jet_energyRing_dR4_neut_Jet_rawEnergy","Jet_energyRing_dR0_ch_Jet_rawEnergy","Jet_energyRing_dR1_ch_Jet_rawEnergy","Jet_energyRing_dR2_ch_Jet_rawEnergy","Jet_energyRing_dR3_ch_Jet_rawEnergy","Jet_energyRing_dR4_ch_Jet_rawEnergy","Jet_energyRing_dR0_mu_Jet_rawEnergy","Jet_energyRing_dR1_mu_Jet_rawEnergy","Jet_energyRing_dR2_mu_Jet_rawEnergy","Jet_energyRing_dR3_mu_Jet_rawEnergy","Jet_energyRing_dR4_mu_Jet_rawEnergy","Jet_chHEF","Jet_chEmEF","Jet_leptonPtRelInv","isEle","isMu","isOther","Jet_mass","Jet_ptd"]),
    variables = cms.PSet(
    Jet_pt = cms.string("pt*jecFactor('Uncorrected')"),
    Jet_mt = cms.string("mt*jecFactor('Uncorrected')"),
    Jet_eta = cms.string("eta"),
    Jet_mass = cms.string("mass*jecFactor('Uncorrected')"),
    Jet_ptd = cms.string("userFloat('ptD')"),
    Jet_leadTrackPt = cms.string("userFloat('leadTrackPt')"),
    Jet_vtxNtrk = cms.string("userInt('vtxNtrk')"),
    Jet_vtxMass = cms.string("userFloat('vtxMass')"),
    Jet_vtx3dL = cms.string("userFloat('vtx3dL')"),
    Jet_vtx3deL = cms.string("userFloat('vtx3deL')"),
    Jet_vtxPt = cms.string("userFloat('vtxPt')"),
    #Jet_leptonPt = cms.string("userFloat('leptonPt')"),
    Jet_leptonPtRel = cms.string("userFloat('leptonPtRelv0')"),
    Jet_leptonPtRelInv = cms.string("userFloat('leptonPtRelInvv0')*jecFactor('Uncorrected')"),
    Jet_leptonDeltaR = cms.string("userFloat('leptonDeltaR')"),
    #Jet_leptonPdgId = cms.string("userInt('leptonPdgId')"),
    Jet_neHEF = cms.string("neutralHadronEnergyFraction()"),
    Jet_neEmEF = cms.string("neutralEmEnergyFraction()"),
    Jet_chHEF = cms.string("chargedHadronEnergyFraction()"),
    Jet_chEmEF = cms.string("chargedEmEnergyFraction()"),
    isMu = cms.string("?abs(userInt('leptonPdgId'))==13?1:0"),
    isEle = cms.string("?abs(userInt('leptonPdgId'))==11?1:0"),
    isOther = cms.string("?userInt('leptonPdgId')==0?1:0"),
    ),
     inputTensorName = cms.string("ffwd_inp"),
     outputTensorName = cms.string("ffwd_out/BiasAdd"),
     outputNames = cms.vstring(["corr","res"]),
     #outputFormulas = cms.vstring(["at(0)*0.39077115058898926+1.0610932111740112","0.5*(at(2)-at(1))*0.39077115058898926"]),        # 2017 dont know the source
     outputFormulas = cms.vstring(["at(0)*0.27912887930870056+1.0545977354049683","0.5*(at(2)-at(1))*0.27912887930870056"]),       # 2018 training
     #nThreads = cms.uint32(1),
     #singleThreadPool = cms.string("no_threads"),
)


##### Soft Activity tables
#saJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#    src = cms.InputTag("softActivityJets"),
#    cut = cms.string(""),
#    maxLen = cms.uint32(6),
#    name = cms.string("SoftActivityJet"),
#    doc  = cms.string("jets clustered from charged candidates compatible with primary vertex (" + chsForSATkJets.cut.value()+")"),
#    singleton = cms.bool(False), # the number of entries is variable
#    extension = cms.bool(False), # this is the main table for the jets
#    variables = cms.PSet(P3Vars,
#  )
#)

#saJetTable.variables.pt.precision=10
#saJetTable.variables.eta.precision=8
#saJetTable.variables.phi.precision=8

#saTable = cms.EDProducer("GlobalVariablesTableProducer",
#    variables = cms.PSet(
#        SoftActivityJetHT = ExtVar( cms.InputTag("softActivityJets"), "candidatescalarsum", doc = "scalar sum of soft activity jet pt, pt>1" ),
#        SoftActivityJetHT10 = ExtVar( cms.InputTag("softActivityJets10"), "candidatescalarsum", doc = "scalar sum of soft activity jet pt , pt >10"  ),
#        SoftActivityJetHT5 = ExtVar( cms.InputTag("softActivityJets5"), "candidatescalarsum", doc = "scalar sum of soft activity jet pt, pt>5"  ),
#        SoftActivityJetHT2 = ExtVar( cms.InputTag("softActivityJets2"), "candidatescalarsum", doc = "scalar sum of soft activity jet pt, pt >2"  ),
#        SoftActivityJetNjets10 = ExtVar( cms.InputTag("softActivityJets10"), "candidatesize", doc = "number of soft activity jet pt, pt >2"  ),
#        SoftActivityJetNjets5 = ExtVar( cms.InputTag("softActivityJets5"), "candidatesize", doc = "number of soft activity jet pt, pt >5"  ),
#        SoftActivityJetNjets2 = ExtVar( cms.InputTag("softActivityJets2"), "candidatesize", doc = "number of soft activity jet pt, pt >10"  ),
#
#    )
#)



## BOOSTED STUFF #################
#fatJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#    src = cms.InputTag("finalJetsAK8"),
#    cut = cms.string(" pt > 170"), #probably already applied in miniaod
#    name = cms.string("FatJet"),
#    doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
#    singleton = cms.bool(False), # the number of entries is variable
#    extension = cms.bool(False), # this is the main table for the jets
#    variables = cms.PSet(P4Vars,
#        jetId = Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')",int,doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
#        area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
#        rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),
#        tau1 = Var("userFloat('NjettinessAK8Puppi:tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
#        tau2 = Var("userFloat('NjettinessAK8Puppi:tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
#        tau3 = Var("userFloat('NjettinessAK8Puppi:tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
#        tau4 = Var("userFloat('NjettinessAK8Puppi:tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
#       # n2b1 = Var("userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2')", float, doc="N2 with beta=1", precision=10),
#       # n3b1 = Var("userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3')", float, doc="N3 with beta=1", precision=10),
#        msoftdrop = Var("groomedMass('SoftDropPuppi')",float, doc="Corrected soft drop mass with PUPPI",precision=10),
#        btagCMVA = Var("bDiscriminator('pfCombinedMVAV2BJetTags')",float,doc="CMVA V2 btag discriminator",precision=10),
#        btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
#        btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
#        btagHbb = Var("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')",float,doc="Higgs to BB tagger discriminator",precision=10),
#        btagDDBvL = Var("bDiscriminator('pfMassIndependentDeepDoubleBvLJetTags:probHbb')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->bb vs QCD",precision=10),
#        btagDDCvL = Var("bDiscriminator('pfMassIndependentDeepDoubleCvLJetTags:probHcc')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->cc vs QCD",precision=10),
#        btagDDCvB = Var("bDiscriminator('pfMassIndependentDeepDoubleCvBJetTags:probHcc')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->cc vs H(Z)->bb",precision=10),
#        deepTag_TvsQCD = Var("bDiscriminator('pfDeepBoostedDiscriminatorsJetTags:TvsQCD')",float,doc="DeepBoostedJet tagger top vs QCD discriminator",precision=10),
#        deepTag_WvsQCD = Var("bDiscriminator('pfDeepBoostedDiscriminatorsJetTags:WvsQCD')",float,doc="DeepBoostedJet tagger W vs QCD discriminator",precision=10),
#        deepTag_ZvsQCD = Var("bDiscriminator('pfDeepBoostedDiscriminatorsJetTags:ZvsQCD')",float,doc="DeepBoostedJet tagger Z vs QCD discriminator",precision=10),
#        deepTag_H = Var("bDiscriminator('pfDeepBoostedJetTags:probHbb')+bDiscriminator('pfDeepBoostedJetTags:probHcc')+bDiscriminator('pfDeepBoostedJetTags:probHqqqq')",float,doc="DeepBoostedJet tagger H(bb,cc,4q) sum",precision=10),
#        deepTag_QCD = Var("bDiscriminator('pfDeepBoostedJetTags:probQCDbb')+bDiscriminator('pfDeepBoostedJetTags:probQCDcc')+bDiscriminator('pfDeepBoostedJetTags:probQCDb')+bDiscriminator('pfDeepBoostedJetTags:probQCDc')+bDiscriminator('pfDeepBoostedJetTags:probQCDothers')",float,doc="DeepBoostedJet tagger QCD(bb,cc,b,c,others) sum",precision=10),
#        deepTag_QCDothers = Var("bDiscriminator('pfDeepBoostedJetTags:probQCDothers')",float,doc="DeepBoostedJet tagger QCDothers value",precision=10),
#	deepTagMD_TvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger top vs QCD discriminator",precision=10),
#        deepTagMD_WvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger W vs QCD discriminator",precision=10),
#        deepTagMD_ZvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z vs QCD discriminator",precision=10),
#        deepTagMD_ZHbbvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z/H->bb vs QCD discriminator",precision=10),
#        deepTagMD_ZbbvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z->bb vs QCD discriminator",precision=10),
#        deepTagMD_HbbvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger H->bb vs QCD discriminator",precision=10),
#        deepTagMD_ZHccvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z/H->cc vs QCD discriminator",precision=10),
#        deepTagMD_H4qvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger H->4q vs QCD discriminator",precision=10),
#        deepTagMD_bbvsLight = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->bb vs light flavour discriminator",precision=10),
#        deepTagMD_ccvsLight = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight')",float,doc="Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->cc vs light flavour discriminator",precision=10),
#        subJetIdx1 = Var("?nSubjetCollections()>0 && subjets('SoftDropPuppi').size()>0?subjets('SoftDropPuppi')[0].key():-1", int,
#		     doc="index of first subjet"),
#        subJetIdx2 = Var("?nSubjetCollections()>0 && subjets('SoftDropPuppi').size()>1?subjets('SoftDropPuppi')[1].key():-1", int,
#		     doc="index of second subjet"),
#
##        btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",float,doc="CMVA V2 btag discriminator",precision=10),
##puIdDisc = Var("userFloat('pileupJetId:fullDiscriminant')",float,doc="Pilup ID discriminant",precision=10),
##        nConstituents = Var("numberOfDaughters()",int,doc="Number of particles in the jet"),
##        rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),
#    )
#)
#### Era dependent customization
##run2_miniAOD_80XLegacy.toModify( bjetNN,outputFormulas = cms.vstring(["at(0)*0.31628304719924927+1.0454729795455933","0.5*(at(2)-at(1))*0.31628304719924927"]))
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, msoftdrop_chs = Var("userFloat('ak8PFJetsCHSSoftDropMass')",float, doc="Legacy uncorrected soft drop mass with CHS",precision=10))
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables.tau1, expr = cms.string("userFloat(\'ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1\')"),)
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables.tau2, expr = cms.string("userFloat(\'ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2\')"),)
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables.tau3, expr = cms.string("userFloat(\'ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3\')"),)
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, tau4 = None)
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, n2b1 = None)
#run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, n3b1 = None)
#for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
#    modifier.toModify( fatJetTable.variables, jetId = Var("userInt('tightId')*2+userInt('looseId')",int,doc="Jet ID flags bit1 is loose, bit2 is tight"))
#
#
#
#
#subJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#    src = cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked","SubJets"),
#    cut = cms.string(""), #probably already applied in miniaod
#    name = cms.string("SubJet"),
#    doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
#    singleton = cms.bool(False), # the number of entries is variable
#    extension = cms.bool(False), # this is the main table for the jets
#    variables = cms.PSet(P4Vars,
#        btagCMVA = Var("bDiscriminator('pfCombinedMVAV2BJetTags')",float,doc="CMVA V2 btag discriminator",precision=10),
#        btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
#        btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
#        rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),                 
#        tau1 = Var("userFloat('NjettinessAK8Subjets:tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
#        tau2 = Var("userFloat('NjettinessAK8Subjets:tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
#        tau3 = Var("userFloat('NjettinessAK8Subjets:tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
#        tau4 = Var("userFloat('NjettinessAK8Subjets:tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
#        n2b1 = Var("userFloat('nb1AK8PuppiSoftDropSubjets:ecfN2')", float, doc="N2 with beta=1", precision=10),
#        n3b1 = Var("userFloat('nb1AK8PuppiSoftDropSubjets:ecfN3')", float, doc="N3 with beta=1", precision=10),
#    )
#)
#
##jets are not as precise as muons
#fatJetTable.variables.pt.precision=10
#subJetTable.variables.pt.precision=10
#
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, tau1 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, tau2 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, tau3 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, tau4 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, n2b1 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, n3b1 = None)
#run2_miniAOD_80XLegacy.toModify( subJetTable.variables, btagCMVA = None, btagDeepB = None)


corrT1METJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("corrT1METJets"),
    cut = cms.string(""),
    name = cms.string("CorrT1METJet"),
    doc  = cms.string("Additional low-pt jets for Type-1 MET re-correction"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the jets
    variables = cms.PSet(
        rawPt = Var("pt()*jecFactor('Uncorrected')",float,precision=10),
        eta  = Var("eta",  float,precision=12),
        phi = Var("phi", float, precision=12),
        area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
    )
)





## MC STUFF ######################
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJets #packed GenParts hold all particles with status 1 including visible
genParticlesForJets = genParticlesForJets.clone(
    src = cms.InputTag("packedGenParticles")
)

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
ak4GenJetsWithNu = ak4GenJets.clone(
    src = "genParticlesForJets"
)




jetMCTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("linkedObjectsNew","jets"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Jet"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), # this is an extension  table for the jets
    variables = cms.PSet(
        partonFlavour = Var("partonFlavour()", int, doc="flavour from parton matching"),
        hadronFlavour = Var("hadronFlavour()", int, doc="flavour from hadron ghost clustering"),
        genJetIdx = Var("?genJetFwdRef().backRef().isNonnull()?genJetFwdRef().backRef().key():-1", int, doc="index of matched gen jet"),
        
    )
)

#slimmedGenJets flavor info and association
genJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedGenJets"), #slimmed GenJets are built with slimmedGenParticles, holding only visible status =1 part + intermediate states 
    cut = cms.string("pt"),
    name = cms.string("GenJet"),
    doc  = cms.string("slimmedGenJets, i.e. ak4 Jets made with visible genparticles"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the genjets
    variables = cms.PSet(P4Vars,
	#anything else?
    )
)
patJetPartons = cms.EDProducer('HadronAndPartonSelector',
    src = cms.InputTag("generator"),
    particles = cms.InputTag("prunedGenParticles"),
    partonMode = cms.string("Auto"),
    fullChainPhysPartons = cms.bool(True)
)
genJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
    jets = genJetTable.src,
    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
    partons = cms.InputTag("patJetPartons","physicsPartons"),
    leptons = cms.InputTag("patJetPartons","leptons"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(False)
)
genJetFlavourTable = cms.EDProducer("GenJetFlavourTableProducer",
    name = genJetTable.name,
    src = genJetTable.src,
    cut = genJetTable.cut,
    deltaR = cms.double(0.1),
    jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
)
#GenJetsWithNuInfo
#print("after flav n")
genWNuJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("ak4GenJetsWithNu"), #slimmed GenJets are built with slimmedGenParticles, holding only visible status =1 part + intermediate states 
    cut = cms.string("pt"),
    name = cms.string("GenJetNu"),
    doc  = cms.string("packedGenJets, i.e. ak4 Jets made with all status =1  genparticles"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the genjets
    variables = cms.PSet(P4Vars,
	#anything else?
    )
)
print("after flav new jets association")
genWNuJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
    jets = genWNuJetTable.src,
    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
    partons = cms.InputTag("patJetPartons","physicsPartons"),
    leptons = cms.InputTag("patJetPartons","leptons"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(False)
)

genWNuJetFlavourTable = cms.EDProducer("GenJetFlavourTableProducer",
    name = genWNuJetTable.name,
    src = genWNuJetTable.src,
    cut = genWNuJetTable.cut,
    deltaR = cms.double(0.1),
    jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
)

#genJetAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
#    src = cms.InputTag("slimmedGenJetsAK8"),
#    cut = cms.string("pt > 100."),
#    name = cms.string("GenJetAK8"),
#    doc  = cms.string("slimmedGenJetsAK8, i.e. ak8 Jets made with visible genparticles"),
#    singleton = cms.bool(False), # the number of entries is variable
#    extension = cms.bool(False), # this is the main table for the genjets
#    variables = cms.PSet(P4Vars,
#	#anything else?
#    )
#)
#genJetAK8FlavourAssociation = cms.EDProducer("JetFlavourClustering",
#    jets = genJetAK8Table.src,
#    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
#    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
#    partons = cms.InputTag("patJetPartons","physicsPartons"),
#    leptons = cms.InputTag("patJetPartons","leptons"),
#    jetAlgorithm = cms.string("AntiKt"),
#    rParam = cms.double(0.8),
#    ghostRescaling = cms.double(1e-18),
#    hadronFlavourHasPriority = cms.bool(False)
#)
#genJetAK8FlavourTable = cms.EDProducer("GenJetFlavourTableProducer",
#    name = genJetAK8Table.name,
#    src = genJetAK8Table.src,
#    cut = genJetAK8Table.cut,
#    deltaR = cms.double(0.1),
#    jetFlavourInfos = cms.InputTag("genJetAK8FlavourAssociation"),
#)
#genSubJetAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
#    src = cms.InputTag("slimmedGenJetsAK8SoftDropSubJets"),
#    cut = cms.string(""),  ## These don't get a pt cut, but in miniAOD only subjets from fat jets with pt > 100 are kept
#    name = cms.string("SubGenJetAK8"),
#    doc  = cms.string("slimmedGenJetsAK8SoftDropSubJets, i.e. subjets of ak8 Jets made with visible genparticles"),
#    singleton = cms.bool(False), # the number of entries is variable
#    extension = cms.bool(False), # this is the main table for the genjets
#    variables = cms.PSet(P4Vars,
#	#anything else?
#    )
#)
### Era dependent customization
run2_miniAOD_80XLegacy.toModify( genJetFlavourTable, jetFlavourInfos = cms.InputTag("genJetFlavourAssociation"),)
run2_miniAOD_80XLegacy.toModify( genWNuJetFlavourTable, jetFlavourInfos = cms.InputTag("genJetFlavourAssociation"),)

from RecoJets.JetProducers.QGTagger_cfi import  QGTagger
qgtagger=QGTagger.clone(srcJets="updatedJets",srcVertexCollection="offlineSlimmedPrimaryVertices")

#before cross linking
jetSequence = cms.Sequence( jetCorrFactorsNano+updatedJets+ 
				tightJetId+tightJetIdLepVeto+bJetVars+jercVars+qgtagger+updatedJetsWithUserData+
				pfParticleNetAK4LastJetTagInfos+pfParticleNetAK4LastJetTags+
				pfParticleNetAK4LastNegativeJetTagInfos+pfParticleNetAK4LastNegativeJetTags+
				pfParTAK4LastJetTagInfos+pfParTAK4LastJetTags+
				pfParTAK4LastNegativeJetTagInfos+pfParTAK4LastNegativeJetTags+
				slimmedJetsUpdated+
                slimmedJetsJESUpdated+
				#   softActivityJets+softActivityJets2+softActivityJets5+softActivityJets10+
                          	finalJets)


#after cross linkining
jetTables = cms.Sequence(bjetNN2018+jetTable)#+subJetTable)

#MC only producers and tables
jetMC = cms.Sequence(jetMCTable+genJetTable+patJetPartons+genJetFlavourTable+genParticlesForJets+ak4GenJetsWithNu+genWNuJetTable+genWNuJetFlavourTable)