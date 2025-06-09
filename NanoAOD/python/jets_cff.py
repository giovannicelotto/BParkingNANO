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
from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
# Note: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
#      (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CMSSW_7_6_4_and_above )
jetCorrFactorsNano = patJetCorrFactors.clone(src='slimmedJets',
                                            levels = cms.vstring(   'L1FastJet',
                                                'L2Relative',
                                                'L3Absolute',
	                                            'L2L3Residual'),
                                            primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *

# JES nominal correction
# https://github.com/cms-sw/cmssw/blob/7a8aff4ce3c54f73b53e996aeb0c6369cd226d7b/RecoJets/JetProducers/python/PileupJetID_cfi.py#L18
#pileupJetId = cms.EDProducer('PileupJetIdProducer',
#     produceJetIds = cms.bool(True),
#     jetids = cms.InputTag(""),
#     runMvas = cms.bool(True),
#     jets = cms.InputTag("ak4PFJetsCHS"),
#     vertexes = cms.InputTag("offlinePrimaryVertices"),
#     algos = cms.VPSet(_stdalgos),
#     rho     = cms.InputTag("fixedGridRhoFastjetAll"),
#     jec     = cms.string("AK4PFchs"),
#     applyJec = cms.bool(True),
#     inputIsCorrected = cms.bool(False),
#     residualsFromTxt = cms.bool(False),
#     srcConstituentWeights = cms.InputTag(""),
##     residualsTxt     = cms.FileInPath("RecoJets/JetProducers/data/download.url") # must be an existing file
#)

# pileupJetIdNano=pileupJetId.clone(jets="updatedJets",algos = cms.VPSet(_chsalgos_106X_UL18),inputIsCorrected=True,applyJec=False,vertexes="offlineSlimmedPrimaryVertices")
#puIdNanoId = cms.InputTag('pileupJetIdNano:fullId'),
# JEc corretto
updatedJets = updatedPatJets.clone(
	addBTagInfo=False,
	jetSource='slimmedJets',
	jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsNano") ),
)

jecSources = [
    "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
    "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta",
    "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
    "RelativeBal", "RelativeSample", "RelativeFSR", "PileUpDataMC",
    "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"
]

# Jet ID to be applied on uncorrected Jets (neither JEC nor JES)
looseJetId = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('RUN2ULCHS'),
			    quality = cms.string('LOOSE'),
			  ),
                src = cms.InputTag("updatedJets")
)
tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('RUN2ULCHS'),
			    quality = cms.string('TIGHT'),
			  ),
            src = cms.InputTag("updatedJets")
)
tightJetIdLepVeto = cms.EDProducer("PatJetIDValueMapProducer",
			  filterParams=cms.PSet(
			    version = cms.string('RUN2ULCHS'),
			    quality = cms.string('TIGHTLEPVETO'),
			  ),
                src = cms.InputTag("updatedJets")
)



for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
    print("Here running change of versions")
    # check if running
    modifier.toModify( tightJetId.filterParams, version = "WINTER16" )
    modifier.toModify( tightJetIdLepVeto.filterParams, version = "WINTER16" )


# Updated Jets = JES corrected not JER
# Look here bjet applied to updated Jets
# https://github.com/cms-sw/cmssw/blob/CMSSW_13_3_X/PhysicsTools/NanoAOD/python/jetsAK4_CHS_cff.py
bJetVars = cms.EDProducer("JetRegressionVarProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("updatedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
    #musrc = cms.InputTag("slimmedMuons"),
    #elesrc = cms.InputTag("slimmedElectrons")
)

# Same for JercVars
jercVars = cms.EDProducer("BetaStarPackedCandidateVarProducer",
    srcJet = cms.InputTag("updatedJets"),    
    srcPF = cms.InputTag("packedPFCandidates"),
    maxDR = cms.double(0.4)
)

# variables for 2018 regression
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
         ),
     userInts = cms.PSet(
        tightId = cms.InputTag("tightJetId"),
        tightIdLepVeto = cms.InputTag("tightJetIdLepVeto"),
        vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
        leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
        
     ),
)     
    
## Update final jet collection                                                                                                                                                    
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
slimmedJetsUpdated = updatedPatJets.clone(
    jetSource = "updatedJetsWithUserData",
    addJetCorrFactors = False,
)


for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
    modifier.toModify( updatedJetsWithUserData.userInts,
            looseId = cms.InputTag("looseJetId"),
    )


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
slimmedSmearedJetsDown = slimmedSmearedJets.clone(variation=cms.int32(-1))
slimmedSmearedJetsUp = slimmedSmearedJets.clone(variation=cms.int32(1))

finalJets = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("slimmedSmearedJets"),
    #src = cms.InputTag("slimmedJetsJESUpdated"),
    cut = cms.string("pt > 20 && abs(eta) < 2.5")
)



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
        btagDeepFlavCvB = Var("?(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb'))>0?bDiscriminator('pfDeepFlavourJetTags:probc')/(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')):-1",float,doc="DeepJet c vs b+bb+lepb discriminator",precision=10),
        btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
        btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",float,doc="DeepCSV charm btag discriminator",precision=10),
        btagDeepFlavC = Var("bDiscriminator('pfDeepFlavourJetTags:probc')",float,doc="DeepFlavour charm tag discriminator",precision=10),

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
        # add jec systematics 

    # compact sources:
    # 
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
     outputFormulas = cms.vstring(["at(0)*0.27912887930870056+1.0545977354049683","0.5*(at(2)-at(1))*0.27912887930870056"]),       # 2018 training
)

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
        #genJetNuIdx = Var("?genJetNuFwdRef().backRef().isNonnull()?genJetNuFwdRef().backRef().key():-1", int, doc="index of matched gen jet Nu"),
        
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



### Era dependent customization
run2_miniAOD_80XLegacy.toModify( genJetFlavourTable, jetFlavourInfos = cms.InputTag("genJetFlavourAssociation"),)
run2_miniAOD_80XLegacy.toModify( genWNuJetFlavourTable, jetFlavourInfos = cms.InputTag("genJetFlavourAssociation"),)

from RecoJets.JetProducers.QGTagger_cfi import  QGTagger
qgtagger=QGTagger.clone(srcJets="updatedJets",srcVertexCollection="offlineSlimmedPrimaryVertices")

#before cross linking
jetSequence = cms.Sequence( jetCorrFactorsNano+updatedJets+ 
				tightJetId+tightJetIdLepVeto+bJetVars+jercVars+qgtagger+updatedJetsWithUserData+
				#pfParticleNetAK4LastJetTagInfos+pfParticleNetAK4LastJetTags+
				#pfParticleNetAK4LastNegativeJetTagInfos+pfParticleNetAK4LastNegativeJetTags+
				#pfParTAK4LastJetTagInfos+pfParTAK4LastJetTags+
				#pfParTAK4LastNegativeJetTagInfos+pfParTAK4LastNegativeJetTags+
				slimmedJetsUpdated+
                slimmedJetsJESUpdated+
                slimmedSmearedJets +
                #slimmedSmearedJetsUp +
                #slimmedSmearedJetsDown +
                

				#   softActivityJets+softActivityJets2+softActivityJets5+softActivityJets10+
                finalJets)


#after cross linkining
jetTables = cms.Sequence(bjetNN2018+jetTable)#+subJetTable)

#MC only producers and tables
jetMC = cms.Sequence(jetMCTable+genJetTable+patJetPartons+genJetFlavourTable+genParticlesForJets+ak4GenJetsWithNu+genWNuJetTable+genWNuJetFlavourTable)