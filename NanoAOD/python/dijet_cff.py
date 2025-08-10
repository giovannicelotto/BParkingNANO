import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var


dijetBuilder = cms.EDProducer("DijetBuilderProducer",
    jets = cms.InputTag("linkedObjectsNew","jets"),
    muons = cms.InputTag("linkedObjectsNew","muons"),
    bRegCorr = cms.InputTag("bjetNN2018", "corr")
)
dijetPtFilter = cms.EDFilter("DijetPtFilter",
    src = cms.InputTag("dijetBuilder"),
    minPt = cms.double(90.0)
)
dijetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("dijetBuilder"),
    cut = cms.string(""),
    name = cms.string("dijet"),
    doc = cms.string("Dijet system variables"),
    singleton = cms.bool(True),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt  = Var("pt",  float, doc="pt of dijet", precision=10),
        eta = Var("eta", float, doc="eta of dijet", precision=10),
        phi = Var("phi", float, doc="phi of dijet", precision=10),
        mass= Var("mass",float, doc="mass of dijet", precision=10),
    )
)

# Export dictionary so importing file can attach to process
dijetSequence = cms.Sequence(dijetBuilder  + dijetPtFilter+  dijetTable)
