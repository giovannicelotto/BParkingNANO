import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

beamSpotTable = cms.EDProducer("SimpleBeamspotFlatTableProducer",
    src = cms.InputTag("offlineBeamSpot"),
    name = cms.string("BeamSpot"),
    doc = cms.string("offlineBeamSpot, the offline reconstructed beamspot"),
    singleton = cms.bool(True),  # there's always exactly one MET per event
    extension = cms.bool(False), # this is the main table for the MET
    variables = cms.PSet(
       type = Var("type()","int8",doc="BeamSpot type (Unknown = -1, Fake = 0, LHC = 1, Tracker = 2)"),
       z = Var("position().z()",float,doc="BeamSpot center, z coordinate (cm)",precision=-1),
       zError = Var("z0Error()",float,doc="Error on BeamSpot center, z coordinate (cm)",precision=-1),
       sigmaZ = Var("sigmaZ()",float,doc="Width of BeamSpot in z (cm)",precision=-1),
       sigmaZError = Var("sigmaZ0Error()",float,doc="Error on width of BeamSpot in z (cm)",precision=-1),
    ),
)

rhoTable = cms.EDProducer("GlobalVariablesTableProducer",
    name = cms.string("Rho"),
    variables = cms.PSet(
        # removed GC
        #fixedGridRhoAll = ExtVar( cms.InputTag("fixedGridRhoAll"), "double", doc = "rho from all PF Candidates, no foreground removal (for isolation of prompt photons)" ),
        #fixedGridRhoFastjetAll = ExtVar( cms.InputTag("fixedGridRhoFastjetAll"), "double", doc = "rho from all PF Candidates, used e.g. for JECs" ),
        #fixedGridRhoFastjetCentralNeutral = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralNeutral"), "double", doc = "rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations" ),
        #fixedGridRhoFastjetCentralCalo = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralCalo"), "double", doc = "rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation" ),
        #fixedGridRhoFastjetCentral = ExtVar( cms.InputTag("fixedGridRhoFastjetCentral"), "double", doc = "rho from all PF Candidates for central region, used e.g. for JECs" ),
        #fixedGridRhoFastjetCentralChargedPileUp = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"), "double", doc = "rho from charged PF Candidates for central region, used e.g. for JECs" ),
    )
)

puTable = cms.EDProducer("NPUTablesProducer",
        src = cms.InputTag("slimmedAddPileupInfo"),
        pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
        zbins = cms.vdouble( [0.0,1.7,2.6,3.0,3.5,4.2,5.2,6.0,7.5,9.0,12.0] ),
        savePtHatMax = cms.bool(False), 
)

genTable  = cms.EDProducer("SimpleGenEventFlatTableProducer",
        src = cms.InputTag("generator"),
        cut = cms.string(""), 
        name= cms.string("Generator"),
        doc = cms.string("Generator information"),
        singleton = cms.bool(True), 
        extension = cms.bool(False),
    variables = cms.PSet(
        # commented GC
        #x1 = Var( "?hasPDF?pdf().x.first:-1", float, doc="x1 fraction of proton momentum carried by the first parton",precision=14 ),
        #x2 = Var( "?hasPDF?pdf().x.second:-1", float, doc="x2 fraction of proton momentum carried by the second parton",precision=14 ),
        #xpdf1 = Var( "?hasPDF?pdf().xPDF.first:-1", float, doc="x*pdf(x) for the first parton", precision=14 ),
        #xpdf2 = Var( "?hasPDF?pdf().xPDF.second:-1", float, doc="x*pdf(x) for the second parton", precision=14 ),
        #id1 = Var( "?hasPDF?pdf().id.first:-1", int, doc="id of first parton", precision=6 ),
        #id2 = Var( "?hasPDF?pdf().id.second:-1", int, doc="id of second parton", precision=6 ),
        #scalePDF = Var( "?hasPDF?pdf().scalePDF:-1", float, doc="Q2 scale for PDF", precision=14 ),
        #binvar = Var("?hasBinningValues()?binningValues()[0]:-1", float, doc="MC generation binning value", precision=14),
        weight = Var("weight()", float,doc="MC generator weight", precision=14),
        ),
)

genFilterTable = cms.EDProducer("SimpleGenFilterFlatTableProducerLumi",
        src = cms.InputTag("genFilterEfficiencyProducer"),
        cut = cms.string(""), 
        name= cms.string("GenFilter"),
        doc = cms.string("Generator filter information"),
        singleton = cms.bool(True), 
        extension = cms.bool(False),
    variables = cms.PSet(
        #sumPassWeights        = Var("sumPassWeights()",            "float", doc="Sum of passed weights of all visited events",      precision=23),
        #sumFailWeights        = Var("sumFailWeights()",            "float", doc="Sum of failed weights of all visited events",      precision=23),
        #sumTotalWeights_      = Var("sumTotalWeights_()",            "float", doc="Sum of weights of all visited events",      precision=23),
        sumw                  = Var("sumWeights()",            "float", doc="Sum of weights of all visited events",      precision=23),
        numEventsTotal        = Var("numEventsTotal()",        int,   doc="generator filter: total number of events",  precision=6),
        numEventsPassed       = Var("numEventsPassed()",       int,   doc="generator filter: passed number of events", precision=6),
        filterEfficiency      = Var("filterEfficiency()",      float, doc="generator filter: efficiency",              precision=23),
        filterEfficiencyError = Var("filterEfficiencyError()", float, doc="generator filter: efficiency error",        precision=23),
        ),
)

globalTablesTask = cms.Task(beamSpotTable, rhoTable)
globalTablesMCTask = cms.Task(puTable,genTable,genFilterTable)
globalTablesMC = puTable+genTable+genFilterTable
globalTables = beamSpotTable + rhoTable