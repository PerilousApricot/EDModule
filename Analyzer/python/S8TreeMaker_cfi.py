import FWCore.ParameterSet.Config as cms

S8TreeMaker = cms.EDAnalyzer(
    'S8TreeMaker',

    primaryVertices = cms.string("offlinePrimaryVertices"),
    muons = cms.string("selectedPatMuonsForPtRel"),
    jets = cms.string("selectedPatJetsAK5PF"),

    isPythia = cms.bool(False)
)
