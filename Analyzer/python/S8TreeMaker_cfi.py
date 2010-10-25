import FWCore.ParameterSet.Config as cms

S8TreeMaker = cms.EDAnalyzer(
    'S8TreeMaker',

    primaryVertices = cms.string("offlinePrimaryVertices"),
    muons = cms.string("selectedPatMuonsForPtRel"),
    electrons = cms.string("selectedPatElectronsForS8"),
    jets = cms.string("selectedPatJetsAK5PF"),
    triggers = cms.string("TriggerResults::REDIGI"),

    inputType = cms.string("undefined")
)
