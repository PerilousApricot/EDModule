import FWCore.ParameterSet.Config as cms

NtupleMaker = cms.EDAnalyzer(
    "NtupleMaker",

    metCollection      = cms.string("patMETs"),
    electronCollection = cms.string("selectedPatElectrons"),
    muonCollection     = cms.string("selectedPatMuons"),
    jetCollection      = cms.string("selectedPatJets"),
    pvCollection       = cms.string("offlinePrimaryVertices"),
    beamSpot           = cms.string("offlineBeamSpot")
)
