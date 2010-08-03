import FWCore.ParameterSet.Config as cms

NtupleMaker = cms.EDAnalyzer(
    "NtupleMaker",

    metCollection      = cms.InputTag("selectedPatJets"),
    electronCollection = cms.InputTag("selectedPatElectrons"),
    muonCollection     = cms.InputTag("selectedPatMuons"),
    jetCollection      = cms.InputTag("selectedPatJets"),
    pvCollection       = cms.InputTag("offlinePrimaryVertices"),
    beamSpot           = cms.InputTag("offlineBeamSpot")
)
