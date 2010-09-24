import FWCore.ParameterSet.Config as cms

PATTreeMaker = cms.EDFilter(
    'PATTreeMaker',

    metTag = cms.string("patMETs"),
    muonTag = cms.string("selectedPatMuons"),
    jetTag = cms.string("selectedPatJets"),
    electronTag = cms.string("selectedPatElectrons"),
    beamSpotTag = cms.string("offlineBeamSpot"),
)
