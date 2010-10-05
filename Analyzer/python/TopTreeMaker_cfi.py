import FWCore.ParameterSet.Config as cms

TreeMaker = cms.EDFilter(
    'TreeMaker',

    mets = cms.string("metMuonJESCorAK5"),
    muons = cms.string("muons"),
    jets = cms.string("ak5CaloJetsL2L3"),
    electrons = cms.string("gsfElectrons"),
    beamSpots = cms.string("offlineBeamSpot")
)
