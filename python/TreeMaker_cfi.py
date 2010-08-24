import FWCore.ParameterSet.Config as cms

TreeMaker = cms.EDFilter(
    'TreeMaker',

    metTag = cms.string("metMuonJESCorAK5"),
    muonTag = cms.string("muons"),
    jetTag = cms.string("ak5CaloJetsL2L3"),
    electronTag = cms.string("gsfElectrons"),
    beamSpotTag = cms.string("offlineBeamSpot"),

    jetIDParams = cms.PSet(
        useRecHits = cms.bool(True),
        hbheRecHitsColl = cms.InputTag("hbhereco"),
        hoRecHitsColl   = cms.InputTag("horeco"),
        hfRecHitsColl   = cms.InputTag("hfreco"),
        ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
        eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
    )
)
