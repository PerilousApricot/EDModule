import FWCore.ParameterSet.Config as cms

HLTFilter = cms.EDFilter(
    'HLTFilter',

    hltTag = cms.string("TriggerResults::REDIGI"),
    hltName = cms.string("HLT_Mu9")
)
