import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:./ttmuj.root'
    )
)

process.load("Top.EDAnalyzers.NtupleMaker_cfi")

process.p = cms.Path(process.NtupleMaker)
