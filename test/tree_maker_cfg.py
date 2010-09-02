import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/data/Run2010A/Mu/RECO/v4/000/144/114/9C954151-32B4-DF11-BB88-001D09F27003.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/114/5C9CA515-20B4-DF11-9D62-0030487A3DE0.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/114/00CA69A9-1EB4-DF11-B869-0030487CD6D2.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/FE61BD50-CAB3-DF11-B7FD-001D09F29321.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/FCEBD143-D1B3-DF11-854A-001D09F2423B.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TFileService = cms.Service(
    "TFileService",

    fileName = cms.string("top_tree.root")
)

process.load("Top.EDAnalyzers.TreeMaker_cfi")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

process.load("Top.Filters.HLTFilter_cfi")
process.HLTFilter.hltTag = cms.string("TriggerResults::HLT")
process.load("Top.Filters.PVFilter_cfi")

process.p = cms.Path(
    process.HLTFilter *
    process.PVFilter *
    process.ak5CaloJetsL2L3 *
    process.metCorSequence *
    process.TreeMaker
)

import FWCore.ParameterSet.printPaths as pp      
pp.printPaths(process)
