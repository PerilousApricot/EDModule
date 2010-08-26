import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/Spring10/ZJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0014/F0D822C6-4D47-DF11-9590-0030487E5179.root',
        '/store/mc/Spring10/ZJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0014/E4D01BAA-B647-DF11-A2C8-0030488A1188.root',
        '/store/mc/Spring10/ZJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0014/9AF1B5DE-5247-DF11-8125-003048D4363C.root',
        '/store/mc/Spring10/ZJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0014/4075725A-5547-DF11-AFC6-0030487D5D67.root',
        '/store/mc/Spring10/ZJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0014/3C600AD4-5247-DF11-A57E-0030487D5D67.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("Top.EDAnalyzers.TreeMaker_cfi")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

process.load("Top.Filters.HLTFilter_cfi")
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
