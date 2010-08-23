/**
 * TreeMaker
 * 
 *
 * Created by Samvel Khalatian on August 20, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Top/Tree/interface/Electron.h"
#include "Top/Tree/interface/Jet.h"
#include "Top/Tree/interface/Muon.h"

#include "Top/EDAnalyzers/interface/TreeMaker.h"

using std::cout;
using std::endl;
using std::string;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;

void setP4(top::LorentzVector *topP4,
           const math::XYZTLorentzVector *cmsswP4)
{
    topP4->SetPxPyPzE(cmsswP4->px(),
                      cmsswP4->py(),
                      cmsswP4->pz(),
                      cmsswP4->energy());
}

void setP4(top::LorentzVector *topP4,
           const math::XYZTLorentzVector &cmsswP4)
{
    setP4(topP4, &cmsswP4);
}

TreeMaker::TreeMaker(const edm::ParameterSet &config)
{
    _metTag = config.getParameter<string>("metTag");
    _muonTag = config.getParameter<string>("muonTag");
    _jetTag = config.getParameter<string>("jetTag");
    _electronTag = config.getParameter<string>("electronTag");
}

TreeMaker::~TreeMaker()
{
}

void TreeMaker::beginJob()
{
    _topFile.reset(new TFile("top_tree.root", "RECREATE"));

    if (!_topFile->IsOpen())
    {
        throw edm::Exception(ErrorCodes(8001))
            << "failed to create output file." << endl;

        _topFile.reset();

        return;
    }

    _topEvent.reset(new top::Event());
    _topTree.reset(new TTree("top", "Top TTmuj Tree."));
    _topTree->Branch("event", _topEvent.get(), 32000, 0);
}

void TreeMaker::endJob()
{
    if (!_topEvent.get())
        return;

    _topFile->Write();

    _topTree.reset();

    _topFile.reset();

    // Note: Event should be destroyed after ROOT file is written and closed.
    _topEvent.reset();
}

void TreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    if (!_topEvent.get())
        return;

    // Extract MET
    Handle<CaloMETCollection> caloMets;
    event.getByLabel(InputTag(_metTag), caloMets);

    if (!caloMets.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Calo METs.";

        return;
    }

    // Extract Muons
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muonTag), muons);

    if (!muons.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Muons.";

        return;
    }

    // Extract Jets
    Handle<CaloJetCollection> jets;
    event.getByLabel(InputTag(_jetTag), jets);

    if (!jets.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Jets.";

        return;
    }

    // Extract Electrons
    Handle<GsfElectronCollection> electrons;
    event.getByLabel(InputTag(_electronTag), electrons);

    if (!electrons.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract Electrons.";

        return;
    }

    _topEvent->reset();

    _topEvent->id()->setRun(event.id().run());
    _topEvent->id()->setLumiBlock(event.id().luminosityBlock());
    _topEvent->id()->setEvent(event.id().event());

    CaloMETCollection::const_iterator met = caloMets->begin();
    setP4(_topEvent->met()->p4(), met->p4());

    // Process all Muons
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        top::Muon topMuon;
        topMuon.setIsGlobal(muon->isGlobalMuon());
        topMuon.setIsTracker(muon->isTrackerMuon());

        setP4(topMuon.p4(), muon->p4());

        _topEvent->muons()->push_back(topMuon);
    }

    // Process All Jets
    for(CaloJetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        top::Jet topJet;

        setP4(topJet.p4(), jet->p4());

        _topEvent->jets()->push_back(topJet);
    }

    // Process All Electrons
    for(GsfElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        top::Electron topElectron;

        setP4(topElectron.p4(), electron->p4());

        _topEvent->electrons()->push_back(topElectron);
    }

    _topTree->Fill();
}

DEFINE_FWK_MODULE(TreeMaker);
