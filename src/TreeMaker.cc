/**
 * TreeMaker
 * 
 *
 * Created by Samvel Khalatian on August 20, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"

#include "Top/EDAnalyzers/interface/Tools.h"
#include "Top/Tree/interface/Electron.h"
#include "Top/Tree/interface/ElectronIsolation.h"
#include "Top/Tree/interface/EventID.h"
#include "Top/Tree/interface/Jet.h"
#include "Top/Tree/interface/JetEnergy.h"
#include "Top/Tree/interface/Muon.h"
#include "Top/Tree/interface/MuonIsolation.h"

#include "Top/EDAnalyzers/interface/TreeMaker.h"

using std::cout;
using std::endl;
using std::string;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;
using edm::ParameterSet;

using reco::BeamSpot;
using reco::helper::JetIDHelper;

TreeMaker::TreeMaker(const edm::ParameterSet &config)
{
    _metTag = config.getParameter<string>("metTag");
    _muonTag = config.getParameter<string>("muonTag");
    _jetTag = config.getParameter<string>("jetTag");
    _electronTag = config.getParameter<string>("electronTag");
    _beamSpotTag = config.getParameter<string>("beamSpotTag");

    _jetID = new JetIDHelper(config.getParameter<ParameterSet>("jetIDParams"));
}

TreeMaker::~TreeMaker()
{
    delete _jetID;
}

void TreeMaker::beginJob()
{
    edm::Service<TFileService> fileService;

    _topTree = fileService->make<TTree>("top", "Top ttmuj tree.");

    _topEvent.reset(new top::Event());
    _topTree->Branch("event", _topEvent.get(), 32000, 0);
}

void TreeMaker::endJob()
{
    if (!_topEvent.get())
        return;

    // Note: Event should be destroyed after ROOT file is written and closed.
    _topEvent.reset();
}

void TreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    if (!_topEvent.get())
        return;

    using namespace top::tools;

    // Extract BeamSpot
    Handle<BeamSpot> beamSpot;
    event.getByLabel(InputTag(_beamSpotTag), beamSpot);

    if (!beamSpot.isValid())
    {
        LogWarning("TreeMaker")
            << "failed to extract BeamSpot.";

        return;
    }

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

    {
        top::EventID id;
        id.setRun(event.id().run());
        id.setLumiBlock(event.id().luminosityBlock());
        id.setEvent(event.id().event());

        _topEvent->setID(id);
    }

    CaloMETCollection::const_iterator met = caloMets->begin();
    setP4(_topEvent->met().p4(), met->p4());

    // Process all Muons
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // Only Muons with basic cuts are saved:
        if (muon->isGlobalMuon() &&
            10 <= muon->pt() &&
            2.5 >= fabs(muon->eta()))
        {
            top::Muon topMuon;
            topMuon.setIsGlobal(muon->isGlobalMuon());
            topMuon.setIsTracker(muon->isTrackerMuon());

            topMuon.setEta(muon->eta());
            topMuon.setPhi(muon->phi());

            topMuon.setMatches(muon->numberOfMatches());

            setP4(topMuon.p4(), muon->p4());

            setIsolation(topMuon, top::Muon::R03, muon->isolationR03());
            setIsolation(topMuon, top::Muon::R05, muon->isolationR05());

            // Inner Track is only available for the Tracker Muons
            if (muon->isTrackerMuon())
            {
                top::ImpactParameter ip;

                ip.setValue(muon->innerTrack()->dxy(beamSpot->position()));

                ip.setError(sqrt(muon->innerTrack()->d0Error() *
                                 muon->innerTrack()->d0Error() +
                                 0.5 *
                                 beamSpot->BeamWidthX() *
                                 beamSpot->BeamWidthX() +
                                 0.5 *
                                 beamSpot->BeamWidthY() *
                                 beamSpot->BeamWidthY()));

                topMuon.setImpactParameter(top::Muon::BS2D, ip);

                topMuon.setInnerValidHits(muon->innerTrack()->numberOfValidHits());
            }

            // Next properties are only defined for the Global Muon
            topMuon.setOuterValidHits(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
            topMuon.setChi2(muon->globalTrack()->chi2());
            topMuon.setNdof(muon->globalTrack()->ndof());

            _topEvent->muons().push_back(topMuon);
        }
    } // End loop over muons

    // Process All Jets
    for(CaloJetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        // Select only energetic jets
        if (30 <= jet->pt() &&
            2.4 >= fabs(jet->eta()))
        {
            _jetID->calculate(event, *jet);

            top::Jet topJet;

            setP4(topJet.p4(), jet->p4());
            setEnergy(topJet, jet->getSpecific());

            topJet.setHits90(_jetID->n90Hits());
            topJet.setHpd(_jetID->fHPD());

            _topEvent->jets().push_back(topJet);
        }
    }

    // Process All Electrons
    for(GsfElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        // Select electrons to be stored
        if (15 <= electron->et() &&
            2.5 >= fabs(electron->eta()))
        {
            top::Electron topElectron;

            setP4(topElectron.p4(), electron->p4());

            setIsolation(topElectron, top::Electron::R03, electron->isolationVariables03());
            setIsolation(topElectron, top::Electron::R04, electron->isolationVariables04());

            _topEvent->electrons().push_back(topElectron);
        }
    }

    _topTree->Fill();
}

DEFINE_FWK_MODULE(TreeMaker);
