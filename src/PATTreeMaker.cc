/**
 * PATTreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Top/EDAnalyzers/interface/Tools.h"
#include "Top/Tree/interface/Electron.h"
#include "Top/Tree/interface/ElectronIsolation.h"
#include "Top/Tree/interface/Event.h"
#include "Top/Tree/interface/EventID.h"
#include "Top/Tree/interface/Jet.h"
#include "Top/Tree/interface/JetEnergy.h"
#include "Top/Tree/interface/Muon.h"
#include "Top/Tree/interface/MuonIsolation.h"

#include "Top/EDAnalyzers/interface/PATTreeMaker.h"

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

PATTreeMaker::PATTreeMaker(const edm::ParameterSet &config)
{
    _topEvent = 0;
    _topTree = 0;

    _metTag = config.getParameter<string>("metTag");
    _muonTag = config.getParameter<string>("muonTag");
    _jetTag = config.getParameter<string>("jetTag");
    _electronTag = config.getParameter<string>("electronTag");
    _beamSpotTag = config.getParameter<string>("beamSpotTag");
}

PATTreeMaker::~PATTreeMaker()
{
    if (_topEvent)
    {
        // Note: Event should be destroyed after ROOT file is written and closed.
        delete _topEvent;
    }
}

void PATTreeMaker::beginJob()
{
    edm::Service<TFileService> fileService;

    _topTree = fileService->make<TTree>("top", "Top ttmuj tree.");

    _topEvent = new top::Event();
    _topTree->Branch("event", _topEvent, 32000, 0);
}

void PATTreeMaker::endJob()
{
}

void PATTreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    using pat::ElectronCollection;
    using pat::JetCollection;
    using pat::MuonCollection;
    using pat::METCollection;

    using namespace top::tools;

    using top::EventID;

    if (!_topEvent)
        return;

    // Extract BeamSpot
    Handle<BeamSpot> beamSpot;
    event.getByLabel(InputTag(_beamSpotTag), beamSpot);

    if (!beamSpot.isValid())
    {
        LogWarning("PATTreeMaker")
            << "failed to extract BeamSpot.";

        return;
    }

    // Extract MET
    Handle<METCollection> mets;
    event.getByLabel(InputTag(_metTag), mets);

    if (!mets.isValid())
    {
        LogWarning("PATTreeMaker")
            << "failed to extract Calo METs.";

        return;
    }

    // Extract Muons
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muonTag), muons);

    if (!muons.isValid())
    {
        LogWarning("PATTreeMaker")
            << "failed to extract Muons.";

        return;
    }

    // Extract Jets
    Handle<JetCollection> jets;
    event.getByLabel(InputTag(_jetTag), jets);

    if (!jets.isValid())
    {
        LogWarning("PATTreeMaker")
            << "failed to extract Jets.";

        return;
    }

    // Extract Electrons
    Handle<ElectronCollection> electrons;
    event.getByLabel(InputTag(_electronTag), electrons);

    if (!electrons.isValid())
    {
        LogWarning("PATTreeMaker")
            << "failed to extract Electrons.";

        return;
    }

    _topEvent->reset();

    {
        EventID id;
        id.setRun(event.id().run());
        id.setLumiBlock(event.id().luminosityBlock());
        id.setEvent(event.id().event());

        _topEvent->setID(id);
    }

/*
    // Use only first MET
    METCollection::const_iterator met = mets->begin();
    setP4(_topEvent->met().p4(), met->p4());

    // Process all Muons
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // At lease very loose Global Muons are saved
        //
        if (muon->isGlobalMuon() &&
            10 <= muon->pt() &&
            2.5 >= fabs(muon->eta()))
        {
            // Create top::Muon
            top::Muon topMuon;
            topMuon.setIsGlobal(muon->isGlobalMuon());
            topMuon.setIsTracker(muon->isTrackerMuon());

            topMuon.setMatches(muon->numberOfMatches());

            setP4(topMuon.p4(), muon->p4());

            setIsolation(topMuon.isolation(top::Muon::R03), muon->isolationR03());
            setIsolation(topMuon.isolation(top::Muon::R05), muon->isolationR05());

            // Inner Track is only available for the Tracker Muons
            if (muon->isTrackerMuon())
            {
                top::ImpactParameter *ip = topMuon.impactParameter(top::Muon::BS2D);

                ip->setValue(muon->innerTrack()->dxy(beamSpot->position()));

                ip->setError(sqrt(muon->innerTrack()->d0Error() *
                                  muon->innerTrack()->d0Error()
                                  +
                                  .5 *
                                  beamSpot->BeamWidthX() *
                                  beamSpot->BeamWidthX()
                                  +
                                  .5 * beamSpot->BeamWidthY() *
                                  beamSpot->BeamWidthY()));

                topMuon.setInnerValidHits(muon->innerTrack()->numberOfValidHits());
            }

            // Next properties are only defined for the Global Muon
            topMuon.setOuterValidHits(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
            topMuon.setChi2(muon->globalTrack()->chi2());
            topMuon.setNdof(muon->globalTrack()->ndof());

            _topEvent->muons().push_back(topMuon);
        }
    }

    // Process All Jets
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        // Select only energetic jets
        if (30 <= jet->pt() &&
            2.4 >= fabs(jet->eta()))
        {
            top::Jet topJet;

            setP4(topJet.p4(), jet->p4());
            setEnergy(topJet.energy(), jet->caloSpecific());

            topJet.setHits90(jet->jetID().n90Hits);
            topJet.setHpd(jet->jetID().fHPD);

            _topEvent->jets().push_back(topJet);
        }
    }

    // Process All Electrons
    for(ElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        // Select electrons to be stored
        if (15 <= electron->et() &&
            2.5 >= fabs(electron->eta()))
        {
            top::Electron topElectron;

            setP4(topElectron.p4(), electron->p4());

            setIsolation(topElectron.isolation(top::Electron::R03), electron->isolationVariables03());
            setIsolation(topElectron.isolation(top::Electron::R04), electron->isolationVariables04());

            _topEvent->electrons().push_back(topElectron);
        }
    }
    */

    _topTree->Fill();
}

DEFINE_FWK_MODULE(PATTreeMaker);
