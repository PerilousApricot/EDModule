/**
 * NtupleMaker
 * 
 *
 * Created by Samvel Khalatian on August 3, 2010
 * Based on by Jian Wang's ABCD code.
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Top/EDAnalyzers/interface/NtupleMaker.h"

double relIso(const pat::Muon &muon)
{
    return (muon.trackIso() + muon.ecalIso() + muon.hcalIso()) / muon.pt();
}

NtupleMaker::NtupleMaker(const edm::ParameterSet &iConfig)
{
}

NtupleMaker::~NtupleMaker()
{
}

// ------------ method called to for each event  ------------
void NtupleMaker::analyze(const edm::Event& event,
                          const edm::EventSetup& iSetup)
{
    using std::vector;
    using std::clog;
    using std::vector;
    using pat::ElectronCollection;
    using pat::JetCollection;
    using pat::MuonCollection;
    using pat::METCollection;

    // -[ Trigger ]-
    if (!processTrigger(event))
        return;

    // -[ Primary Vertex ]-
    if (!processPrimaryVertex(event))
        return;

    // -[ Electrons ]-
    edm::Handle<ElectronCollection> electrons;
    event.getByLabel(edm::InputTag(_electronCollection), electrons);

    if (!electrons.isValid())
    {
        std::clog << "[NtupleMaker] "
            << "failed to extract electrons."
            << std::endl;

        return;
    }

    // -[ Muons ]-
    edm::Handle<MuonCollection> muons;
    event.getByLabel(edm::InputTag(_muonCollection), muons);

    if (!muons.isValid())
    {
        clog << "[NtupleMaker] "
            << "failed to extract muons." << std::endl;

        return;
    }

    // -[ Jets ]-
    edm::Handle<JetCollection> jets;
    event.getByLabel(edm::InputTag(_jetCollection), jets);

    if (!jets.isValid())
    {
        clog << "[NtupleMaker] "
            << "failed to extract jets." << std::endl;

        return;
    }

    // -[ MET ]-
    edm::Handle<METCollection> mets;
    event.getByLabel(edm::InputTag("patMETs"), mets);

    if (!mets.isValid())
    {
        clog << "[WAnalyzer] "
            << "failed to extract patMETs." << std::endl;

        return;
    }

    // Find Good jets
    typedef vector<const pat::Jet *>  Jets;
    typedef vector<const pat::Muon *> Muons;

    Jets goodJets;
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        if (30 < jet->pt()
            && 2.4 > fabs(jet->eta())
            && .01 < jet->emEnergyFraction()
            && 1 < jet->jetID().n90Hits
            && .98 > jet->jetID().fHPD)

            goodJets.push_back(&(*jet));
    }

    // -[ Process Muons ]-
    Muons looseMuons;
    Muons tightMuons;
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // Loose Cuts
        if (muon->isGlobalMuon()
            && 10 < muon->pt()
            && 2.5 > fabs(muon->eta())
            && .2 > relIso(*muon))
        {
            looseMuons.push_back(&(*muon));

            // Tight Cuts
            if (muon->isTrackerMuon()
                && 20 < muon->pt()
                && 2.1 > fabs(muon->eta())
                && .05 > relIso(*muon)
                && 10 > muon->globalTrack()->normalizedChi2()
                && 0 < muon->globalTrack()->hitPattern().numberOfValidMuonHits()
                && 10 < muon->innerTrack()->numberOfValidHits()
                && .02 > muon->dB())
            {
                const pat::Muon *tightMuon = &(*muon);

                // We've got a candidate for the Tight muon. Check if it is
                // isolated with respect to all goodJets
                for(Jets::const_iterator jet = goodJets.begin();
                    0 != tightMuon && goodJets.end() != jet;
                    ++jet)
                {
                    if (.3 >= deltaR(muon->eta(), muon->phi(),
                                     (*jet)->eta(), (*jet)->phi()))

                        tightMuon = 0;
                }

                if (tightMuon)
                    tightMuons.push_back(tightMuon);
            }
        }
    } // end loop over muons

/*
    size_t n_loose=0;
    size_t n_tight=0;
    size_t n_muon=0;

    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        double reliso = (muon->isolationR03().hadEt +
                         muon->isolationR03().emEt +
                         muon->isolationR03().sumPt) / muon->pt();


        // find Minimum DeltaR between muon and any good jet
        double DeltaR = 3.;
        for(Collection::const_iterator jet = jets->begin();
            jets->end() != jet;
            ++jet)
        {
            if (30 < jet->pt()
                && 2.4 > fabs(jet->eta())
                && .01 < jet->emEnergyFraction()
                && 1 < jet->jetID().n90Hits
                && .98 > jet->jetID().fHPD)
            {
                double dr = deltaR(muon->eta(),  muon->phi(),
                                   jet->eta(), jet->phi());
                if(dr < DeltaR)
                    DeltaR = dr;
            }
        }


        // loose muon ID
        if (muon->isGlobalMuon() &&
            abs(muon->eta()) < 2.5 &&
            muon->pt() > 10.)
        {
            if (reliso < 0.2)
                ++n_loose;

            ++n_muon;
            if (n_tight == 0)
            {
                muon_d0_ = -1.* muon->innerTrack()->dxy(point);
                muon_d0Error_ = sqrt(muon->innerTrack()->d0Error() *
                                     muon->innerTrack()->d0Error() +
                                     beamSpot.BeamWidthX() *
                                     beamSpot.BeamWidthX());

                muon_old_reliso_= muon->pt() / (muon->pt() +
                                              muon->isolationR03().sumPt +
                                              muon->isolationR03().emEt +
                                              muon->isolationR03().hadEt);

                muon_pt_ = muon->pt();
                muon_eta_ = muon->eta();
                muon_phi_ = muon->phi();
                muon_chi2_ = muon->globalTrack()->normalizedChi2();
                muon_muonhits_ = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
                muon_trackerhits_ = muon->innerTrack()->numberOfValidHits();
                TrackerMu_ = muon->isTrackerMuon();
                muon_jet_dr_ = DeltaR;
                double w_et = cmet->et() + muon->pt();
                double w_px = cmet->px() + muon->px();
                double w_py = cmet->py() + muon->py();
                w_mt_ = sqrt(w_et * w_et - w_px * w_px - w_py * w_py);
            }
        }

        
        // tight muon cuts
        if (muon->isGlobalMuon() &&
            muon->isTrackerMuon() &&
            abs(muon->eta()) < 2.1 &&
            muon->pt() > 20. &&
            reliso < 0.05 &&
            muon->innerTrack()->numberOfValidHits() >= 11 &&
            muon->globalTrack()->normalizedChi2() < 10. &&
            muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
            DeltaR > 0.3 &&
            abs(muon->innerTrack()->dxy(point)) < 0.02)
        {
            ++n_tight;
            muon_d0_ = -1. * muon->innerTrack()->dxy(point);
            muon_d0Error_ = sqrt(muon->innerTrack()->d0Error() *
                                 muon->innerTrack()->d0Error() +
                                 beamSpot.BeamWidthX() *
                                 beamSpot.BeamWidthX());

            muon_old_reliso_= muon->pt() / (muon->pt() +
                                          muon->isolationR03().sumPt +
                                          muon->isolationR03().emEt +
                                          muon->isolationR03().hadEt);

            muon_pt_ = muon->pt();
            muon_eta_ = muon->eta();
            muon_phi_ = muon->phi();
            muon_chi2_ = muon->globalTrack()->normalizedChi2();
            muon_muonhits_ = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
            muon_trackerhits_ = muon->innerTrack()->numberOfValidHits();
            TrackerMu_ = muon->isTrackerMuon();
            muon_jet_dr_ = DeltaR;
            double w_et = cmet->et() + muon->pt();
            double w_px = cmet->px() + muon->px();
            double w_py = cmet->py() + muon->py();
            w_mt_ = sqrt(w_et * w_et - w_px * w_px - w_py * w_py);
        }
    }

    if (!(n_muon == 1 ||
          (n_tight == 1 && n_loose == 1)))
        return;
        */

    if (1 != tightMuons.size() ||
        1 != looseMuons.size())

        return;

    // Electron veto
    for(ElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        if (15 < electron->et()
            && 2.5 > fabs(electron->eta())
            && .2 > (electron->dr03TkSumPt()
                        + electron->dr03EcalRecHitSumEt()
                        + electron->dr03HcalTowerSumEt())
                        / electron->et())
            return;
    }

    const pat::Muon *tightMuon = tightMuons[0];

    // Store Jets Pt
    double minDeltaR = 1000000;
    for(Jets::const_iterator jet = goodJets.begin();
        goodJets.end() != jet;
        ++jet)
    {
        _jet_pt[jet - goodJets.begin()] = (*jet)->pt();

        double tmpDeltaR = deltaR(tightMuon->eta(), tightMuon->phi(),
                                  (*jet)->eta(), (*jet)->phi());

        if (tmpDeltaR < minDeltaR)
            minDeltaR = tmpDeltaR;
    }

    METCollection::const_iterator met = mets->begin();

    // Store Tight Muon values
    _muon.pt = tightMuon->pt();
    _muon.eta = tightMuon->eta();
    _muon.phi = tightMuon->phi();
    _muon.d0 = tightMuon->dB();
    //_muon.d0Error_ = tightMuon->dBerr();
    _muon.oldRelIso = tightMuon->pt() / (tightMuon->pt() +
                                         tightMuon->trackIso() +
                                         tightMuon->ecalIso() +
                                         tightMuon->hcalIso());

    _muon.chi2 = tightMuon->globalTrack()->normalizedChi2();
    _muon.muonHits = tightMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
    _muon.trackerHits = tightMuon->innerTrack()->numberOfValidHits();
    _muon.deltaR = minDeltaR;
    _muon.isTracker = tightMuon->isTrackerMuon();
    _muon.isGlobal = tightMuon->isGlobalMuon();

    // Calculate momentum 4-vector: Met + Muon
    _njets = goodJets.size();
    _met = met->et();
    _w_mt = (tightMuon->p4() + met->p4()).mass();

    _tree->Fill();
}

void NtupleMaker::beginJob()
{
    _ntuple = new TFile("ttmuj_ntuple.root", "RECREATE");
    _tree = new TTree("top", "top");

    _tree->Branch("muon_pt",          &(_muon.pt),          "muon_pt/F");
    _tree->Branch("muon_eta",         &(_muon.eta),         "muon_eta/F");
    _tree->Branch("muon_phi",         &(_muon.phi),         "muon_phi/F"); 
    _tree->Branch("muon_d0",          &(_muon.d0),          "muon_d0/F");
    _tree->Branch("muon_d0Error",     &(_muon.d0err),       "muon_d0Error/F");
    _tree->Branch("muon_old_reliso",  &(_muon.oldRelIso),   "muon_old_reliso/F");
    _tree->Branch("muon_chi2",        &(_muon.chi2),        "muon_chi2/F");
    _tree->Branch("muon_muonhits",    &(_muon.muonHits),    "muon_muonhits/I");
    _tree->Branch("muon_trackerhits", &(_muon.trackerHits), "muon_trackerhits/I");
    _tree->Branch("muon_jet_dr",      &(_muon.deltaR),      "muon_jet_dr/F");
    _tree->Branch("TrackerMu",        &(_muon.isTracker),   "TrackerMu/I");
    _tree->Branch("njets",            &(_njets), "njets/I");
    _tree->Branch("jet_pt",           &(_jet_pt), "jet_pt[njets]/F");
    _tree->Branch("met",              &(_met), "met/F");
    _tree->Branch("w_mt",             &(_w_mt), "w_mt/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob()
{
    _ntuple->cd();
    _tree->Write();

    delete _tree;
    delete _ntuple;
}

bool NtupleMaker::processTrigger(const edm::Event &event)
{
    edm::Handle<pat::TriggerEvent> trigger;
    event.getByLabel(edm::InputTag("patTriggerEvent"), trigger);

    if (!trigger.isValid())
    {
        std::clog << "[NtupleMaker] "
            << "failed to extract trigger results. Skip event."
            << std::endl;

        return false;
    }

    // Check if trigger was run and at least one accepted.
    if (!trigger->wasRun() || !trigger->wasAccept())
        return false;

    // Extract Muon Path
    const pat::TriggerPath *path = trigger->path("HLT_Mu9");

    return path && path->wasAccept();
}

bool NtupleMaker::processPrimaryVertex(const edm::Event &event)
{
    using reco::VertexCollection;

    edm::Handle<VertexCollection> primaryVertex;
    event.getByLabel(edm::InputTag("offlinePrimaryVertices"),
                     primaryVertex);

    if (!primaryVertex.isValid())
    {
        std::clog << "[NtupleMaker] "
            << "failed to extract primary Vertecies."
            << std::endl;

        return false;
    }

    if (!primaryVertex->size())
        return false;

    VertexCollection::const_iterator firstVertex = primaryVertex->begin();

    // Validate Primary vertex
    return !firstVertex->isFake()
           && 4 < firstVertex->ndof()
           && 15 > fabs(firstVertex->z())
           && 2.0 > fabs(firstVertex->position().Rho());
}

DEFINE_FWK_MODULE(NtupleMaker);
