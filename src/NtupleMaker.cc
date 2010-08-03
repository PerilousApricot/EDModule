/**
 * NtupleMaker
 * 
 *
 * Created by Samvel Khalatian on August 3, 2010
 * Based on by Jian Wang's ABCD code.
 * Copyright 2010, All rights reserved
 */

#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"

#include "Top/EDAnalyzers/interface/NtupleMaker.h"

using std::vector;
using std::string;

using namespace edm;
using namespace reco;
using reco::helper::JetIDHelper;

NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig):
  hltTag_(iConfig.getParameter<InputTag>("hltTag"))
{
   //now do what ever initialization is needed
    jetID = new JetIDHelper(iConfig.getParameter<ParameterSet>("JetIDParams"));
}

NtupleMaker::~NtupleMaker()
{
    delete jetID;
}

// ------------ method called to for each event  ------------
void NtupleMaker::analyze(const edm::Event& iEvent,
                          const edm::EventSetup& iSetup)
{
    Handle<TriggerResults> hlt;
    iEvent.getByLabel(hltTag_,hlt);
    const TriggerNames & hltNames_ = iEvent.triggerNames(*hlt);

    Handle<MuonCollection> muons;
    iEvent.getByLabel("muons",muons);

    Handle<CaloJetCollection> jets;
    iEvent.getByLabel("ak5CaloJetsL2L3",jets);

    Handle<vector<reco::Vertex> > pvtx;
    iEvent.getByLabel("offlinePrimaryVertices",pvtx);

    BeamSpot beamSpot;
    Handle<BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    beamSpot = *beamSpotHandle;
    math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());

    Handle<GsfElectronCollection> electrons;
    iEvent.getByLabel("gsfElectrons",electrons);

    const CaloMET *cmet;
    Handle<CaloMETCollection> cmetHandle;
    iEvent.getByLabel("metMuonJESCorAK5", cmetHandle);
    const CaloMETCollection *cmetCol = cmetHandle.product();
    cmet = &(cmetCol->front());

    bool hlt_mu_ = false; 
    for(size_t itrig = 0; itrig != hlt->size(); ++itrig)
    {
        string hltName = hltNames_.triggerName(itrig);

        if (hltName == "HLT_Mu9" &&
            hlt->accept(itrig))
        {
            hlt_mu_ = true;
            break;
        }
    }

    if (!hlt_mu_)
        return;

    if (pvtx->size() < 1)
        return;

    Vertex const & pv = pvtx->at(0);
    if (!(!pv.isFake() &&
            pv.ndof() > 4 &&
            fabs(pv.z()) < 15. &&
            pv.position().Rho() < 2.0))
        return;

    size_t n_loose=0;
    size_t n_tight=0;
    size_t n_muon=0;
    size_t n_electron=0;

    for(MuonCollection::const_iterator mu = muons->begin();
        mu != muons->end();
        ++mu)
    {
        double reliso = (mu->isolationR03().hadEt +
                         mu->isolationR03().emEt +
                         mu->isolationR03().sumPt) / mu->pt();

        double DeltaR = 3.;
        for(CaloJetCollection::const_iterator jet = jets->begin();
            jet != jets->end();
            ++jet)
        {
            jetID->calculate(iEvent, *jet);
            if (jet->pt() > 30. &&
                abs(jet->eta()) < 2.4 &&
                jet->emEnergyFraction() > 0.01 &&
                jetID->n90Hits() > 1 &&
                jetID->fHPD() < 0.98)
            {
                double dr = deltaR(mu->eta(),  mu->phi(),
                                   jet->eta(), jet->phi());
                if(dr < DeltaR)
                    DeltaR=dr;
            }
        }


        if (mu->isGlobalMuon() &&
            abs(mu->eta()) < 2.5 &&
            mu->pt() > 10.)
        {
            if (reliso < 0.2)
                ++n_loose;

            ++n_muon;
            if (n_tight == 0)
            {
                muon_d0_ = -1.* mu->innerTrack()->dxy(point);
                muon_d0Error_ = sqrt(mu->innerTrack()->d0Error() *
                                     mu->innerTrack()->d0Error() +
                                     beamSpot.BeamWidthX() *
                                     beamSpot.BeamWidthX());

                muon_old_reliso_= mu->pt() / (mu->pt() +
                                              mu->isolationR03().sumPt +
                                              mu->isolationR03().emEt +
                                              mu->isolationR03().hadEt);

                muon_pt_ = mu->pt();
                muon_eta_ = mu->eta();
                muon_phi_ = mu->phi();
                muon_chi2_ = mu->globalTrack()->normalizedChi2();
                muon_muonhits_ = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
                muon_trackerhits_ = mu->innerTrack()->numberOfValidHits();
                TrackerMu_ = mu->isTrackerMuon();
                muon_jet_dr_ = DeltaR;
                double w_et = cmet->et() + mu->pt();
                double w_px = cmet->px() + mu->px();
                double w_py = cmet->py() + mu->py();
                w_mt_ = sqrt(w_et * w_et - w_px * w_px - w_py * w_py);
            }
        }

        if (mu->isGlobalMuon() &&
            mu->isTrackerMuon() &&
            abs(mu->eta()) < 2.1 &&
            mu->pt() > 20. &&
            reliso < 0.05 &&
            mu->innerTrack()->numberOfValidHits() >= 11 &&
            mu->globalTrack()->normalizedChi2() < 10. &&
            mu->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
            DeltaR > 0.3 &&
            abs(mu->innerTrack()->dxy(point)) < 0.02)
        {
            ++n_tight;
            muon_d0_ = -1. * mu->innerTrack()->dxy(point);
            muon_d0Error_ = sqrt(mu->innerTrack()->d0Error() *
                                 mu->innerTrack()->d0Error() +
                                 beamSpot.BeamWidthX() *
                                 beamSpot.BeamWidthX());

            muon_old_reliso_= mu->pt() / (mu->pt() +
                                          mu->isolationR03().sumPt +
                                          mu->isolationR03().emEt +
                                          mu->isolationR03().hadEt);

            muon_pt_ = mu->pt();
            muon_eta_ = mu->eta();
            muon_phi_ = mu->phi();
            muon_chi2_ = mu->globalTrack()->normalizedChi2();
            muon_muonhits_ = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
            muon_trackerhits_ = mu->innerTrack()->numberOfValidHits();
            TrackerMu_ = mu->isTrackerMuon();
            muon_jet_dr_ = DeltaR;
            double w_et = cmet->et() + mu->pt();
            double w_px = cmet->px() + mu->px();
            double w_py = cmet->py() + mu->py();
            w_mt_ = sqrt(w_et * w_et - w_px * w_px - w_py * w_py);
        }
    }

    if (!(n_muon == 1 ||
          (n_tight == 1 && n_loose == 1)))
        return;

    for(GsfElectronCollection::const_iterator elec = electrons->begin();
        elec != electrons->end();
        ++elec)
    {
        double RelIso = (elec->dr03TkSumPt()+elec->dr03EcalRecHitSumEt()+elec->dr03HcalTowerSumEt())/elec->et();
        if (elec->et() > 15. &&
            abs(elec->eta()) < 2.5 &&
            RelIso < 0.2)

            ++n_electron;
    }

    if (n_electron != 0)
        return;

    njets_ = 0;
    for(CaloJetCollection::const_iterator jet = jets->begin();
        jet != jets->end(); 
        ++jet)
    {
        jetID->calculate(iEvent, *jet);
        if (jet->pt() > 30. &&
            abs(jet->eta()) < 2.4 &&
            jet->emEnergyFraction() > 0.01 &&
            jetID->n90Hits() > 1 &&
            jetID->fHPD() < 0.98)
        {
            jet_pt_[njets_] = jet->pt();
            ++njets_;
        }
    }

    met_ = cmet->et();

    ftree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void NtupleMaker::beginJob()
{
    theFile = new TFile("ttmuj_ntuple.root", "RECREATE");
    ftree = new TTree("top", "top");

    ftree->Branch("muon_pt", &muon_pt_, "muon_pt/F");
    ftree->Branch("muon_eta", &muon_eta_, "muon_eta/F");
    ftree->Branch("muon_phi", &muon_phi_, "muon_phi/F"); 
    ftree->Branch("muon_d0", &muon_d0_, "muon_d0/F");
    ftree->Branch("muon_d0Error", &muon_d0Error_, "muon_d0Error/F");
    ftree->Branch("muon_old_reliso", &muon_old_reliso_, "muon_old_reliso/F");
    ftree->Branch("muon_chi2", &muon_chi2_, "muon_chi2/F");
    ftree->Branch("muon_muonhits", &muon_muonhits_, "muon_muonhits/I");
    ftree->Branch("muon_trackerhits", &muon_trackerhits_, "muon_trackerhits/I");
    ftree->Branch("muon_jet_dr", &muon_jet_dr_, "muon_jet_dr/F");
    ftree->Branch("TrackerMu", &TrackerMu_, "TrackerMu/I");
    ftree->Branch("njets", &njets_, "njets/I");
    ftree->Branch("jet_pt", jet_pt_, "jet_pt[njets]/F");
    ftree->Branch("met", &met_, "met/F");
    ftree->Branch("w_mt", &w_mt_, "w_mt/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob()
{
    theFile->cd();
    ftree->Write();

    delete ftree;
    delete theFile;
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
