/**
 * Tools
 * top
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#include "TLorentzVector.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "Top/Tree/interface/Jet.h"
#include "Top/Tree/interface/JetEnergy.h"
#include "Top/Tree/interface/ElectronIsolation.h"
#include "Top/Tree/interface/MuonIsolation.h"

#include "Top/EDAnalyzers/interface/Tools.h"

void top::tools::setP4(TLorentzVector *topP4,
                       const math::XYZTLorentzVector *cmsswP4)
{
    topP4->SetPxPyPzE(cmsswP4->px(),
                      cmsswP4->py(),
                      cmsswP4->pz(),
                      cmsswP4->energy());
}

void top::tools::setEnergy(top::Jet &jet,
                           const reco::CaloJet::Specific &specific)
{
    top::JetEnergy energy;

    energy.setEcalMax(specific.mMaxEInEmTowers);
    energy.setHcalMax(specific.mMaxEInHadTowers);

    energy.setHcalInHO(specific.mHadEnergyInHO);
    energy.setHcalInHB(specific.mHadEnergyInHB);
    energy.setHcalInHF(specific.mHadEnergyInHF);
    energy.setHcalInHE(specific.mHadEnergyInHE);

    energy.setEcalInEB(specific.mEmEnergyInEB);
    energy.setEcalInEE(specific.mEmEnergyInEE);
    energy.setEcalInHF(specific.mEmEnergyInHF);

    energy.setEcalFraction(specific.mEnergyFractionEm);
    energy.setHcalFraction(specific.mEnergyFractionHadronic);

    jet.setEnergy(energy);
}

void top::tools::setIsolation(top::Muon &muon,
                              const top::Muon::ISO &iso,
                              const reco::MuonIsolation &recoIso)
{
    top::MuonIsolation topIso;

    topIso.setTrackPt(recoIso.sumPt);
    topIso.setEcalEt(recoIso.emEt);
    topIso.setHcalEt(recoIso.hadEt);

    topIso.setTracks(recoIso.nTracks);
    topIso.setJets(recoIso.nJets);

    topIso.setTrackPtVeto(recoIso.trackerVetoPt);
    topIso.setEcalEtVeto(recoIso.emVetoEt);
    topIso.setHcalEtVeto(recoIso.hadVetoEt);

    muon.setIsolation(iso, topIso);
}

void top::tools::setIsolation(top::Electron &electron,
                              const top::Electron::ISO &iso,
                              const reco::GsfElectron::IsolationVariables &recoIso)
{
    top::ElectronIsolation topIso;

    topIso.setTrackPt(recoIso.tkSumPt);
    topIso.setEcalEt(recoIso.ecalRecHitSumEt);
    topIso.setHcalEt(recoIso.hcalDepth1TowerSumEt + recoIso.hcalDepth2TowerSumEt);

    electron.setIsolation(iso, topIso);
}
