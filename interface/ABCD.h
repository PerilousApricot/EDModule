// -*- C++ -*-
//
// Package:    ABCD
// Class:      ABCD
// 
/**\class ABCD ABCD.cc Exercise/ABCD/src/ABCD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  "Jian Wang"
//        Modified:  Samvel Khalatian
//         Created:  Fri Jun 11 12:14:21 CDT 2010
// $Id$
//
//

#ifndef TOP_EDANALYZERS_ABCD
#define TOP_EDANALYZERS_ABCD

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"

#include "TFile.h"
#include "TTree.h"

class TH1;

//
// class declaration
//

class ABCD : public edm::EDAnalyzer {
    public:
        explicit ABCD(const edm::ParameterSet&);
        ~ABCD();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        TFile *theFile;
        TTree *ftree;

        int _run;
        int _lumi;

        float muon_pt_;
        float muon_eta_;
        float muon_phi_;
        float muon_chi2_;
        int   muon_muonhits_;
        int   muon_trackerhits_;
        float muon_jet_dr_;
        float muon_d0_;
        float _muon_d0pv2d;
        float muon_d0Error_;
        float muon_old_reliso_;
        float _muon_coord[3];

        float _muon_iso03_track;
        float _muon_iso03_ecal;
        float _muon_iso03_hcal;
        float _muon_iso03_ecal_veto;
        float _muon_iso03_hcal_veto;

        float _muon_iso05_track;
        float _muon_iso05_ecal;
        float _muon_iso05_hcal;
        float _muon_iso05_ecal_veto;
        float _muon_iso05_hcal_veto;

        float _muon_mustations;
        int   TrackerMu_;
        int   GlobalMu_;

        float jet_px_[20];
        float jet_py_[20];
        float jet_pz_[20];
        float jet_e_[20];

        size_t _npvs;
        float  _pv_coord[3];
        size_t njets_;
        float met_;
        float w_mt_;

        edm::TriggerNames hltNames_;
        helper::JetIDHelper *jetID;
        edm::InputTag hltTag_;

        TH1 *_cutflow;
        bool _isDataInput;
};

#endif
