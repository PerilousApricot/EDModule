/**
 * NtupleMaker
 * 
 *
 * Created by Samvel Khalatian on August 3, 2010
 * Based on by Jian Wang's ABCD code.
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_EDANALYZERS_NTUPLE_MAKER
#define TOP_EDANALYZERS_NTUPLE_MAKER

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace reco
{
    namespace helper
    {
        class JetIDHelper;
    }
}

class TFile;
class TTree;

class NtupleMaker: public edm::EDAnalyzer
{
    public:
        explicit NtupleMaker(const edm::ParameterSet &);
        virtual ~NtupleMaker();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob() ;

        TFile *theFile;
        TTree *ftree;
        float muon_pt_;
        float muon_eta_;
        float muon_phi_;
        float muon_chi2_;
        int muon_muonhits_;
        int muon_trackerhits_;
        float muon_jet_dr_;
        float muon_d0_;
        float muon_d0Error_;
        float muon_old_reliso_;
        int TrackerMu_;
        float jet_pt_[20];
        size_t njets_;
        float met_;
        float w_mt_;
        edm::TriggerNames hltNames_;
        reco::helper::JetIDHelper *jetID;
        edm::InputTag hltTag_;
};

#endif
