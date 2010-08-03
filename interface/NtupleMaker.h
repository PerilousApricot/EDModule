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

#include <string>

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

class TFile;
class TTree;

class NtupleMaker : public edm::EDAnalyzer
{
    /*
     * Produce ttmuj Ntuple
     */
    public:
        NtupleMaker(const edm::ParameterSet &);
        virtual ~NtupleMaker();

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        bool processTrigger(const edm::Event &);
        bool processPrimaryVertex(const edm::Event &);

        TFile *_ntuple;
        TTree *_tree;

        std::string _electronCollection;
        std::string _muonCollection;
        std::string _jetCollection;
        std::string _metCollection;
        std::string _pvCollection;

        struct Muon
        {
            float pt;
            float eta;
            float phi;
            float d0;
            float d0err;
            float oldRelIso;
            float chi2;
            int   muonHits;
            int   trackerHits;
            float deltaR;
            int   isTracker;
            int   isGlobal;
        };

        Muon   _muon;
        size_t _njets;
        float  _jet_pt[20];
        float  _met;
        float  _w_mt;
};

#endif
