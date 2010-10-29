/**
 * S8TreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TREEMAKER
#define S8_TREEMAKER

#include <memory>
#include <string>

#include "TFile.h"

#include "Tree/System8/interface/S8Event.h"
#include "Tree/System8/interface/S8TreeInfo.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

class TTree;

namespace reco
{
    class Vertex;
}

class S8TreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce S8 ROOT Tree
     */
    public:
        S8TreeMaker(const edm::ParameterSet &);
        virtual ~S8TreeMaker();

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        void processEventID(const edm::Event &);
        void processElectrons(const edm::Event &);
        void processJets(const edm::Event &);
        void processMuons(const edm::Event &);
        void processPrimaryVertices(const edm::Event &);
        void processTriggers(const edm::Event &);

        bool isGoodPrimaryVertex(const reco::Vertex &, const bool & = false); 

        std::auto_ptr<s8::Event>     _event;
        std::auto_ptr<s8::TreeInfo>  _treeInfo;
        TTree                       *_tree;
        PFJetIDSelectionFunctor      _jetSelector;

        std::string _primaryVertices;
        std::string _jets;
        std::string _muons;
        std::string _electrons;
        std::string _triggers;
        bool        _isPythia;
};

#endif
