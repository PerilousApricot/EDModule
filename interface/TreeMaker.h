/**
 * TreeMaker
 * 
 *
 * Created by Samvel Khalatian on August 20, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_TREE_TREEMAKER
#define TOP_TREE_TREEMAKER

#include <string>

#include "TFile.h"
#include "TTree.h"

#include "Top/Tree/interface/Event.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace reco
{
    namespace helper
    {
        class JetIDHelper;
    }
}

class TreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce Top ROOT Tree
     */
    public:
        TreeMaker(const edm::ParameterSet &);
        virtual ~TreeMaker();

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        std::auto_ptr<top::Event> _topEvent;
        std::auto_ptr<TFile>      _topFile;
        std::auto_ptr<TTree>      _topTree;

        std::string _metTag;
        std::string _muonTag;
        std::string _jetTag;
        std::string _electronTag;

        reco::helper::JetIDHelper *_jetID;
};

#endif
