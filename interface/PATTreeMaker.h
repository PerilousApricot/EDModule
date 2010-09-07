/**
 * PATTreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_TREE_PATTREEMAKER
#define TOP_TREE_PATTREEMAKER

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class TTree;

namespace top
{
    class Event;
}

class PATTreeMaker : public edm::EDAnalyzer
{
    /*
     * Produce Top ROOT Tree
     */
    public:
        PATTreeMaker(const edm::ParameterSet &);
        virtual ~PATTreeMaker();

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        top::Event *_topEvent;
        TTree      *_topTree;

        std::string _metTag;
        std::string _muonTag;
        std::string _jetTag;
        std::string _electronTag;
        std::string _beamSpotTag;
};

#endif
