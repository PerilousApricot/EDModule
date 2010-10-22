/**
 * S8TreeMaker
 * 
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "Tree/System8/interface/S8EventID.h"
#include "Tree/System8/interface/S8GenEvent.h"
#include "Tree/System8/interface/S8GenParticle.h"
#include "Tree/System8/interface/S8Jet.h"
#include "Tree/System8/interface/S8Lepton.h"
#include "EDModule/Analyzer/interface/Tools.h"

#include "EDModule/Analyzer/interface/S8TreeMaker.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using edm::errors::ErrorCodes;
using edm::Handle;
using edm::InputTag;
using edm::LogWarning;
using edm::LogInfo;
using edm::ParameterSet;
using edm::TriggerResults;
using edm::TriggerNames;

using reco::Vertex;

using namespace top::tools;

S8TreeMaker::S8TreeMaker(const edm::ParameterSet &config)
{
    _primaryVertices = config.getParameter<string>("primaryVertices");
    _jets = config.getParameter<string>("jets");
    _muons = config.getParameter<string>("muons");
    _electrons = config.getParameter<string>("electrons");
    _triggers = config.getParameter<string>("triggers");

    _isPythia = config.getParameter<bool>("isPythia");
}

S8TreeMaker::~S8TreeMaker()
{
}

void S8TreeMaker::beginJob()
{
    edm::Service<TFileService> fileService;

    _tree = fileService->make<TTree>("s8", "System8 tree.");

    _event.reset(new s8::Event());
    _tree->Branch("event", _event.get(), 32000, 0);
}

void S8TreeMaker::endJob()
{
    if (!_event.get())
        return;

    // Note: Event should be destroyed after ROOT file is written and closed.
    //
    _event.reset();
}

void S8TreeMaker::analyze(const edm::Event &event, const edm::EventSetup &)
{
    if (!_event.get())
        return;

    using pat::JetCollection;
    using pat::MuonCollection;
    using pat::ElectronCollection;

    // Extract Muons
    //
    Handle<MuonCollection> muons;
    event.getByLabel(InputTag(_muons), muons);

    // Extract Electrons
    //
    Handle<ElectronCollection> electrons;
    event.getByLabel(InputTag(_electrons), electrons);

    if (!muons.isValid() && !electrons.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Leptons.";

        return;
    }

    // Extract Jets
    //
    Handle<JetCollection> jets;
    event.getByLabel(InputTag(_jets), jets);

    if (!jets.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Jets.";

        return;
    }

    // Primary Vertices
    //
    typedef vector<Vertex> PVCollection;

    Handle<PVCollection> primaryVertices;
    event.getByLabel(InputTag(_primaryVertices), primaryVertices);

    if (!primaryVertices.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Primary Vertices.";

        return;
    }

    if (primaryVertices->empty())
    {
        LogWarning("S8TreeMaker")
            << "primary vertices collection is empty.";

        return;
    }

    // Triggers
    //
    Handle<TriggerResults> triggers;
    event.getByLabel(InputTag(_triggers), triggers);

    if (!triggers.isValid())
    {
        LogWarning("S8TreeMaker")
            << "failed to extract Triggers";

        return;
    }

    _event->reset();

    // check if Event is Pythia
    //
    if (_isPythia && !event.isRealData())
    {
        Handle<GenEventInfoProduct> generator;
        event.getByLabel(InputTag("generator"), generator);

        if (!generator.isValid())
        {
            LogWarning("S8TreeMaker")
                << "failed to extract Generator";

            return;
        }

        _event->gen().setPtHat(generator->qScale());
    }

    {
        s8::EventID &id = _event->id();

        id.setRun(event.id().run());
        id.setLumiBlock(event.id().luminosityBlock());
        id.setEvent(event.id().event());
    }

    processTriggers(event, *triggers);

    // Process Primary Vertices
    //
    for(PVCollection::const_iterator vertex = primaryVertices->begin();
        primaryVertices->end() != vertex;
        ++vertex)
    {
        if (!isGoodPrimaryVertex(*vertex, event.isRealData()))
            continue;

        s8::PrimaryVertex s8Vertex;

        setVertex(s8Vertex.vertex(), vertex->position());
        s8Vertex.setNdof(vertex->ndof());
        s8Vertex.setRho(vertex->position().Rho());

        _event->primaryVertices().push_back(s8Vertex);
    }

    // Process all Muons
    //
    for(MuonCollection::const_iterator muon = muons->begin();
        muons->end() != muon;
        ++muon)
    {
        // Only Muons with basic cuts are saved:
        //
        if (1 >= muon->numberOfMatches())
            continue;

        s8::Lepton s8Muon;

        setP4(s8Muon.p4(), muon->p4());
        setVertex(s8Muon.vertex(), muon->vertex());

        s8Muon.impactParameter().first = muon->dB();
        s8Muon.impactParameter().second = muon->edB();

        // Extract GenParticle information
        //
        if (muon->genLepton())
        {
            s8::GenParticle &s8GenParticle = s8Muon.genParticle();

            setP4(s8GenParticle.p4(), muon->genLepton()->p4());
            setVertex(s8GenParticle.vertex(), muon->genLepton()->vertex());

            s8GenParticle.setId(muon->genLepton()->pdgId());
            if (muon->genLepton()->mother())
                s8GenParticle.setParentId(muon->genLepton()->mother()->pdgId());
        }

        _event->muons().push_back(s8Muon);
    } // End loop over muons

    // Process all Electrons
    //
    for(ElectronCollection::const_iterator electron = electrons->begin();
        electrons->end() != electron;
        ++electron)
    {
        s8::Lepton s8Electron;

        setP4(s8Electron.p4(), electron->p4());
        setVertex(s8Electron.vertex(), electron->vertex());

        s8Electron.impactParameter().first = electron->dB();
        s8Electron.impactParameter().second = electron->edB();

        // Extract GenParticle information
        //
        if (electron->genLepton())
        {
            s8::GenParticle &s8GenParticle = s8Electron.genParticle();

            setP4(s8GenParticle.p4(), electron->genLepton()->p4());
            setVertex(s8GenParticle.vertex(), electron->genLepton()->vertex());

            s8GenParticle.setId(electron->genLepton()->pdgId());
            if (electron->genLepton()->mother())
                s8GenParticle.setParentId(electron->genLepton()->mother()->pdgId());
        }

        _event->electrons().push_back(s8Electron);
    } // End loop over electrons

    // Process All Jets
    //
    for(JetCollection::const_iterator jet = jets->begin();
        jets->end() != jet;
        ++jet)
    {
        using s8::Jet;
        Jet s8Jet;

        setP4(s8Jet.p4(), jet->p4());

        s8Jet.setFlavour(jet->partonFlavour());
        s8Jet.setTracks(jet->associatedTracks().size());

        // Save b-taggers
        //
        s8Jet.setBTag(Jet::TCHE,
                      jet->bDiscriminator("trackCountingHighEffBJetTags"));

        s8Jet.setBTag(Jet::TCHP,
                      jet->bDiscriminator("trackCountingHighPurBJetTags"));

        s8Jet.setBTag(Jet::JP,
                      jet->bDiscriminator("jetProbabilityBJetTags"));

        s8Jet.setBTag(Jet::SSV,
                      jet->bDiscriminator("simpleSecondaryVertexBJetTags"));

        s8Jet.setBTag(Jet::SSVHE,
                      jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));

        s8Jet.setBTag(Jet::SSVHP,
                      jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));

        _event->jets().push_back(s8Jet);
    }

    _tree->Fill();
}

bool S8TreeMaker::isGoodPrimaryVertex(const Vertex &vertex,
                                    const bool &isData)
{
    return !vertex.isFake() &&
            4 <= vertex.ndof() &&
            (isData ? 24 : 15) >= fabs(vertex.z()) &&
            2 >= fabs(vertex.position().Rho());
}

void S8TreeMaker::processTriggers(const edm::Event &event,
                                  const TriggerResults &triggers)
{
    // Get list of triggers
    //
    typedef std::vector<std::string> Triggers;
    const Triggers &triggerNames =
        event.triggerNames(triggers).triggerNames();

    for(Triggers::const_iterator trigger = triggerNames.begin();
        triggerNames.end() != trigger;
        ++trigger)
    {
        // check if this is any trigger version
        //
        s8::Trigger s8Trigger;
        
        smatch matches;
        if (!regex_match(*trigger, matches, regex("^(.*)(?:_v(\\d+))?$")))
        {
            LogWarning("S8TreeMaker")
                << "Do not understand trigger: " << *trigger;
        }

        s8Trigger.setName(matches[1]);

        if (matches[2].matched)
            s8Trigger.setVersion(lexical_cast<int>(matches[2]));

        s8Trigger.setIsPass(triggers.accept(distance(triggerNames.begin(),
                                                     trigger)));

        _event->triggers().push_back(s8Trigger);
    }
}

DEFINE_FWK_MODULE(S8TreeMaker);
