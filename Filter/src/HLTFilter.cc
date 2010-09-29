/**
 * HLTFilter
 * top
 *
 * Created by Samvel Khalatian on August 19, 2010
 * Copyright 2010, All rights reserved
 */

#include <algorithm>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "EDModule/Filter/interface/HLTFilter.h"

using top::HLTFilter;

HLTFilter::HLTFilter(const edm::ParameterSet &parameters)
{
    using std::string;

    _hltTag  = parameters.getParameter<string>("hltTag");
    _hltName = parameters.getParameter<string>("hltName");
}

HLTFilter::~HLTFilter()
{
}

bool HLTFilter::filter(edm::Event &event, const edm::EventSetup &)
{
    using std::distance;
    using std::find;

    using edm::Handle;
    using edm::TriggerResults;
    using edm::TriggerNames;

    // Extract Trigger Restuls
    Handle<TriggerResults> triggerResults;
    event.getByLabel(edm::InputTag(_hltTag), triggerResults);

    if (!triggerResults.isValid())
        return false;

    // Get list of triggers
    typedef std::vector<std::string> Triggers;
    const Triggers &triggerNames =
        event.triggerNames(*triggerResults).triggerNames();

    // Find specified trigger in the list
    Triggers::const_iterator trigger = find(triggerNames.begin(),
                                            triggerNames.end(),
                                            _hltName);

    // Return result of (not)found trigger
    return triggerNames.end() != trigger ?
            triggerResults->accept(distance(triggerNames.begin(), trigger)) :
            false;
}

DEFINE_FWK_MODULE(HLTFilter);
