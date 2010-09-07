/**
 * Tools
 * top
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_EDANALYZER_TOOLS
#define TOP_EDANALYZER_TOOLS

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class TLorentzVector;

namespace reco
{
    class MuonIsolation;
}

namespace top
{
    class JetEnergy;
    class ElectronIsolation;
    class MuonIsolation;

    namespace tools
    {
        void setP4(TLorentzVector *,
                   const math::XYZTLorentzVector *);

        inline void setP4(TLorentzVector *p4_1,
                          const math::XYZTLorentzVector &p4_2)
        {
            setP4(p4_1, &p4_2);
        }



        void setEnergy(top::JetEnergy *,
                        const reco::CaloJet::Specific *);

        inline void setEnergy(top::JetEnergy *eng_1,
                        const reco::CaloJet::Specific &eng_2)
        {
            setEnergy(eng_1, &eng_2);
        }



        void setIsolation(top::MuonIsolation *,
                          const reco::MuonIsolation *);


        inline void setIsolation(top::MuonIsolation *topIso,
                                 const reco::MuonIsolation &recoIso)
        {
            setIsolation(topIso, &recoIso);
        }



        void setIsolation(top::ElectronIsolation *,
                          const reco::GsfElectron::IsolationVariables *);

        inline void setIsolation(top::ElectronIsolation *topIso,
                                 const reco::GsfElectron::IsolationVariables &recoIso)
        {
            setIsolation(topIso, &recoIso);
        }
    }
}

#endif
