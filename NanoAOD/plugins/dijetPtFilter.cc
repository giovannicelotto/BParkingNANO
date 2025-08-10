#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

class DijetPtFilter : public edm::global::EDFilter<> {
public:
  explicit DijetPtFilter(const edm::ParameterSet& iConfig)
    : srcToken_(consumes<edm::View<pat::CompositeCandidate>>(iConfig.getParameter<edm::InputTag>("src"))),
      minPt_(iConfig.getParameter<double>("minPt")) {}

  bool filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    edm::Handle<edm::View<pat::CompositeCandidate>> dijetHandle;
    iEvent.getByToken(srcToken_, dijetHandle);
    
    if (dijetHandle->empty()) return false; // no candidate, reject event

    const auto& dijet = dijetHandle->at(0);
    return (dijet.pt() > minPt_);
  }

private:
  const edm::EDGetTokenT<edm::View<pat::CompositeCandidate>> srcToken_;
  const double minPt_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DijetPtFilter);
