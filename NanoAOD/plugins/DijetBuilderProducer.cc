#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/View.h"

#include "TLorentzVector.h"


class DijetBuilderProducer : public edm::global::EDProducer<> {
public:
  explicit DijetBuilderProducer(const edm::ParameterSet&);
  ~DijetBuilderProducer() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<edm::View<pat::Jet>> jets_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  edm::EDGetTokenT<edm::ValueMap<float>> bRegCorr_;
};

DijetBuilderProducer::DijetBuilderProducer(const edm::ParameterSet& iConfig) {
  jets_ = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  muons_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  bRegCorr_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("bRegCorr"));

  produces<std::vector<pat::CompositeCandidate>>();
}

void DijetBuilderProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  edm::Handle<edm::View<pat::Jet>> jets;
  edm::Handle<edm::View<pat::Muon>> muons;
  edm::Handle<edm::ValueMap<float>> bRegCorrHandle;
  iEvent.getByToken(bRegCorr_, bRegCorrHandle);
  iEvent.getByToken(jets_, jets);
  iEvent.getByToken(muons_, muons);

  auto output = std::make_unique<std::vector<pat::CompositeCandidate>>();

  std::vector<unsigned> jetsWithTriggeringMuon;
  std::vector<unsigned> muonIdxs;

  for (size_t i = 0; i < jets->size(); ++i) {
    const auto& jet = jets->at(i);

    if (std::abs(jet.eta()) > 2.5) continue;
    if (jet.pt() < 20) continue;
    //std::cout<<"Jet pT "<< jet.pt() << "jet id "<< jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") <<std::endl;
    if ((jet.pt() < 50) && jet.userInt("pileupJetId:fullId") < 4) continue;
    if (jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") < 6) continue;

    // Muon matching
    int muonIdx1 = -1;
    if (jet.hasOverlaps("muons") && !jet.overlaps("muons").empty()) {
        const auto & muPtr = jet.overlaps("muons")[0];  // this is edm::Ptr<reco::Candidate>
        muonIdx1 = muPtr.key();  // key() is here
    }

int muonIdx2 = -1;
if (jet.hasOverlaps("muons") && jet.overlaps("muons").size() > 1) {
    const auto & muPtr = jet.overlaps("muons")[1];  // edm::Ptr<reco::Candidate>
    muonIdx2 = muPtr.key();  // get index of 2nd muon
}

    //int muonIdx1 = -1;
    //for (size_t jj = 0; jj < muons->size(); ++jj) {
    //    if (deltaR(jet, muons->at(jj)) < 0.4) { // define your overlap criterion
    //        muonIdx1 = jj;
    //        break;
    //    }
    //}
    //int muonIdx2 = -1;
    //for (size_t jj = 0; jj < muons->size(); ++jj) {
    //    if ((deltaR(jet, muons->at(jj)) < 0.4) && (int(jj) != muonIdx1)) {  
    //        muonIdx2 = jj;
    //        break;
    //    }
    //}
      //int muonIdx1 = jet.userInt("muonIdx1");
      //int muonIdx2 = jet.userInt("muonIdx2");

      auto isTriggering = [&](int idx) {
  bool res = (idx >= 0 && idx < (int)muons->size() && muons->at(idx).userInt("isTriggering") == 1);
  //std::cout << "Muon idx: " << idx;
  //if (idx >= 0 && idx < (int)muons->size()) {
    //std::cout << ", isTriggering flag: " << muons->at(idx).userInt("isTriggering");
  //} else {
    //std::cout << " (index out of range)";
  //}
  //std::cout << ", returns: " << res << std::endl;
  return res;
};

if (isTriggering(muonIdx1)) {
  //std::cout << "MuonIdx1 is triggering, jet index: " << i << std::endl;
  jetsWithTriggeringMuon.push_back(i);
  muonIdxs.push_back(muonIdx1);
}
      else if (isTriggering(muonIdx2)) {
        jetsWithTriggeringMuon.push_back(i);
        muonIdxs.push_back(muonIdx2);
      }
  }

  int selected1 = 999, selected2 = 999;

  if (jetsWithTriggeringMuon.size() >= 2) {
    //std::cout<<"2 muons jets "<<jetsWithTriggeringMuon[0]<<" "<<jetsWithTriggeringMuon[1]<<std::endl;
    selected1 = jetsWithTriggeringMuon[0];
    selected2 = jetsWithTriggeringMuon[1];
  }
  else if (jetsWithTriggeringMuon.size() == 1) {
    //std::cout<<"1 muons jets "<<jetsWithTriggeringMuon[0]<<std::endl;
    selected1 = jetsWithTriggeringMuon[0];

    // Now find a second jet based on b-tag
    for (size_t jj = 0; jj < jets->size(); ++jj) {
      if ((int)jj == selected1) continue;
      const auto& jet = jets->at(jj);
      if (jet.pt() < 20) continue;
      if (std::abs(jet.eta()) > 2.5) continue;
      if (jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") < 6) continue;
      if ((jet.pt() < 50) && jet.userInt("pileupJetId:fullId") < 4) continue;
      //std::cout<<"Jet pT "<< jet.pt() << "jet id "<< jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") <<std::endl;

      float btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")+jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+jet.bDiscriminator("pfDeepFlavourJetTags:problepb");

      if (btag > 0.71) {
        selected2 = jj;
        //std::cout<<"Sel2 Tight"<<selected2<<std::endl;
        break;
      }
    }

    // fallback to medium / loose...
    if (selected2 == 999) {
      for (size_t jj = 0; jj < jets->size(); ++jj) {
        if ((int)jj == selected1) continue;
        const auto& jet = jets->at(jj);
        if (jet.pt() < 20) continue;
        if (std::abs(jet.eta()) > 2.5) continue;
        if (jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") < 6) continue;
        if ((jet.pt() < 50) && jet.userInt("pileupJetId:fullId") < 4) continue;
        float btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")+jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        if (btag > 0.2783) {
          selected2 = jj;
          //std::cout<<"Sel2 Medium"<<selected2<<std::endl;
          break;
        }
      }
    }

    if (selected2 == 999) {
      for (size_t jj = 0; jj < jets->size(); ++jj) {
        if ((int)jj == selected1) continue;
        const auto& jet = jets->at(jj);
        if (jet.pt() < 20) continue;
        if (std::abs(jet.eta()) > 2.5) continue;
        if (jet.userInt("tightId")*2+4*jet.userInt("tightIdLepVeto") < 6) continue;
        if ((jet.pt() < 50) && jet.userInt("pileupJetId:fullId") < 4) continue;
        float btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")+jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        if (btag > 0.0490) {
          selected2 = jj;
          //std::cout<<"Sel2 Loose"<<selected2<<std::endl;
          break;
        }
      }
    }

    if (selected2 == 999){
      //std::cout<<"second jets not chosen with btag"<<std::endl;
      selected1 = 999;}
  }

  if (selected1 != 999 && selected2 != 999) {
    const auto& jet0 = jets->at(selected1);
    const auto& jet1 = jets->at(selected2);

    float corr0 = (*bRegCorrHandle)[jets->refAt(selected1)];
    float corr1 = (*bRegCorrHandle)[jets->refAt(selected2)];

    // Apply only to pt
    TLorentzVector vec0, vec1;
    vec0.SetPtEtaPhiE(jet0.pt() * corr0, jet0.eta(), jet0.phi(), jet0.p4().E() * corr0);
    vec1.SetPtEtaPhiE(jet1.pt() * corr1, jet1.eta(), jet1.phi(), jet1.p4().E() * corr1);



    reco::Candidate::LorentzVector p4_0(vec0.Px(), vec0.Py(), vec0.Pz(), vec0.E());
    reco::Candidate::LorentzVector p4_1(vec1.Px(), vec1.Py(), vec1.Pz(), vec1.E());
    //std::cout<<"pt corr : "<<jet0.pt() * corr0<<"  "<<jet1.pt() * corr1<<std::endl;
    //std::cout<<"pt  : "<<jet0.pt() <<"  "<<jet1.pt() <<std::endl;
    //std::cout<<"eta  : "<<jet0.eta() <<"  "<<jet1.eta() <<std::endl;
    //std::cout<<"index :"<<selected1<<" "<<selected2<<std::endl;

    pat::CompositeCandidate dijet;
    dijet.setP4(p4_0 + p4_1);

    output->push_back(dijet);
  }
  else{
    //std::cout<<"Not even one muon"<<std::endl;
    pat::CompositeCandidate dijet;
    dijet.setP4(reco::Candidate::LorentzVector(0., 0., 0., 0.));
    output->push_back(dijet);
  }

  iEvent.put(std::move(output));
}
DEFINE_FWK_MODULE(DijetBuilderProducer);