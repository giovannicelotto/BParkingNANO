#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "correction.h"

class CorrectionLibJECUpdater : public edm::global::EDProducer<> {
	public:
		explicit CorrectionLibJECUpdater(const edm::ParameterSet& iConfig)
			: jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
			correctionFile_(iConfig.getParameter<edm::FileInPath>("correctionFile").fullPath()) {

				// Carica tutte le correzioni dal file JSON
				correctionSet_ = correction::CorrectionSet::from_file(correctionFile_);

				// Ottieni i nomi di tutte le correzioni disponibili
				 for (const auto& correction : *correctionSet_) {
         			   corrections_.push_back(correction.first); // First is the name of the correction
       				 }

				produces<std::vector<pat::Jet>>();
			}

		void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const override {
			edm::Handle<std::vector<pat::Jet>> jets;
			iEvent.getByToken(jetsToken_, jets);
			auto updatedJets = std::make_unique<std::vector<pat::Jet>>(*jets);
			auto systematics = {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample","RelativeStatEC", "RelativeStatFSR","RelativeFSR","RelativeStatHF","PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "Total", "TotalNoFlavor", "TotalNoTime", "TotalNoFlavorNoTime"};
			//auto systematics = {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation","SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta","RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF","RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF","RelativeBal", "RelativeSample", "RelativeFSR", "PileUpDataMC","PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"};

			for (auto& jet : *updatedJets) {
				float pt = jet.pt();
				float eta = jet.eta();
			//std::cout << "Jet Pt: " << pt << ", Jet Eta: " << eta << std::endl;

				for (const auto& corrName : corrections_) {
					auto correction = correctionSet_->at(corrName);
					int count =0;
					std::string label = "";
					for (const auto& sys : systematics){
					if (corrName.find(sys) != std::string::npos) {
						   count++;
						   label = sys;
						   continue;
						}
					}
					if (count !=1) continue;

					float unc = correction->evaluate({eta, pt});
					jet.addUserFloat("jecSys" + label + "Up", unc);
					jet.addUserFloat("jecSys" + label + "Down", -unc);
				}
			}

			iEvent.put(std::move(updatedJets));
		}

	private:
		edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
		std::string correctionFile_;
		std::shared_ptr<const correction::CorrectionSet> correctionSet_;
		std::vector<std::string> corrections_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CorrectionLibJECUpdater);
