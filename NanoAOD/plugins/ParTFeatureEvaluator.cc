#include <TVector3.h>
#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

using namespace btagbtvdeep;

class ParTFeatureEvaluator : public edm::stream::EDProducer<> {
	public:
		explicit ParTFeatureEvaluator(const edm::ParameterSet &);
		~ParTFeatureEvaluator() override;

		static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

	private:

		void beginStream(edm::StreamID) override {}
		void produce(edm::Event &, const edm::EventSetup &) override;
		void endStream() override {}

		void fillParticleFeatures(DeepBoostedJetFeatures & fts, 
				const reco::Jet & jet, 
				const std::vector<math::XYZTLorentzVector> & tau_pfcandidates,  
				const pat::MuonCollection & muons, 
				const pat::ElectronCollection & electrons,
				const pat::PhotonCollection & photons
				);
		void fillSVFeatures(DeepBoostedJetFeatures & fts, 
				const reco::Jet & jet);
		void fillKVFeatures(DeepBoostedJetFeatures & fts, 
				const reco::Jet & jet);
		void fillLVFeatures(DeepBoostedJetFeatures & fts, 
				const reco::Jet & jet);
		void fillLostTrackFeatures(DeepBoostedJetFeatures & fts, 
				const reco::Jet & jet);
		bool useTrackProperties(const pat::PackedCandidate* cand);

		bool isNanOrInf (float i){
			if(std::isnan(i) or std::isinf(i)) return true;
			else return false;
		}

		const double jet_radius_;
		const double min_jet_pt_;
		const double max_jet_eta_;
		const double min_jet_eta_;
		const double min_pt_for_track_properties_;
		const double min_pt_for_pfcandidates_;
		const double min_pt_for_losttrack_;
		const double max_dr_for_losttrack_;
		const double min_pt_for_taus_;
		const double max_eta_for_taus_;
		const bool   use_puppiP4_;
		const double min_puppi_wgt_;
		const bool   flip_ip_sign_;
		const int    ip_sign_;

		edm::EDGetTokenT<pat::MuonCollection > muon_token_;
		edm::EDGetTokenT<pat::ElectronCollection > electron_token_;
		edm::EDGetTokenT<pat::PhotonCollection > photon_token_;
		edm::EDGetTokenT<pat::TauCollection > tau_token_;
		edm::EDGetTokenT<edm::View<reco::Jet> > jet_token_;
		edm::EDGetTokenT<pat::PackedCandidateCollection > losttrack_token_;
		edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
		edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;
		edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> kv_token_;
		edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> lv_token_;
		edm::EDGetTokenT<edm::ValueMap<float> > puppi_weights_token_;
		edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

		edm::Handle<reco::VertexCollection> vtxs_;
		edm::Handle<reco::VertexCompositePtrCandidateCollection> svs_;
		edm::Handle<reco::VertexCompositePtrCandidateCollection> kvs_;
		edm::Handle<reco::VertexCompositePtrCandidateCollection> lvs_;
		edm::Handle<edm::View<reco::Candidate> > pfcands_;
		edm::Handle<pat::PackedCandidateCollection > losttracks_;
		edm::Handle<edm::ValueMap<float> > puppi_weights_;
		edm::ESHandle<TransientTrackBuilder> track_builder_;

		const static std::vector<std::string> particle_ch_features_;
		const static std::vector<std::string> particle_neu_features_;
		const static std::vector<std::string> particle_muon_features_;
		const static std::vector<std::string> particle_electron_features_;
		const static std::vector<std::string> particle_photon_features_;
		const static std::vector<std::string> sv_features_;
		const static std::vector<std::string> kv_features_;
		const static std::vector<std::string> lv_features_;
		const static std::vector<std::string> losttrack_features_;
		const reco::Vertex *pv_ = nullptr;

};

const std::vector<std::string> ParTFeatureEvaluator::particle_ch_features_{
	"jet_pfcand_ch_pt_log","jet_pfcand_ch_energy_log","jet_pfcand_ch_eta","jet_pfcand_ch_dphi","jet_pfcand_ch_deta","jet_pfcand_ch_ptrel","jet_pfcand_ch_etarel","jet_pfcand_ch_pperp","jet_pfcand_ch_ppara","jet_pfcand_ch_id","jet_pfcand_ch_charge","jet_pfcand_ch_frompv","jet_pfcand_ch_dz","jet_pfcand_ch_dzsig","jet_pfcand_ch_dxy","jet_pfcand_ch_dxysig","jet_pfcand_ch_hcalfraction","jet_pfcand_ch_calofraction","jet_pfcand_ch_track_chi2","jet_pfcand_ch_track_qual","jet_pfcand_ch_track_pterr","jet_pfcand_ch_track_etaerr","jet_pfcand_ch_track_phierr","jet_pfcand_ch_trackjet_sip2d","jet_pfcand_ch_trackjet_sip2dsig","jet_pfcand_ch_trackjet_sip3d","jet_pfcand_ch_trackjet_sip3dsig","jet_pfcand_ch_trackjet_dist","jet_pfcand_ch_trackjet_decayL","jet_pfcand_ch_nhits","jet_pfcand_ch_npixhits","jet_pfcand_ch_npixbarrelhits","jet_pfcand_ch_nstriphits","jet_pfcand_ch_nstriptibhits","jet_pfcand_ch_nstriptobhits","jet_pfcand_ch_nlosthits","jet_pfcand_ch_nlayers","jet_pfcand_ch_npixlayers","jet_pfcand_ch_nstriplayers","jet_pfcand_ch_tau_signal","jet_pfcand_ch_pt","jet_pfcand_ch_eta_v","jet_pfcand_ch_phi","jet_pfcand_ch_energy","jet_pfcand_ch_mask"};


const std::vector<std::string> ParTFeatureEvaluator::particle_neu_features_{
	"jet_pfcand_neu_pt_log","jet_pfcand_neu_energy_log","jet_pfcand_neu_eta","jet_pfcand_neu_dphi","jet_pfcand_neu_deta","jet_pfcand_neu_ptrel","jet_pfcand_neu_etarel","jet_pfcand_neu_pperp","jet_pfcand_neu_ppara","jet_pfcand_neu_id","jet_pfcand_neu_frompv","jet_pfcand_neu_dz","jet_pfcand_neu_dxy","jet_pfcand_neu_hcalfraction","jet_pfcand_neu_calofraction","jet_pfcand_neu_tau_signal","jet_pfcand_neu_pt","jet_pfcand_neu_eta_v","jet_pfcand_neu_phi","jet_pfcand_neu_energy","jet_pfcand_neu_mask"};

const std::vector<std::string> ParTFeatureEvaluator::particle_muon_features_{
	"jet_pfcand_muon_pt_log","jet_pfcand_muon_energy_log","jet_pfcand_muon_eta","jet_pfcand_muon_dphi","jet_pfcand_muon_deta","jet_pfcand_muon_chi2","jet_pfcand_muon_dxy","jet_pfcand_muon_dz","jet_pfcand_muon_segcomp","jet_pfcand_muon_validfraction","jet_pfcand_muon_trkKink","jet_pfcand_muon_pterr","jet_pfcand_muon_type","jet_pfcand_muon_nvalidhits","jet_pfcand_muon_nstations","jet_pfcand_muon_nlayers","jet_pfcand_muon_nhits","jet_pfcand_muon_pt","jet_pfcand_muon_eta_v","jet_pfcand_muon_phi","jet_pfcand_muon_energy","jet_pfcand_muon_mask"};

const std::vector<std::string> ParTFeatureEvaluator::particle_electron_features_{
	"jet_pfcand_electron_pt_log","jet_pfcand_electron_energy_log","jet_pfcand_electron_eta","jet_pfcand_electron_dphi","jet_pfcand_electron_deta","jet_pfcand_electron_dxy","jet_pfcand_electron_dz","jet_pfcand_electron_eOverP","jet_pfcand_electron_detaIn","jet_pfcand_electron_dphiIn","jet_pfcand_electron_r9","jet_pfcand_electron_hOverE","jet_pfcand_electron_sigIetaIeta","jet_pfcand_electron_convProb","jet_pfcand_electron_pt","jet_pfcand_electron_eta_v","jet_pfcand_electron_phi","jet_pfcand_electron_energy","jet_pfcand_electron_mask"};

const std::vector<std::string> ParTFeatureEvaluator::particle_photon_features_{
	"jet_pfcand_photon_pt_log","jet_pfcand_photon_energy_log","jet_pfcand_photon_eta","jet_pfcand_photon_dphi","jet_pfcand_photon_deta","jet_pfcand_photon_r9","jet_pfcand_photon_hOverE","jet_pfcand_photon_sigIetaIeta","jet_pfcand_photon_eVeto","jet_pfcand_photon_pt","jet_pfcand_photon_eta_v","jet_pfcand_photon_phi","jet_pfcand_photon_energy","jet_pfcand_photon_mask"};

const std::vector<std::string> ParTFeatureEvaluator::sv_features_{"jet_sv_pt_log","jet_sv_energy_log","jet_sv_eta","jet_sv_mass","jet_sv_deta","jet_sv_dphi", "jet_sv_eta", "jet_sv_ntrack", "jet_sv_chi2", "jet_sv_dxy", "jet_sv_dxysig", "jet_sv_d3d", "jet_sv_d3dsig","jet_sv_d3d_sign","jet_sv_costheta","jet_sv_pt","jet_sv_eta_v","jet_sv_phi","jet_sv_energy","jet_sv_mask"};

const std::vector<std::string> ParTFeatureEvaluator::kv_features_{"jet_kaon_pt_log","jet_kaon_energy_log","jet_kaon_eta","jet_kaon_mass","jet_kaon_deta","jet_kaon_dphi", "jet_kaon_eta", "jet_kaon_ntrack", "jet_kaon_chi2", "jet_kaon_dxy", "jet_kaon_dxysig", "jet_kaon_d3d", "jet_kaon_d3dsig","jet_kaon_d3d_sign","jet_kaon_costheta","jet_kaon_pt","jet_kaon_eta_v","jet_kaon_phi","jet_kaon_energy","jet_kaon_mask"};

const std::vector<std::string> ParTFeatureEvaluator::lv_features_{"jet_lambda_pt_log","jet_lambda_energy_log","jet_lambda_eta","jet_lambda_mass","jet_lambda_deta","jet_lambda_dphi", "jet_lambda_eta", "jet_lambda_ntrack", "jet_lambda_chi2", "jet_lambda_dxy", "jet_lambda_dxysig", "jet_lambda_d3d", "jet_lambda_d3dsig","jet_lambda_d3d_sign","jet_lambda_costheta","jet_lambda_pt","jet_lambda_eta_v","jet_lambda_phi","jet_lambda_energy","jet_lambda_mask"};

const std::vector<std::string> ParTFeatureEvaluator::losttrack_features_{
	"jet_losttrack_pt_log","jet_losttrack_energy_log","jet_losttrack_eta","jet_losttrack_dphi","jet_losttrack_deta","jet_losttrack_ptrel","jet_losttrack_etarel","jet_losttrack_charge","jet_losttrack_frompv","jet_losttrack_dz","jet_losttrack_dzsig","jet_losttrack_dxy","jet_losttrack_dxysig","jet_losttrack_track_chi2","jet_losttrack_track_qual","jet_losttrack_track_pterr","jet_losttrack_track_etaerr","jet_losttrack_track_phierr","jet_losttrack_trackjet_sip2d","jet_losttrack_trackjet_sip2dsig","jet_losttrack_trackjet_sip3d","jet_losttrack_trackjet_sip3dsig","jet_losttrack_trackjet_dist","jet_losttrack_trackjet_decayL","jet_losttrack_nhits","jet_losttrack_npixhits","jet_losttrack_npixbarrelhits","jet_losttrack_nstriphits","jet_losttrack_nstriptibhits","jet_losttrack_nstriptobhits","jet_losttrack_nlosthits","jet_losttrack_nlayers","jet_losttrack_npixlayers","jet_losttrack_nstriplayers","jet_losttrack_pt","jet_losttrack_eta_v","jet_losttrack_phi","jet_losttrack_energy","jet_losttrack_mask"};

ParTFeatureEvaluator::ParTFeatureEvaluator(const edm::ParameterSet &iConfig):
	jet_radius_(iConfig.getParameter<double>("jet_radius")),
	min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
	max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
	min_jet_eta_(iConfig.getParameter<double>("min_jet_eta")),
	min_pt_for_track_properties_(iConfig.getParameter<double>("min_pt_for_track_properties")),
	min_pt_for_pfcandidates_(iConfig.getParameter<double>("min_pt_for_pfcandidates")),
	min_pt_for_losttrack_(iConfig.getParameter<double>("min_pt_for_losttrack")),
	max_dr_for_losttrack_(iConfig.getParameter<double>("max_dr_for_losttrack")),
	min_pt_for_taus_(iConfig.getParameter<double>("min_pt_for_taus")),
	max_eta_for_taus_(iConfig.getParameter<double>("max_eta_for_taus")),
	use_puppiP4_(iConfig.getParameter<bool>("use_puppiP4")),
	min_puppi_wgt_(iConfig.getParameter<double>("min_puppi_wgt")),
	flip_ip_sign_(iConfig.getParameter<bool>("flip_ip_sign")),
	ip_sign_((flip_ip_sign_) ? -1 : 1),
	muon_token_(consumes<pat::MuonCollection >(iConfig.getParameter<edm::InputTag>("muons"))),
	electron_token_(consumes<pat::ElectronCollection >(iConfig.getParameter<edm::InputTag>("electrons"))),
	photon_token_(consumes<pat::PhotonCollection >(iConfig.getParameter<edm::InputTag>("photons"))),
	tau_token_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jet_token_(consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
	losttrack_token_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("losttracks"))),
	vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
	kv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("kaon_vertices"))),
	lv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("lambda_vertices"))),
	puppi_weights_token_(mayConsume<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppi_weights"))),
	track_builder_token_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))) {

		produces<std::vector<reco::DeepBoostedJetTagInfo>>();

	}

ParTFeatureEvaluator::~ParTFeatureEvaluator() {}

void ParTFeatureEvaluator::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
	// pfDeepBoostedJetTagInfos
	edm::ParameterSetDescription desc;
	desc.add<double>("jet_radius", 0.8);
	desc.add<double>("min_jet_pt", 150);
	desc.add<double>("max_jet_eta", 99);
	desc.add<double>("min_jet_eta", 0.0);
	desc.add<double>("min_pt_for_track_properties", -1);
	desc.add<double>("min_pt_for_pfcandidates", -1);
	desc.add<double>("min_pt_for_losttrack", 1);
	desc.add<double>("max_dr_for_losttrack", 0.4);
	desc.add<double>("min_pt_for_taus",20.);
	desc.add<double>("max_eta_for_taus",2.4);
	desc.add<bool>("use_puppiP4", false);
	desc.add<bool>("flip_ip_sign", false);
	desc.add<double>("min_puppi_wgt", -1.);
	desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
	desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
	desc.add<edm::InputTag>("kaon_vertices", edm::InputTag("slimmedKshortVertices"));
	desc.add<edm::InputTag>("lambda_vertices", edm::InputTag("slimmedLambdaVertices"));
	desc.add<edm::InputTag>("losttracks", edm::InputTag("lostTracks"));
	desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
	desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
	desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
	desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
	desc.add<edm::InputTag>("photons", edm::InputTag("slimmedPhotons"));
	desc.add<edm::InputTag>("puppi_weights", edm::InputTag(""));
	descriptions.add("ParTFeatureEvaluator", desc);
}

void ParTFeatureEvaluator::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {

	// output collection
	auto output_tag_infos = std::make_unique<std::vector<reco::DeepBoostedJetTagInfo>>();
	// Input jets
	auto jets = iEvent.getHandle(jet_token_);
	// Input muons
	auto muons = iEvent.getHandle(muon_token_);  
	// Input taus
	auto taus = iEvent.getHandle(tau_token_);  
	// Input electrons
	auto electrons = iEvent.getHandle(electron_token_);
	// Input photons
	auto photons = iEvent.getHandle(photon_token_);
	// Input lost tracks
	iEvent.getByToken(losttrack_token_, losttracks_);
	// Primary vertexes
	iEvent.getByToken(vtx_token_, vtxs_);
	if (vtxs_->empty()) {
		// produce empty TagInfos in case no primary vertex
		iEvent.put(std::move(output_tag_infos));
		return;  // exit event
	}
	// Leading vertex
	pv_ = &vtxs_->at(0);
	// Secondary vertexs
	iEvent.getByToken(sv_token_, svs_);
	// kaons vertexs
	iEvent.getByToken(kv_token_, kvs_);
	// lambdas vertexs
	iEvent.getByToken(lv_token_, lvs_);
	// puppi value map
	iEvent.getByToken(puppi_weights_token_, puppi_weights_);
	// Track builder
	track_builder_ = iSetup.getHandle(track_builder_token_);

	// tau signal candidates
	std::vector<math::XYZTLorentzVector> tau_pfcandidates;
	for (size_t itau = 0; itau < taus->size(); itau++) {
		if(taus->at(itau).pt() < min_pt_for_taus_) continue;
		if(fabs(taus->at(itau).eta()) > max_eta_for_taus_) continue;
		for(unsigned ipart = 0; ipart < taus->at(itau).signalCands().size(); ipart++){
			const pat::PackedCandidate* pfcand = dynamic_cast<const pat::PackedCandidate*> (taus->at(itau).signalCands()[ipart].get());
			tau_pfcandidates.push_back(pfcand->p4());
		}
	}

	// Loop over jet
	for (std::size_t jet_n = 0; jet_n < jets->size(); jet_n++) {

		const auto & jet = (*jets)[jet_n];        
		edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);

		// create jet features
		DeepBoostedJetFeatures features;
		for (const auto &name : particle_ch_features_) 
			features.add(name);
		for (const auto &name : particle_neu_features_) 
			features.add(name);
		for (const auto &name : particle_muon_features_) 
			features.add(name);
		for (const auto &name : particle_electron_features_) 
			features.add(name);
		for (const auto &name : particle_photon_features_) 
			features.add(name);
		for (const auto &name : sv_features_)
			features.add(name);
		for (const auto &name : kv_features_)
			features.add(name);
		for (const auto &name : lv_features_)
			features.add(name);
		for (const auto &name : losttrack_features_)
			features.add(name);

		// fill values only if above pt threshold and has daughters, otherwise left
		bool fill_vars = true;
		if ((jet.pt() < min_jet_pt_ and dynamic_cast<const pat::Jet *>(&jet)->correctedJet("Uncorrected").pt() < min_jet_pt_) or 
				std::abs(jet.eta()) > max_jet_eta_ or 
				std::abs(jet.eta()) < min_jet_eta_ )
			fill_vars = false;
		if (jet.numberOfDaughters() == 0)
			fill_vars = false;

		// fill features
		if (fill_vars) {
			fillParticleFeatures(features, jet, tau_pfcandidates, *muons, *electrons, *photons);
			fillSVFeatures(features, jet);
			fillKVFeatures(features, jet);
			fillLVFeatures(features, jet);
			fillLostTrackFeatures(features, jet);
			features.check_consistency(particle_ch_features_);
			features.check_consistency(particle_neu_features_);
			features.check_consistency(particle_muon_features_);
			features.check_consistency(particle_electron_features_);
			features.check_consistency(particle_photon_features_);
			features.check_consistency(sv_features_);
			features.check_consistency(kv_features_);
			features.check_consistency(lv_features_);
			features.check_consistency(losttrack_features_);
		}

		// this should always be done even if features are not filled
		output_tag_infos->emplace_back(features, jet_ref);

	}
	// move output collection
	iEvent.put(std::move(output_tag_infos));
}

bool ParTFeatureEvaluator::useTrackProperties(const pat::PackedCandidate* cand) {
	const auto *track = cand->bestTrack();
	return track != nullptr and track->pt() > min_pt_for_track_properties_;
};

void ParTFeatureEvaluator::fillParticleFeatures(DeepBoostedJetFeatures & fts, 
		const reco::Jet & jet,
		const std::vector<math::XYZTLorentzVector> & tau_pfcandidates,
		const pat::MuonCollection & muons, 
		const pat::ElectronCollection & electrons,
		const pat::PhotonCollection & photons
		) {
	// some jet properties
	math::XYZVector jet_dir = jet.momentum().Unit();
	TVector3 jet_direction (jet.momentum().Unit().x(), jet.momentum().Unit().y(), jet.momentum().Unit().z());
	GlobalVector jet_ref_track_dir (jet.px(), jet.py(), jet.pz());

	// vertexes
	reco::VertexRef pv_ass = reco::VertexRef(vtxs_, 0);
	math::XYZPoint  pv_ass_pos = pv_ass->position();

	// make list of pf-candidates to be considered
	std::vector<const pat::PackedCandidate*> daughters;
	std::vector<float> puppiWeights;

	// build the list of interesting pf-candidates in the jet
	for (const auto & dau : jet.daughterPtrVector()) {
		const pat::PackedCandidate* cand = dynamic_cast<const pat::PackedCandidate*>(&(*dau));
		if(not cand)
			throw edm::Exception(edm::errors::InvalidReference)
				<< "Cannot convert to either pat::PackedCandidate";
		// remove low puppi weight 
		if(puppi_weights_.isValid() and (*puppi_weights_)[dau] < min_puppi_wgt_) continue;
		else if(not puppi_weights_.isValid() and cand->puppiWeight() < min_puppi_wgt_) continue;
		// base requirements on PF candidates
		if (cand->pt() < min_pt_for_pfcandidates_) continue;
		// filling daughters
		daughters.push_back(cand);
		puppiWeights.push_back(puppi_weights_.isValid() ? (*puppi_weights_)[dau] : cand->puppiWeight());
	}

	// sort by Puppi-weighted pt
	std::vector<size_t> daughterIndices (daughters.size());
	std::iota(daughterIndices.begin(), daughterIndices.end(), 0);
	if(use_puppiP4_)
		std::sort(daughterIndices.begin(), daughterIndices.end(), [&](const auto & i1, const auto & i2) {
				return puppiWeights[i1]*daughters[i1]->pt() > puppiWeights[i2]*daughters[i2]->pt();
				});    
	else
		std::sort(daughterIndices.begin(), daughterIndices.end(), [&](const auto &i1, const auto & i2) {
				return daughters[i1]->pt() > daughters[i2]->pt(); 
				});


	// build sub daughter list of particles
	std::vector<size_t> daughters_ch_index;
	std::vector<size_t> daughters_neu_index;
	std::vector<size_t> daughters_muon_index;
	std::vector<size_t> daughters_electron_index;
	std::vector<size_t> daughters_photon_index;

	for (const auto & icand : daughterIndices) {    
		auto cand = daughters[icand];
		if(cand->charge() != 0){
			daughters_ch_index.push_back(icand);
			// muon specific
			if(abs(cand->pdgId()) == 13){
				int ipos = -1;
				float minDR = 1000;
				for (size_t imu = 0; imu < muons.size(); imu++) {
					if(std::find(daughters_muon_index.begin(),daughters_muon_index.end(),imu) != daughters_muon_index.end()) continue;
					if(not muons[imu].isPFMuon()) continue;
					float dR = reco::deltaR(muons[imu].p4(),cand->p4());
					if(dR < jet_radius_ and dR < minDR){
						minDR = dR;
						ipos = imu;
					}
				}	
				if(ipos >= 0)
					daughters_muon_index.push_back(ipos);
			}
			// electron specific
			else if(abs(cand->pdgId()) == 11){	
				int ipos = -1;
				for (size_t iel = 0; iel < electrons.size(); iel++) {
					if(std::find(daughters_electron_index.begin(),daughters_electron_index.end(),iel) != daughters_electron_index.end()) continue;
					if(electrons[iel].isPF()){
						for(const auto & element : electrons[iel].associatedPackedPFCandidates()){
							if(abs(element->pdgId()) == 11 and element->p4() == cand->p4()){
								ipos = iel;
								break;
							}
							if(ipos != -1) break;
						}
					}
				}	  
				if(ipos >= 0)
					daughters_electron_index.push_back(ipos);
			}
		}
		else{
			daughters_neu_index.push_back(icand);
			// photon specific
			if(abs(cand->pdgId()) == 22){	
				int ipos = -1;
				for (size_t iph = 0; iph < photons.size(); iph++) {
					if(std::find(daughters_photon_index.begin(),daughters_photon_index.end(),iph) != daughters_photon_index.end()) continue;
					for(const auto & element : photons[iph].associatedPackedPFCandidates()){
						if(abs(element->pdgId()) == 22 and element->p4() == cand->p4()){
							ipos = iph;
							break;
						}
						if(ipos != -1) break;
					}
				}
				if(ipos >= 0)
					daughters_photon_index.push_back(ipos);
			}
		}
	}

	// reserve space
	for (const auto &name : particle_ch_features_)
		fts.reserve(name, daughters_ch_index.size());
	for (const auto &name : particle_neu_features_)
		fts.reserve(name, daughters_neu_index.size());
	for (const auto &name : particle_muon_features_)
		fts.reserve(name, daughters_muon_index.size());
	for (const auto &name : particle_electron_features_)
		fts.reserve(name, daughters_electron_index.size());
	for (const auto &name : particle_photon_features_)
		fts.reserve(name, daughters_photon_index.size());

	// Build observables ch pf
	for (const auto & icand : daughters_ch_index){
		auto cand = daughters[icand];
		auto puppiw = puppiWeights[icand];

		// input particle is a packed PF candidate
		auto candP4 = use_puppiP4_ ? puppiw * cand->p4() : cand->p4();
		auto candP3 = use_puppiP4_ ? puppiw * cand->momentum() : cand->momentum();
		TVector3 cand_direction(candP3.x(), candP3.y(), candP3.z());

		// candidate track
		const reco::Track *track = nullptr;
		if(useTrackProperties(cand))
			track = cand->bestTrack();

		fts.fill("jet_pfcand_ch_pt", isNanOrInf(cand->pt()) ? 0 : cand->pt());
		fts.fill("jet_pfcand_ch_pt_log", isNanOrInf(std::log(cand->pt())) ? 0 : std::log(cand->pt()));
		fts.fill("jet_pfcand_ch_eta", isNanOrInf(cand->eta()) ? 0 : cand->eta());
		fts.fill("jet_pfcand_ch_eta_v", isNanOrInf(cand->eta()) ? 0 : cand->eta());
		fts.fill("jet_pfcand_ch_phi", isNanOrInf(cand->phi()) ? 0 : cand->phi());
		fts.fill("jet_pfcand_ch_energy",isNanOrInf(cand->energy()) ? 0 : cand->energy());
		fts.fill("jet_pfcand_ch_energy_log", isNanOrInf(std::log(cand->energy())) ? 0 : std::log(cand->energy()));
		fts.fill("jet_pfcand_ch_deta",isNanOrInf(cand->eta()-jet.eta()) ? 0 : cand->eta()-jet.eta());
		fts.fill("jet_pfcand_ch_dphi",isNanOrInf(reco::deltaPhi(cand->phi(),jet.phi())) ? 0 : reco::deltaPhi(cand->phi(),jet.phi()));
		fts.fill("jet_pfcand_ch_calofraction",isNanOrInf(cand->caloFraction()) ? 0 : cand->caloFraction());
		fts.fill("jet_pfcand_ch_hcalfraction",isNanOrInf(cand->hcalFraction()) ? 0 : cand->hcalFraction());
		fts.fill("jet_pfcand_ch_dxy",ip_sign_*(isNanOrInf(cand->dxy(pv_ass_pos)) ? 0 : cand->dxy(pv_ass_pos)));
		fts.fill("jet_pfcand_ch_dz",ip_sign_*(isNanOrInf(cand->dz(pv_ass_pos)) ? 0 : cand->dz(pv_ass_pos)));
		fts.fill("jet_pfcand_ch_frompv",isNanOrInf(cand->fromPV()) ? 0 : cand->fromPV());
		fts.fill("jet_pfcand_ch_pperp",isNanOrInf(jet_direction.Perp(cand_direction)/cand_direction.Mag()) ? 0 : jet_direction.Perp(cand_direction)/cand_direction.Mag());
		fts.fill("jet_pfcand_ch_ppara",isNanOrInf(jet_direction.Dot(cand_direction)/cand_direction.Mag()) ? 0 : jet_direction.Dot(cand_direction)/cand_direction.Mag());
		fts.fill("jet_pfcand_ch_ptrel",isNanOrInf(cand_direction.Perp(jet_direction)) ? 0 : cand_direction.Perp(jet_direction));
		fts.fill("jet_pfcand_ch_etarel",isNanOrInf(reco::btau::etaRel(jet_dir,cand->momentum())) ? 0 : reco::btau::etaRel(jet_dir,cand->momentum()));
		fts.fill("jet_pfcand_ch_charge",isNanOrInf(cand->charge()) ? 0 : cand->charge());
		fts.fill("jet_pfcand_ch_nlosthits",isNanOrInf(cand->lostInnerHits()) ? 0 : cand->lostInnerHits());

		if(std::find(tau_pfcandidates.begin(),tau_pfcandidates.end(),candP4) != tau_pfcandidates.end())
			fts.fill("jet_pfcand_ch_tau_signal",1);
		else
			fts.fill("jet_pfcand_ch_tau_signal",0);

		if(abs(cand->pdgId()) == 11)
			fts.fill("jet_pfcand_ch_id",0);
		else if(abs(cand->pdgId()) == 13)
			fts.fill("jet_pfcand_ch_id",1);
		else if(abs(cand->pdgId()) != 11 and abs(cand->pdgId()) != 13 and cand->charge()!=0)
			fts.fill("jet_pfcand_ch_id",2);
		else
			fts.fill("jet_pfcand_ch_id",-1);

		if(track){
			fts.fill("jet_pfcand_ch_dxysig",ip_sign_*(isNanOrInf(cand->dxy(pv_ass_pos)/cand->dxyError()) ? 0 : cand->dxy(pv_ass_pos)/cand->dxyError()));
			fts.fill("jet_pfcand_ch_dzsig",ip_sign_*(isNanOrInf(cand->dz(pv_ass_pos)/cand->dzError()) ? 0 : cand->dz(pv_ass_pos)/cand->dzError()));
			fts.fill("jet_pfcand_ch_track_chi2",isNanOrInf(track->normalizedChi2()) ? 0 : track->normalizedChi2());
			fts.fill("jet_pfcand_ch_track_qual",isNanOrInf(track->qualityMask()) ? 0 : track->qualityMask());
			fts.fill("jet_pfcand_ch_track_pterr",isNanOrInf(track->ptError()/track->pt()) ? 0 : track->ptError()/track->pt());
			fts.fill("jet_pfcand_ch_track_etaerr",isNanOrInf(track->etaError()) ? 0 : track->etaError());
			fts.fill("jet_pfcand_ch_track_phierr",isNanOrInf(track->phiError()) ? 0 : track->phiError());
			fts.fill("jet_pfcand_ch_nhits",isNanOrInf(track->hitPattern().numberOfValidHits()) ? 0 : track->hitPattern().numberOfValidHits());
			fts.fill("jet_pfcand_ch_npixhits",isNanOrInf(track->hitPattern().numberOfValidPixelHits()) ? 0 : track->hitPattern().numberOfValidPixelHits());
			fts.fill("jet_pfcand_ch_npixbarrelhits",isNanOrInf(track->hitPattern().numberOfValidPixelBarrelHits()) ? 0 : track->hitPattern().numberOfValidPixelBarrelHits());
			fts.fill("jet_pfcand_ch_nstriphits",isNanOrInf(track->hitPattern().numberOfValidStripHits()) ? 0 : track->hitPattern().numberOfValidStripHits());
			fts.fill("jet_pfcand_ch_nstriptibhits",isNanOrInf(track->hitPattern().numberOfValidStripTIBHits()+track->hitPattern().numberOfValidStripTIDHits()) ? 0 : track->hitPattern().numberOfValidStripTIBHits()+track->hitPattern().numberOfValidStripTIDHits());
			fts.fill("jet_pfcand_ch_nstriptobhits",isNanOrInf(track->hitPattern().numberOfValidStripTOBHits()) ? 0 : track->hitPattern().numberOfValidStripTOBHits());
			fts.fill("jet_pfcand_ch_nlayers",isNanOrInf(track->hitPattern().trackerLayersWithMeasurement()) ? 0 : track->hitPattern().trackerLayersWithMeasurement());
			fts.fill("jet_pfcand_ch_npixlayers",isNanOrInf(track->hitPattern().pixelLayersWithMeasurement()) ? 0 : track->hitPattern().pixelLayersWithMeasurement());
			fts.fill("jet_pfcand_ch_nstriplayers",isNanOrInf(track->hitPattern().stripLayersWithMeasurement()) ? 0 : track->hitPattern().stripLayersWithMeasurement());
			reco::TransientTrack transientTrack = track_builder_->build(*track);
			Measurement1D meas_ip2d    = IPTools::signedTransverseImpactParameter(transientTrack,jet_ref_track_dir, *pv_).second;
			Measurement1D meas_ip3d    = IPTools::signedImpactParameter3D(transientTrack, jet_ref_track_dir, *pv_).second;
			Measurement1D meas_jetdist = IPTools::jetTrackDistance(transientTrack, jet_ref_track_dir, *pv_).second;
			Measurement1D meas_decayl  = IPTools::signedDecayLength3D(transientTrack, jet_ref_track_dir, *pv_).second;
			fts.fill("jet_pfcand_ch_trackjet_sip2d",ip_sign_*(isNanOrInf(meas_ip2d.value()) ? 0 : meas_ip2d.value()));
			fts.fill("jet_pfcand_ch_trackjet_sip3d",ip_sign_*(isNanOrInf(meas_ip3d.value()) ? 0 : meas_ip3d.value()));
			fts.fill("jet_pfcand_ch_trackjet_sip2dsig",ip_sign_*(isNanOrInf(meas_ip2d.significance()) ? 0 : meas_ip2d.significance()));
			fts.fill("jet_pfcand_ch_trackjet_sip3dsig",ip_sign_*(isNanOrInf(meas_ip3d.significance()) ? 0 : meas_ip3d.significance()));
			fts.fill("jet_pfcand_ch_trackjet_dist",isNanOrInf(-meas_jetdist.value()) ? 0 : -meas_jetdist.value());
			fts.fill("jet_pfcand_ch_trackjet_decayL",isNanOrInf(meas_decayl.value()) ? 0 : meas_decayl.value());
		}
		else{
			fts.fill("jet_pfcand_ch_dxysig",0);
			fts.fill("jet_pfcand_ch_dzsig",0);
			fts.fill("jet_pfcand_ch_track_chi2",0);
			fts.fill("jet_pfcand_ch_track_qual",0);
			fts.fill("jet_pfcand_ch_track_pterr",0);
			fts.fill("jet_pfcand_ch_track_etaerr",0);
			fts.fill("jet_pfcand_ch_track_phierr",0);
			fts.fill("jet_pfcand_ch_nhits",0);
			fts.fill("jet_pfcand_ch_npixhits",0);
			fts.fill("jet_pfcand_ch_npixbarrelhits",0);
			fts.fill("jet_pfcand_ch_nstriphits",0);
			fts.fill("jet_pfcand_ch_nstriptibhits",0);
			fts.fill("jet_pfcand_ch_nstriptobhits",0);
			fts.fill("jet_pfcand_ch_nlayers",0);
			fts.fill("jet_pfcand_ch_npixlayers",0);
			fts.fill("jet_pfcand_ch_nstriplayers",0);
			fts.fill("jet_pfcand_ch_trackjet_sip2d",0);
			fts.fill("jet_pfcand_ch_trackjet_sip3d",0);
			fts.fill("jet_pfcand_ch_trackjet_sip2dsig",0);
			fts.fill("jet_pfcand_ch_trackjet_sip3dsig",0);
			fts.fill("jet_pfcand_ch_trackjet_dist",0);
			fts.fill("jet_pfcand_ch_trackjet_decayL",0);
		}
		fts.fill("jet_pfcand_ch_mask",1);
	}

	// Build observables neu pf
	for (const auto & icand : daughters_neu_index){
		auto cand = daughters[icand];
		auto puppiw = puppiWeights[icand];

		// input particle is a packed PF candidate
		auto candP3 = use_puppiP4_ ? puppiw * cand->momentum() : cand->momentum();
		TVector3 cand_direction(candP3.x(), candP3.y(), candP3.z());

		fts.fill("jet_pfcand_neu_pt", isNanOrInf(cand->pt()) ? 0 : cand->pt());
		fts.fill("jet_pfcand_neu_pt_log", isNanOrInf(std::log(cand->pt())) ? 0 : std::log(cand->pt()));
		fts.fill("jet_pfcand_neu_eta", isNanOrInf(cand->eta()) ? 0 : cand->eta());
		fts.fill("jet_pfcand_neu_eta_v", isNanOrInf(cand->eta()) ? 0 : cand->eta());
		fts.fill("jet_pfcand_neu_phi", isNanOrInf(cand->phi()) ? 0 : cand->phi());
		fts.fill("jet_pfcand_neu_energy",isNanOrInf(cand->energy()) ? 0 : cand->energy());
		fts.fill("jet_pfcand_neu_energy_log", isNanOrInf(std::log(cand->energy())) ? 0 : std::log(cand->energy()));
		fts.fill("jet_pfcand_neu_deta",isNanOrInf(cand->eta()-jet.eta()) ? 0 : cand->eta()-jet.eta());
		fts.fill("jet_pfcand_neu_dphi",isNanOrInf(reco::deltaPhi(cand->phi(),jet.phi())) ? 0 : reco::deltaPhi(cand->phi(),jet.phi()));
		fts.fill("jet_pfcand_neu_calofraction",isNanOrInf(cand->caloFraction()) ? 0 : cand->caloFraction());
		fts.fill("jet_pfcand_neu_hcalfraction",isNanOrInf(cand->hcalFraction()) ? 0 : cand->hcalFraction());
		fts.fill("jet_pfcand_neu_dxy",ip_sign_*(isNanOrInf(cand->dxy(pv_ass_pos)) ? 0 : cand->dxy(pv_ass_pos)));
		fts.fill("jet_pfcand_neu_dz",ip_sign_*(isNanOrInf(cand->dz(pv_ass_pos)) ? 0 : cand->dz(pv_ass_pos)));
		fts.fill("jet_pfcand_neu_frompv",isNanOrInf(cand->fromPV()) ? 0 : cand->fromPV());
		fts.fill("jet_pfcand_neu_pperp",isNanOrInf(jet_direction.Perp(cand_direction)/cand_direction.Mag()) ? 0 : jet_direction.Perp(cand_direction)/cand_direction.Mag());
		fts.fill("jet_pfcand_neu_ppara",isNanOrInf(jet_direction.Dot(cand_direction)/cand_direction.Mag()) ? 0 : jet_direction.Dot(cand_direction)/cand_direction.Mag());
		fts.fill("jet_pfcand_neu_ptrel",isNanOrInf(cand_direction.Perp(jet_direction)) ? 0 : cand_direction.Perp(jet_direction));
		fts.fill("jet_pfcand_neu_etarel",isNanOrInf(reco::btau::etaRel(jet_dir,cand->momentum())) ? 0 : reco::btau::etaRel(jet_dir,cand->momentum()));

		// tau specific
		if(std::find(tau_pfcandidates.begin(),tau_pfcandidates.end(),cand->p4()) != tau_pfcandidates.end())
			fts.fill("jet_pfcand_neu_tau_signal",1);
		else
			fts.fill("jet_pfcand_neu_tau_signal",0);

		if(abs(cand->pdgId()) == 22)
			fts.fill("jet_pfcand_neu_id",0);
		else if(abs(cand->pdgId()) != 22 and abs(cand->pdgId()) != 1 and abs(cand->pdgId()) != 2)
			fts.fill("jet_pfcand_neu_id",1);
		else if(abs(cand->pdgId()) == 1)
			fts.fill("jet_pfcand_neu_id",2);
		else if(abs(cand->pdgId()) == 2)
			fts.fill("jet_pfcand_neu_id",3);
		else
			fts.fill("jet_pfcand_neu_id",-1);

		fts.fill("jet_pfcand_neu_mask",1.);

	}

	// Build observables muon info
	for (const auto & icand : daughters_muon_index){
		auto cand = muons[icand];
		fts.fill("jet_pfcand_muon_pt", isNanOrInf(cand.pt()) ? 0 : cand.pt());
		fts.fill("jet_pfcand_muon_pt_log", isNanOrInf(std::log(cand.pt())) ? 0 : std::log(cand.pt()));
		fts.fill("jet_pfcand_muon_eta", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_muon_eta_v", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_muon_phi", isNanOrInf(cand.phi()) ? 0 : cand.phi());
		fts.fill("jet_pfcand_muon_energy",isNanOrInf(cand.energy()) ? 0 : cand.energy());
		fts.fill("jet_pfcand_muon_energy_log", isNanOrInf(std::log(cand.energy())) ? 0 : std::log(cand.energy()));
		fts.fill("jet_pfcand_muon_deta",isNanOrInf(cand.eta()-jet.eta()) ? 0 : cand.eta()-jet.eta());
		fts.fill("jet_pfcand_muon_dphi",isNanOrInf(reco::deltaPhi(cand.phi(),jet.phi())) ? 0 : reco::deltaPhi(cand.phi(),jet.phi()));

		int type = 0;
		if(cand.isStandAloneMuon()) type += 1;
		if(cand.isTrackerMuon()) type += 2;
		if(cand.isPFMuon()) type += 4;
		if(cand.isGlobalMuon()) type += 8;
		fts.fill("jet_pfcand_muon_type",type);
		if(cand.isGlobalMuon()){
			fts.fill("jet_pfcand_muon_nvalidhits",isNanOrInf(cand.globalTrack()->hitPattern().numberOfValidMuonHits()) ? 0 : cand.globalTrack()->hitPattern().numberOfValidMuonHits());
			fts.fill("jet_pfcand_muon_chi2",isNanOrInf(cand.globalTrack()->normalizedChi2()) ? 0 : cand.globalTrack()->normalizedChi2());
		}
		else if(cand.innerTrack().isNonnull()){
			fts.fill("jet_pfcand_muon_nvalidhits",isNanOrInf(cand.innerTrack()->hitPattern().numberOfValidMuonHits()) ? 0 : cand.innerTrack()->hitPattern().numberOfValidMuonHits());
			fts.fill("jet_pfcand_muon_chi2",isNanOrInf(cand.innerTrack()->normalizedChi2()) ? 0 : cand.innerTrack()->normalizedChi2());
		}
		else{
			fts.fill("jet_pfcand_muon_nvalidhits",0.);
			fts.fill("jet_pfcand_muon_chi2",0.);
		}

		fts.fill("jet_pfcand_muon_nstations",isNanOrInf(cand.numberOfMatchedStations()) ? 0 : cand.numberOfMatchedStations());
		fts.fill("jet_pfcand_muon_segcomp",isNanOrInf(muon::segmentCompatibility(cand)) ? 0 : muon::segmentCompatibility(cand));
		fts.fill("jet_pfcand_muon_trkKink",isNanOrInf(cand.combinedQuality().trkKink) ? 0 : cand.combinedQuality().trkKink);
		fts.fill("jet_pfcand_muon_validfraction",isNanOrInf((cand.innerTrack().isNonnull()) ? cand.innerTrack()->validFraction() : 0) ? 0 :
				((cand.innerTrack().isNonnull()) ? cand.innerTrack()->validFraction() : 0));
		fts.fill("jet_pfcand_muon_dxy",ip_sign_*(isNanOrInf((cand.muonBestTrack().isNonnull()) ? cand.muonBestTrack()->dxy(pv_ass_pos) : 0) ? 0 :
					((cand.muonBestTrack().isNonnull()) ? cand.muonBestTrack()->dxy(pv_ass_pos) : 0))); 
		fts.fill("jet_pfcand_muon_dz",ip_sign_*(isNanOrInf((cand.muonBestTrack().isNonnull()) ? cand.muonBestTrack()->dz(pv_ass_pos) : 0) ? 0 :
					((cand.muonBestTrack().isNonnull()) ? cand.muonBestTrack()->dz(pv_ass_pos) : 0))); 
		fts.fill("jet_pfcand_muon_pterr",isNanOrInf((cand.tunePMuonBestTrack().isNonnull()) ? cand.tunePMuonBestTrack()->ptError()/cand.tunePMuonBestTrack()->pt() : 0) ? 0 :
				((cand.tunePMuonBestTrack().isNonnull()) ? cand.tunePMuonBestTrack()->ptError()/cand.tunePMuonBestTrack()->pt() : 0));
		fts.fill("jet_pfcand_muon_nlayers",isNanOrInf((cand.innerTrack().isNonnull()) ? cand.innerTrack()->hitPattern().stripLayersWithMeasurement()+cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0) ? 0 : ((cand.innerTrack().isNonnull()) ? cand.innerTrack()->hitPattern().stripLayersWithMeasurement()+cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0));
		fts.fill("jet_pfcand_muon_nhits",isNanOrInf((cand.innerTrack().isNonnull()) ? cand.innerTrack()->hitPattern().numberOfValidStripHits()+cand.innerTrack()->hitPattern().numberOfValidPixelHits() : 0) ? 0 : ((cand.innerTrack().isNonnull()) ? cand.innerTrack()->hitPattern().numberOfValidStripHits()+cand.innerTrack()->hitPattern().numberOfValidPixelHits() : 0));
		fts.fill("jet_pfcand_muon_mask",1);
	}

	// Build observables electron info
	for (const auto & icand : daughters_electron_index){
		auto cand = electrons[icand];
		fts.fill("jet_pfcand_electron_pt", isNanOrInf(cand.pt()) ? 0 : cand.pt());
		fts.fill("jet_pfcand_electron_pt_log", isNanOrInf(std::log(cand.pt())) ? 0 : std::log(cand.pt()));
		fts.fill("jet_pfcand_electron_eta", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_electron_eta_v", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_electron_phi", isNanOrInf(cand.phi()) ? 0 : cand.phi());
		fts.fill("jet_pfcand_electron_energy",isNanOrInf(cand.energy()) ? 0 : cand.energy());
		fts.fill("jet_pfcand_electron_energy_log", isNanOrInf(std::log(cand.energy())) ? 0 : std::log(cand.energy()));
		fts.fill("jet_pfcand_electron_deta",isNanOrInf(cand.eta()-jet.eta()) ? 0 : cand.eta()-jet.eta());
		fts.fill("jet_pfcand_electron_dphi",isNanOrInf(reco::deltaPhi(cand.phi(),jet.phi())) ? 0 : reco::deltaPhi(cand.phi(),jet.phi()));
		fts.fill("jet_pfcand_electron_dxy",ip_sign_*(isNanOrInf((cand.gsfTrack().isNonnull()) ? cand.gsfTrack()->dxy(pv_ass_pos) : 0) ? 0 : ((cand.gsfTrack().isNonnull()) ? cand.gsfTrack()->dxy(pv_ass_pos) : 0)));
		fts.fill("jet_pfcand_electron_dz",ip_sign_*(isNanOrInf((cand.gsfTrack().isNonnull()) ? cand.gsfTrack()->dz(pv_ass_pos) : 0) ? 0 : ((cand.gsfTrack().isNonnull()) ? cand.gsfTrack()->dz(pv_ass_pos) : 0)));
		fts.fill("jet_pfcand_electron_eOverP",isNanOrInf(cand.eSuperClusterOverP()) ? 0 : cand.eSuperClusterOverP());
		fts.fill("jet_pfcand_electron_detaIn",isNanOrInf(cand.deltaEtaSuperClusterTrackAtVtx()) ? 0 : cand.deltaEtaSuperClusterTrackAtVtx());
		fts.fill("jet_pfcand_electron_dphiIn",isNanOrInf(cand.deltaPhiSuperClusterTrackAtVtx()) ? 0 : cand.deltaPhiSuperClusterTrackAtVtx());
		fts.fill("jet_pfcand_electron_hOverE",isNanOrInf(cand.hadronicOverEm()) ? 0 : cand.hadronicOverEm());
		fts.fill("jet_pfcand_electron_sigIetaIeta",isNanOrInf(cand.full5x5_sigmaIetaIeta()) ? 0 : cand.full5x5_sigmaIetaIeta());
		fts.fill("jet_pfcand_electron_r9",isNanOrInf(cand.full5x5_r9()) ? 0 : cand.full5x5_r9());
		fts.fill("jet_pfcand_electron_convProb",isNanOrInf(cand.convVtxFitProb()) ? 0 : cand.convVtxFitProb());
		fts.fill("jet_pfcand_electron_mask",1);
	}

	// Build observables photon info
	for (const auto & icand : daughters_photon_index){
		auto cand = photons[icand];
		fts.fill("jet_pfcand_photon_pt", isNanOrInf(cand.pt()) ? 0 : cand.pt());
		fts.fill("jet_pfcand_photon_pt_log", isNanOrInf(std::log(cand.pt())) ? 0 : std::log(cand.pt()));
		fts.fill("jet_pfcand_photon_eta", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_photon_eta_v", isNanOrInf(cand.eta()) ? 0 : cand.eta());
		fts.fill("jet_pfcand_photon_phi", isNanOrInf(cand.phi()) ? 0 : cand.phi());
		fts.fill("jet_pfcand_photon_energy",isNanOrInf(cand.energy()) ? 0 : cand.energy());
		fts.fill("jet_pfcand_photon_energy_log", isNanOrInf(std::log(cand.energy())) ? 0 : std::log(cand.energy()));
		fts.fill("jet_pfcand_photon_deta",isNanOrInf(cand.eta()-jet.eta()) ? 0 : cand.eta()-jet.eta());
		fts.fill("jet_pfcand_photon_dphi",isNanOrInf(reco::deltaPhi(cand.phi(),jet.phi())) ? 0 : reco::deltaPhi(cand.phi(),jet.phi()));
		fts.fill("jet_pfcand_photon_hOverE",isNanOrInf(cand.hadronicOverEm()) ? 0 : cand.hadronicOverEm());
		fts.fill("jet_pfcand_photon_sigIetaIeta",isNanOrInf(cand.full5x5_sigmaIetaIeta()) ? 0 : cand.full5x5_sigmaIetaIeta());
		fts.fill("jet_pfcand_photon_r9",isNanOrInf(cand.full5x5_r9()) ? 0 : cand.full5x5_r9());
		fts.fill("jet_pfcand_photon_eVeto",isNanOrInf(cand.passElectronVeto()) ? 0 : cand.passElectronVeto());
		fts.fill("jet_pfcand_photon_mask",1);
	}
}

void ParTFeatureEvaluator::fillSVFeatures(DeepBoostedJetFeatures &fts, 
		const reco::Jet & jet) {

	// secondary vertexes matching jet
	std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;
	for (const auto &sv : *svs_) {
		if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
			jetSVs.push_back(&sv);
		}
	}

	// sort by dxy significance
	std::sort(jetSVs.begin(),
			jetSVs.end(),
			[&](const reco::VertexCompositePtrCandidate *sva, 
				const reco::VertexCompositePtrCandidate *svb) {
			return sv_vertex_comparator(*sva, *svb, *pv_);
			});

	// reserve space
	for (const auto &name : sv_features_)
		fts.reserve(name, jetSVs.size());

	GlobalVector jet_global_vec(jet.px(), jet.py(), jet.pz());

	for (const auto *sv : jetSVs) {
		fts.fill("jet_sv_pt", isNanOrInf(sv->pt()) ? 0 : sv->pt());
		fts.fill("jet_sv_pt_log", isNanOrInf(std::log(sv->pt())) ? 0 : std::log(sv->pt()));
		fts.fill("jet_sv_eta", isNanOrInf(sv->eta()) ? 0 : sv->eta());
		fts.fill("jet_sv_eta_v", isNanOrInf(sv->eta()) ? 0 : sv->eta());
		fts.fill("jet_sv_phi", isNanOrInf(sv->phi()) ? 0 : sv->phi());
		fts.fill("jet_sv_energy", isNanOrInf(sv->energy()) ? 0 : sv->energy());
		fts.fill("jet_sv_energy_log", isNanOrInf(std::log(sv->energy())) ? 0 : std::log(sv->energy()));
		fts.fill("jet_sv_mass", isNanOrInf(sv->mass()) ? 0 : sv->mass());
		fts.fill("jet_sv_deta", isNanOrInf(sv->eta() - jet.eta()) ? 0 : sv->eta() - jet.eta());
		fts.fill("jet_sv_dphi", isNanOrInf(reco::deltaPhi(sv->phi(),jet.phi())) ? 0 : reco::deltaPhi(sv->phi(),jet.phi()));
		fts.fill("jet_sv_chi2", isNanOrInf(sv->vertexNormalizedChi2()) ? 0 : sv->vertexNormalizedChi2());
		fts.fill("jet_sv_ntrack", isNanOrInf(sv->numberOfDaughters()) ? 0 : sv->numberOfDaughters());

		reco::Candidate::Vector p = sv->momentum();
		reco::Candidate::Vector d(sv->vx() - pv_->x(), sv->vy() - pv_->y(), sv->vz() - pv_->z());
		auto costheta = p.Unit().Dot(d.Unit());
		fts.fill("jet_sv_costheta", isNanOrInf(costheta) ? 0 : costheta);

		reco::Vertex::CovarianceMatrix csv;
		sv->fillVertexCovariance(csv);
		reco::Vertex svtx(sv->vertex(), csv);    
		VertexDistanceXY dxy;
		auto valxy = dxy.distance(svtx, *pv_);
		VertexDistance3D d3d;
		auto val3d = d3d.distance(svtx, *pv_);    
		auto val3d_signed = d3d.signedDistance(svtx, *pv_, jet_global_vec);

		fts.fill("jet_sv_dxy", isNanOrInf(valxy.value()) ? 0 : valxy.value());
		fts.fill("jet_sv_dxysig", isNanOrInf(valxy.significance()) ? 0 : valxy.significance());    
		fts.fill("jet_sv_d3d", isNanOrInf(val3d.value()) ? 0 : val3d.value());
		fts.fill("jet_sv_d3dsig", isNanOrInf(val3d.significance()) ? 0 : val3d.significance());
		fts.fill("jet_sv_d3d_sign", ip_sign_*(isNanOrInf(val3d_signed.value()/val3d.value()) ? 0 : val3d_signed.value()/val3d.value()));
		fts.fill("jet_sv_mask", 1);
	}
}


void ParTFeatureEvaluator::fillKVFeatures(DeepBoostedJetFeatures &fts, 
		const reco::Jet & jet) {

	// secondary vertexes matching jet
	std::vector<const reco::VertexCompositePtrCandidate*> jetKVs;
	for (const auto &kv : *kvs_) {
		if (reco::deltaR2(kv, jet) < jet_radius_ * jet_radius_) {
			jetKVs.push_back(&kv);
		}
	}

	// sort by dxy significance
	std::sort(jetKVs.begin(),
			jetKVs.end(),
			[&](const reco::VertexCompositePtrCandidate *kva, 
				const reco::VertexCompositePtrCandidate *kvb) {
			return sv_vertex_comparator(*kva, *kvb, *pv_);
			});

	// reserve space
	for (const auto &name : kv_features_)
		fts.reserve(name, jetKVs.size());

	GlobalVector jet_global_vec(jet.px(), jet.py(), jet.pz());

	for (const auto *kv : jetKVs) {
		fts.fill("jet_kaon_pt", isNanOrInf(kv->pt()) ? 0 : kv->pt());
		fts.fill("jet_kaon_pt_log", isNanOrInf(std::log(kv->pt())) ? 0 : std::log(kv->pt()));
		fts.fill("jet_kaon_eta", isNanOrInf(kv->eta()) ? 0 : kv->eta());
		fts.fill("jet_kaon_eta_v", isNanOrInf(kv->eta()) ? 0 : kv->eta());
		fts.fill("jet_kaon_phi", isNanOrInf(kv->phi()) ? 0 : kv->phi());
		fts.fill("jet_kaon_energy", isNanOrInf(kv->energy()) ? 0 : kv->energy());
		fts.fill("jet_kaon_energy_log", isNanOrInf(std::log(kv->energy())) ? 0 : std::log(kv->energy()));
		fts.fill("jet_kaon_mass", isNanOrInf(kv->mass()) ? 0 : kv->mass());
		fts.fill("jet_kaon_deta", isNanOrInf(kv->eta() - jet.eta()) ? 0 : kv->eta() - jet.eta());
		fts.fill("jet_kaon_dphi", isNanOrInf(reco::deltaPhi(kv->phi(),jet.phi())) ? 0 : reco::deltaPhi(kv->phi(),jet.phi()));
		fts.fill("jet_kaon_chi2", isNanOrInf(kv->vertexNormalizedChi2()) ? 0 : kv->vertexNormalizedChi2());
		fts.fill("jet_kaon_ntrack", isNanOrInf(kv->numberOfDaughters()) ? 0 : kv->numberOfDaughters());

		reco::Candidate::Vector p = kv->momentum();
		reco::Candidate::Vector d(kv->vx() - pv_->x(), kv->vy() - pv_->y(), kv->vz() - pv_->z());
		auto costheta = p.Unit().Dot(d.Unit());
		fts.fill("jet_kaon_costheta", isNanOrInf(costheta) ? 0 : costheta);

		reco::Vertex::CovarianceMatrix csv;
		kv->fillVertexCovariance(csv);
		reco::Vertex kvtx(kv->vertex(), csv);    
		VertexDistanceXY dxy;
		auto valxy = dxy.distance(kvtx, *pv_);
		VertexDistance3D d3d;
		auto val3d = d3d.distance(kvtx, *pv_);    
		auto val3d_signed = d3d.signedDistance(kvtx, *pv_, jet_global_vec);

		fts.fill("jet_kaon_dxy", isNanOrInf(valxy.value()) ? 0 : valxy.value());
		fts.fill("jet_kaon_dxysig", isNanOrInf(valxy.significance()) ? 0 : valxy.significance());    
		fts.fill("jet_kaon_d3d", isNanOrInf(val3d.value()) ? 0 : val3d.value());
		fts.fill("jet_kaon_d3dsig", isNanOrInf(val3d.significance()) ? 0 : val3d.significance());
		fts.fill("jet_kaon_d3d_sign", ip_sign_*(isNanOrInf(val3d_signed.value()/val3d.value()) ? 0 : val3d_signed.value()/val3d.value()));
		fts.fill("jet_kaon_mask", 1);
	}
}


void ParTFeatureEvaluator::fillLVFeatures(DeepBoostedJetFeatures &fts, 
		const reco::Jet & jet) {

	// secondary vertexes matching jet
	std::vector<const reco::VertexCompositePtrCandidate*> jetLVs;
	for (const auto &lv : *lvs_) {
		if (reco::deltaR2(lv, jet) < jet_radius_ * jet_radius_) {
			jetLVs.push_back(&lv);
		}
	}

	// sort by dxy significance
	std::sort(jetLVs.begin(),
			jetLVs.end(),
			[&](const reco::VertexCompositePtrCandidate *lva, 
				const reco::VertexCompositePtrCandidate *lvb) {
			return sv_vertex_comparator(*lva, *lvb, *pv_);
			});

	// reserve space
	for (const auto &name : lv_features_)
		fts.reserve(name, jetLVs.size());

	GlobalVector jet_global_vec(jet.px(), jet.py(), jet.pz());

	for (const auto *lv : jetLVs) {
		fts.fill("jet_lambda_pt", isNanOrInf(lv->pt()) ? 0 : lv->pt());
		fts.fill("jet_lambda_pt_log", isNanOrInf(std::log(lv->pt())) ? 0 : std::log(lv->pt()));
		fts.fill("jet_lambda_eta", isNanOrInf(lv->eta()) ? 0 : lv->eta());
		fts.fill("jet_lambda_eta_v", isNanOrInf(lv->eta()) ? 0 : lv->eta());
		fts.fill("jet_lambda_phi", isNanOrInf(lv->phi()) ? 0 : lv->phi());
		fts.fill("jet_lambda_energy", isNanOrInf(lv->energy()) ? 0 : lv->energy());
		fts.fill("jet_lambda_energy_log", isNanOrInf(std::log(lv->energy())) ? 0 : std::log(lv->energy()));
		fts.fill("jet_lambda_mass", isNanOrInf(lv->mass()) ? 0 : lv->mass());
		fts.fill("jet_lambda_deta", isNanOrInf(lv->eta() - jet.eta()) ? 0 : lv->eta() - jet.eta());
		fts.fill("jet_lambda_dphi", isNanOrInf(reco::deltaPhi(lv->phi(),jet.phi())) ? 0 : reco::deltaPhi(lv->phi(),jet.phi()));
		fts.fill("jet_lambda_chi2", isNanOrInf(lv->vertexNormalizedChi2()) ? 0 : lv->vertexNormalizedChi2());
		fts.fill("jet_lambda_ntrack", isNanOrInf(lv->numberOfDaughters()) ? 0 : lv->numberOfDaughters());

		reco::Candidate::Vector p = lv->momentum();
		reco::Candidate::Vector d(lv->vx() - pv_->x(), lv->vy() - pv_->y(), lv->vz() - pv_->z());
		auto costheta = p.Unit().Dot(d.Unit());
		fts.fill("jet_lambda_costheta", isNanOrInf(costheta) ? 0 : costheta);

		reco::Vertex::CovarianceMatrix csv;
		lv->fillVertexCovariance(csv);
		reco::Vertex lvtx(lv->vertex(), csv);    
		VertexDistanceXY dxy;
		auto valxy = dxy.distance(lvtx, *pv_);
		VertexDistance3D d3d;
		auto val3d = d3d.distance(lvtx, *pv_);    
		auto val3d_signed = d3d.signedDistance(lvtx, *pv_, jet_global_vec);

		fts.fill("jet_lambda_dxy", isNanOrInf(valxy.value()) ? 0 : valxy.value());
		fts.fill("jet_lambda_dxysig", isNanOrInf(valxy.significance()) ? 0 : valxy.significance());    
		fts.fill("jet_lambda_d3d", isNanOrInf(val3d.value()) ? 0 : val3d.value());
		fts.fill("jet_lambda_d3dsig", isNanOrInf(val3d.significance()) ? 0 : val3d.significance());
		fts.fill("jet_lambda_d3d_sign", ip_sign_*(isNanOrInf(val3d_signed.value()/val3d.value()) ? 0 : val3d_signed.value()/val3d.value()));
		fts.fill("jet_lambda_mask", 1);
	}
}

void ParTFeatureEvaluator::fillLostTrackFeatures(DeepBoostedJetFeatures &fts, 
		const reco::Jet & jet
		) {

	// some jet properties
	TVector3        jet_direction (jet.momentum().Unit().x(), jet.momentum().Unit().y(), jet.momentum().Unit().z());
	GlobalVector    jet_ref_track_dir (jet.px(), jet.py(), jet.pz());
	math::XYZVector jet_dir = jet.momentum().Unit();

	std::vector<pat::PackedCandidate> jet_lost_tracks;
	for(size_t itrk = 0; itrk < losttracks_->size(); itrk++){
		if(reco::deltaR(losttracks_->at(itrk).p4(),jet.p4()) < max_dr_for_losttrack_ and
				losttracks_->at(itrk).pt() > min_pt_for_losttrack_ ){
			jet_lost_tracks.push_back(losttracks_->at(itrk));
		}
	}
	std::sort(jet_lost_tracks.begin(), jet_lost_tracks.end(), [](const auto &a, const auto &b) { return a.pt() > b.pt();});

	// reserve space
	for (const auto &name : losttrack_features_)
		fts.reserve(name, jet_lost_tracks.size());

	reco::VertexRef pv_ass = reco::VertexRef(vtxs_, 0);
	math::XYZPoint  pv_ass_pos = pv_ass->position();

	for(auto const & ltrack : jet_lost_tracks){

		fts.fill("jet_losttrack_pt",isNanOrInf(ltrack.pt()) ? 0 : ltrack.pt());
		fts.fill("jet_losttrack_pt_log",isNanOrInf(std::log(ltrack.pt())) ? 0 : std::log(ltrack.pt()));
		fts.fill("jet_losttrack_energy",isNanOrInf(ltrack.energy()) ? 0 : ltrack.energy());
		fts.fill("jet_losttrack_energy_log",isNanOrInf(std::log(ltrack.energy())) ? 0 : std::log(ltrack.energy()));
		fts.fill("jet_losttrack_eta",isNanOrInf(ltrack.eta()) ? 0 : ltrack.eta());
		fts.fill("jet_losttrack_eta_v",isNanOrInf(ltrack.eta()) ? 0 : ltrack.eta());
		fts.fill("jet_losttrack_phi",isNanOrInf(ltrack.phi()) ? 0 : ltrack.phi());    
		fts.fill("jet_losttrack_charge",isNanOrInf(ltrack.charge()) ? 0 : ltrack.charge());
		fts.fill("jet_losttrack_frompv",isNanOrInf(ltrack.fromPV()) ? 0 : ltrack.fromPV());
		fts.fill("jet_losttrack_dz",ip_sign_*(isNanOrInf(ltrack.dz(pv_ass_pos)) ? 0 : ltrack.dz(pv_ass_pos)));
		fts.fill("jet_losttrack_dxy",ip_sign_*(isNanOrInf(ltrack.dxy(pv_ass_pos)) ? 0 : ltrack.dxy(pv_ass_pos)));
		TVector3 ltrack_momentum (ltrack.momentum().x(),ltrack.momentum().y(),ltrack.momentum().z());
		fts.fill("jet_losttrack_deta",isNanOrInf(jet_direction.Eta()-ltrack_momentum.Eta()) ? 0 : jet_direction.Eta()-ltrack_momentum.Eta());
		fts.fill("jet_losttrack_dphi",isNanOrInf(jet_direction.DeltaPhi(ltrack_momentum)) ? 0 : jet_direction.DeltaPhi(ltrack_momentum));
		fts.fill("jet_losttrack_etarel",isNanOrInf(reco::btau::etaRel(jet_dir,ltrack.momentum())) ? 0 : reco::btau::etaRel(jet_dir,ltrack.momentum()));
		fts.fill("jet_losttrack_ptrel",isNanOrInf(ltrack_momentum.Perp(jet_direction)) ? 0 : ltrack_momentum.Perp(jet_direction));
		fts.fill("jet_losttrack_nlosthits",isNanOrInf(ltrack.lostInnerHits()) ? 0 : ltrack.lostInnerHits());

		const reco::Track* track = ltrack.bestTrack();
		if(track){
			fts.fill("jet_losttrack_dxysig",ip_sign_*(isNanOrInf(track->dxy(pv_ass_pos)/track->dxyError()) ? 0 : track->dxy(pv_ass_pos)/track->dxyError()));
			fts.fill("jet_losttrack_dzsig",ip_sign_*(isNanOrInf(track->dz(pv_ass_pos)/track->dzError()) ? 0 : track->dz(pv_ass_pos)/track->dzError()));
			fts.fill("jet_losttrack_track_chi2",isNanOrInf(track->normalizedChi2()) ? 0 : track->normalizedChi2());
			fts.fill("jet_losttrack_track_qual",isNanOrInf(track->qualityMask()) ? 0 : track->qualityMask());
			fts.fill("jet_losttrack_track_pterr",isNanOrInf(track->ptError()/track->pt()) ? 0 : track->ptError()/track->pt());
			fts.fill("jet_losttrack_track_etaerr",isNanOrInf(track->etaError()) ? 0 : track->etaError());
			fts.fill("jet_losttrack_track_phierr",isNanOrInf(track->phiError()) ? 0 : track->phiError());
			fts.fill("jet_losttrack_nhits",isNanOrInf(track->hitPattern().numberOfValidHits()) ? 0 : track->hitPattern().numberOfValidHits());
			fts.fill("jet_losttrack_npixhits",isNanOrInf(track->hitPattern().numberOfValidPixelHits()) ? 0 : track->hitPattern().numberOfValidPixelHits());
			fts.fill("jet_losttrack_npixbarrelhits",isNanOrInf(track->hitPattern().numberOfValidPixelBarrelHits()) ? 0 : track->hitPattern().numberOfValidPixelBarrelHits());
			fts.fill("jet_losttrack_nstriphits",isNanOrInf(track->hitPattern().numberOfValidStripHits()) ? 0 : track->hitPattern().numberOfValidStripHits());
			fts.fill("jet_losttrack_nstriptibhits",isNanOrInf(track->hitPattern().numberOfValidStripTIBHits()+track->hitPattern().numberOfValidStripTIDHits()) ? 0 : track->hitPattern().numberOfValidStripTIBHits()+track->hitPattern().numberOfValidStripTIDHits());
			fts.fill("jet_losttrack_nstriptobhits",isNanOrInf(track->hitPattern().numberOfValidStripTOBHits()) ? 0 : track->hitPattern().numberOfValidStripTOBHits());
			fts.fill("jet_losttrack_nlayers",isNanOrInf(track->hitPattern().trackerLayersWithMeasurement()) ? 0 : track->hitPattern().trackerLayersWithMeasurement());
			fts.fill("jet_losttrack_npixlayers",isNanOrInf(track->hitPattern().pixelLayersWithMeasurement()) ? 0 : track->hitPattern().pixelLayersWithMeasurement());
			fts.fill("jet_losttrack_nstriplayers",isNanOrInf(track->hitPattern().stripLayersWithMeasurement()) ? 0 : track->hitPattern().stripLayersWithMeasurement());

			reco::TransientTrack transientTrack = track_builder_->build(*track);
			Measurement1D meas_ip2d    = IPTools::signedTransverseImpactParameter(transientTrack, jet_ref_track_dir, *pv_).second;
			Measurement1D meas_ip3d    = IPTools::signedImpactParameter3D(transientTrack, jet_ref_track_dir, *pv_).second;
			Measurement1D meas_jetdist = IPTools::jetTrackDistance(transientTrack, jet_ref_track_dir, *pv_).second;
			Measurement1D meas_decayl  = IPTools::signedDecayLength3D(transientTrack, jet_ref_track_dir, *pv_).second;

			fts.fill("jet_losttrack_trackjet_sip2d",ip_sign_*(isNanOrInf(meas_ip2d.value()) ? 0 : meas_ip2d.value()));
			fts.fill("jet_losttrack_trackjet_sip2dsig",ip_sign_*(isNanOrInf(meas_ip2d.significance()) ? 0 : meas_ip2d.significance()));
			fts.fill("jet_losttrack_trackjet_sip3d",ip_sign_*(isNanOrInf(meas_ip3d.value()) ? 0 : meas_ip3d.value()));
			fts.fill("jet_losttrack_trackjet_sip3dsig",ip_sign_*(isNanOrInf(meas_ip3d.significance()) ? 0 : meas_ip3d.significance()));
			fts.fill("jet_losttrack_trackjet_dist",isNanOrInf(-meas_jetdist.value()) ? 0 : -meas_jetdist.value());
			fts.fill("jet_losttrack_trackjet_decayL",isNanOrInf(meas_decayl.value()) ? 0 : meas_decayl.value());
		}
		else{
			fts.fill("jet_losttrack_dxysig",0);
			fts.fill("jet_losttrack_dzsig",0);
			fts.fill("jet_losttrack_track_chi2",0);
			fts.fill("jet_losttrack_track_qual",0);
			fts.fill("jet_losttrack_track_pterr",0);
			fts.fill("jet_losttrack_track_etaerr",0);
			fts.fill("jet_losttrack_track_phierr",0);
			fts.fill("jet_losttrack_nhits",0);
			fts.fill("jet_losttrack_npixhits",0);
			fts.fill("jet_losttrack_npixbarrelhits",0);
			fts.fill("jet_losttrack_nstriphits",0);
			fts.fill("jet_losttrack_nstriptibhits",0);
			fts.fill("jet_losttrack_nstriptobhits",0);
			fts.fill("jet_losttrack_nlayers",0);
			fts.fill("jet_losttrack_npixlayers",0);
			fts.fill("jet_losttrack_trackjet_sip2d",0);
			fts.fill("jet_losttrack_trackjet_sip2dsig",0);
			fts.fill("jet_losttrack_trackjet_sip3d",0);
			fts.fill("jet_losttrack_trackjet_sip3dsig",0);
			fts.fill("jet_losttrack_trackjet_dist",0);
			fts.fill("jet_losttrack_trackjet_decayL",0);
		}
		fts.fill("jet_losttrack_mask",1);    
	}
}
// define this as a plug-in
DEFINE_FWK_MODULE(ParTFeatureEvaluator);

