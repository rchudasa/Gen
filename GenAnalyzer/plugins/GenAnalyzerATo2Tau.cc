// -*- C++ -*-
//
// Package: Gen/GenAnalyzerATo2Tau
// Class: GenAnalyzerATo2Tau
//

#include "Gen/GenAnalyzer/interface/GenAnalyzerATo2Tau.h"


using reco::GenParticle;

// Event info
unsigned int run_;
unsigned int lumi_;
unsigned long long event_;

// Kinematics
std::vector<float> v_a_pt_, v_a_eta_, v_a_phi_, v_a_mass_, v_a_mass_vis_;
std::vector<float> v_tau1_pt_, v_tau1_eta_, v_tau1_phi_;
std::vector<float> v_tau2_pt_, v_tau2_eta_, v_tau2_phi_;
std::vector<float> v_dR_tautau_, v_dEta_tautau_, v_dPhi_tautau_;

// Counters
int nTotalEvents_;
int nPassedEvents_;

TLorentzVector SetTs(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}


GenAnalyzerATo2Tau::GenAnalyzerATo2Tau(const edm::ParameterSet& iConfig)
//  : 
{
  isDebug  = iConfig.getParameter<bool>("isDebug");
  genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
   
  tree_ = fs->make<TTree>("RHTree", "a -> tautau Generator Information");

  // Event identification
  tree_->Branch("run",    &run_,    "run/i");
  tree_->Branch("lumi",   &lumi_,   "lumi/i");
  tree_->Branch("event",  &event_,  "event/l");

  // Vector branches for multiple 'a' per event
  tree_->Branch("a_pt",       &v_a_pt_);
  tree_->Branch("a_eta",      &v_a_eta_);
  tree_->Branch("a_phi",      &v_a_phi_);
  tree_->Branch("a_mass",     &v_a_mass_);       // truth mass
  tree_->Branch("a_mass_vis", &v_a_mass_vis_);   // reconstructed from tau pair

  tree_->Branch("tau1_pt",    &v_tau1_pt_);
  tree_->Branch("tau1_eta",   &v_tau1_eta_);
  tree_->Branch("tau1_phi",   &v_tau1_phi_);

  tree_->Branch("tau2_pt",    &v_tau2_pt_);
  tree_->Branch("tau2_eta",   &v_tau2_eta_);
  tree_->Branch("tau2_phi",   &v_tau2_phi_);

  tree_->Branch("dR_tautau",   &v_dR_tautau_);
  tree_->Branch("dEta_tautau", &v_dEta_tautau_);
  tree_->Branch("dPhi_tautau", &v_dPhi_tautau_);
}

GenAnalyzerATo2Tau::~GenAnalyzerATo2Tau(){

}

void GenAnalyzerATo2Tau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  run_   = iEvent.id().run();
  lumi_  = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  // Clear all vectors
  v_a_pt_.clear();       v_a_eta_.clear();       v_a_phi_.clear();
  v_a_mass_.clear();     v_a_mass_vis_.clear();
  v_tau1_pt_.clear();    v_tau1_eta_.clear();    v_tau1_phi_.clear();
  v_tau2_pt_.clear();    v_tau2_eta_.clear();    v_tau2_phi_.clear();
  v_dR_tautau_.clear();  v_dEta_tautau_.clear(); v_dPhi_tautau_.clear();

  bool found = false;

  Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  for (auto iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {

    // Look for pseudoscalar 'a' decaying to two taus
    if (std::abs(iGen->pdgId()) != std::abs(25)) continue;
    if (iGen->numberOfDaughters() != 2) continue;
    std::cout << " apdgId:"<< iGen->pdgId() << " pt:" << iGen->pt() << " mass:" << iGen->mass() ;
 std::cout << "pdg id:" << iGen->daughter(0)->pdgId() << " ," << iGen->daughter(0)->status() << " dau 2: " << iGen->daughter(1)->pdgId() << " , " << iGen->daughter(1)->status()<< std::endl;
    const reco::Candidate* dau0 = iGen->daughter(0);
    const reco::Candidate* dau1 = iGen->daughter(1);
    if (!dau0 || !dau1) continue;

    if (std::abs(dau0->pdgId()) != 15 || std::abs(dau1->pdgId()) != 15) continue;
    if (dau0->status() != 2 || dau1->status() != 2) continue;

    // Build visible di-tau system
    TLorentzVector tau1 = SetTs(dau0->pt(), dau0->eta(), dau0->phi(), dau0->mass());
    TLorentzVector tau2 = SetTs(dau1->pt(), dau1->eta(), dau1->phi(), dau1->mass());
    TLorentzVector a_vis = tau1 + tau2;

    // Order taus: tau1 has higher pT
    if (dau0->pt() < dau1->pt()) std::swap(tau1, tau2);
    std::cout << "Tau1 pt:" << tau1.Pt() << " Tau2 pt:" << tau2.Pt() << std::endl;
    // Fill vectors (exactly like your original code)
    v_a_pt_.push_back(iGen->pt());
    v_a_eta_.push_back(iGen->eta());
    v_a_phi_.push_back(iGen->phi());
    v_a_mass_.push_back(iGen->mass());
    v_a_mass_vis_.push_back(a_vis.M());

    v_tau1_pt_.push_back(tau1.Pt());
    v_tau1_eta_.push_back(tau1.Eta());
    v_tau1_phi_.push_back(tau1.Phi());

    v_tau2_pt_.push_back(tau2.Pt());
    v_tau2_eta_.push_back(tau2.Eta());
    v_tau2_phi_.push_back(tau2.Phi());

    v_dR_tautau_.push_back(deltaR(*dau0, *dau1));
    v_dEta_tautau_.push_back(std::abs(dau0->eta() - dau1->eta()));
    v_dPhi_tautau_.push_back(deltaPhi(dau0->phi(), dau1->phi()));
  }

  tree_->Fill();
  if (!v_a_pt_.empty()) nPassedEvents_++;
  nTotalEvents_++;

}

void GenAnalyzerATo2Tau::beginJob() {
  nTotalEvents_ = nPassedEvents_ = 0;
}

void GenAnalyzerATo2Tau::endJob() {
  std::cout << "\n=== GenAnalyzerATo2Tau Summary ===\n";
  std::cout << "Total events processed      : " << nTotalEvents_ << "\n";
  std::cout << "Events with a -> tautau     : " << nPassedEvents_ << "\n";
  std::cout << "Efficiency                  : "
            << (nTotalEvents_ > 0 ? 100.0 * nPassedEvents_ / nTotalEvents_ : 0)
            << "%\n\n";
}

void GenAnalyzerATo2Tau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenAnalyzerATo2Tau);
