// -*- C++ -*-
//
// Package:     Gen/GenAnalyzerAee
// Class:       GenAnalyzerAee
//
/**\class GenAnalyzerAee GenAnalyzerAee.cc Gen/GenAnalyzer/plugins/GenAnalyzerAee.cc
   Description: Simple generator-level analyzer for A → e⁺ e⁻ (A = pseudoscalar, PDG ID ±25)
*/
//
// Last updated: for single A → ee decay (PDG 25 / -25)
//
#include "Gen/GenAnalyzer/interface/GenAnalyzerAee.h"   // ← update header name if needed
#include <iomanip>

using reco::GenParticle;

int nTotalAee  = 0;
int nPassedAee = 0;

unsigned int runIDAee_;
unsigned int lumiIDAee_;
unsigned long long eventIDAee_;

// ─── Single A → ee quantities ──────────────────────────────────────
vector<float> V_ae_genA_M_inv_;       // invariant mass from e⁺e⁻
vector<float> V_ae_genA_M_;           // generator mass of A
vector<float> V_ae_dR_Ele1_Ele2_;     // ΔR between the two electrons
vector<float> V_ae_A_pt_;
vector<float> V_ae_A_eta_;
vector<float> V_ae_A_phi_;
vector<float> V_ae_Ele1_pt_;
vector<float> V_ae_Ele1_eta_;
vector<float> V_ae_Ele1_phi_;
vector<float> V_ae_Ele2_pt_;
vector<float> V_ae_Ele2_eta_;
vector<float> V_ae_Ele2_phi_;
vector<float> V_ae_Ele1_Ele2_deta_;
vector<float> V_ae_Ele1_Ele2_dphi_;


TLorentzVector SetFourVector(float pt, float eta, float phi, float mass) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, mass);
    return v;
}


GenAnalyzerAee::GenAnalyzerAee(const edm::ParameterSet& iConfig) {

    isDebug     = iConfig.getParameter<bool>("isDebug");
    print_trigger = iConfig.getParameter<bool>("print_trigger");

    RHTree = fs->make<TTree>("AeeTree", "Gen level info - A → ee (A = 25/-25)");

    RHTree->Branch("run",   &runIDAee_);
    RHTree->Branch("lumi",  &lumiIDAee_);
    RHTree->Branch("event", &eventIDAee_);

    RHTree->Branch("GenA_inv",     &V_ae_genA_M_inv_);     // from daughters
    RHTree->Branch("GenA",         &V_ae_genA_M_);         // from gen particle
    RHTree->Branch("dR_Ele1_Ele2", &V_ae_dR_Ele1_Ele2_);

    RHTree->Branch("A_pt",   &V_ae_A_pt_);
    RHTree->Branch("A_eta",  &V_ae_A_eta_);
    RHTree->Branch("A_phi",  &V_ae_A_phi_);

    RHTree->Branch("Ele1_pt",  &V_ae_Ele1_pt_);
    RHTree->Branch("Ele1_eta", &V_ae_Ele1_eta_);
    RHTree->Branch("Ele1_phi", &V_ae_Ele1_phi_);
    RHTree->Branch("Ele2_pt",  &V_ae_Ele2_pt_);
    RHTree->Branch("Ele2_eta", &V_ae_Ele2_eta_);
    RHTree->Branch("Ele2_phi", &V_ae_Ele2_phi_);

    RHTree->Branch("Ele1_Ele2_deta", &V_ae_Ele1_Ele2_deta_);
    RHTree->Branch("Ele1_Ele2_dphi", &V_ae_Ele1_Ele2_dphi_);

    genParticlesToken_ = consumes<std::vector<reco::GenParticle>>(
        iConfig.getParameter<edm::InputTag>("genParticles"));
}

GenAnalyzerAee::~GenAnalyzerAee()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void GenAnalyzerAee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	  using namespace edm;
    eventIDAee_ = iEvent.id().event();
    runIDAee_   = iEvent.id().run();
    lumiIDAee_  = iEvent.id().luminosityBlock();

    // clear vectors
    V_ae_genA_M_inv_.clear();
    V_ae_genA_M_.clear();
    V_ae_dR_Ele1_Ele2_.clear();
    V_ae_A_pt_.clear();   V_ae_A_eta_.clear();   V_ae_A_phi_.clear();
    V_ae_Ele1_pt_.clear(); V_ae_Ele1_eta_.clear(); V_ae_Ele1_phi_.clear();
    V_ae_Ele2_pt_.clear(); V_ae_Ele2_eta_.clear(); V_ae_Ele2_phi_.clear();
    V_ae_Ele1_Ele2_deta_.clear();
    V_ae_Ele1_Ele2_dphi_.clear();

    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    bool found = false;

    for (auto const& p : *genParticles) {

        if (std::abs(p.pdgId()) != 9000036) continue;
        if (p.numberOfDaughters() != 2) continue;

        // Optional: common status filters (adjust according to your sample)
        // if (p.status() != 62 && p.status() != 22) continue;

        const reco::Candidate* d1 = p.daughter(0);
        const reco::Candidate* d2 = p.daughter(1);

        if (std::abs(d1->pdgId()) != 11 || std::abs(d2->pdgId()) != 11) continue;

        // Optional: final state electrons (status 1)
        // if (d1->status() != 1 || d2->status() != 1) continue;

        found = true;

        auto ele1 = SetFourVector(d1->pt(), d1->eta(), d1->phi(), d1->mass());
        auto ele2 = SetFourVector(d2->pt(), d2->eta(), d2->phi(), d2->mass());

        auto vecA = ele1 + ele2;

        V_ae_genA_M_inv_.push_back( vecA.M() );
        V_ae_genA_M_.push_back(     p.mass() );

        float dr_e12 = reco::deltaR(*d1, *d2);
        V_ae_dR_Ele1_Ele2_.push_back(dr_e12);

        V_ae_A_pt_.push_back(  p.pt()  );
        V_ae_A_eta_.push_back( p.eta() );
        V_ae_A_phi_.push_back( p.phi() );

        // We arbitrarily label the first daughter as Ele1, second as Ele2
        // You can sort by pT if desired: if (d2->pt() > d1->pt()) swap(d1,d2);

        V_ae_Ele1_pt_.push_back( d1->pt()  );
        V_ae_Ele1_eta_.push_back(d1->eta() );
        V_ae_Ele1_phi_.push_back(d1->phi() );

        V_ae_Ele2_pt_.push_back( d2->pt()  );
        V_ae_Ele2_eta_.push_back(d2->eta() );
        V_ae_Ele2_phi_.push_back(d2->phi() );

        V_ae_Ele1_Ele2_deta_.push_back( std::abs(d1->eta() - d2->eta()) );
        V_ae_Ele1_Ele2_dphi_.push_back( std::abs(reco::deltaPhi(d1->phi(), d2->phi())) );

        // If you expect only one interesting A per event → break;
        // break;
    }

    nTotalAee++;
    if (found) {
        nPassedAee++;
        RHTree->Fill();
    }
}


void GenAnalyzerAee::beginJob() {
    nTotalAee = nPassedAee = 0;
}


void GenAnalyzerAee::endJob() {
    std::cout << "Events processed: " << nTotalAee << "\n"
              << "Events with A → ee (|pdgId|=25): " << nPassedAee << "  ("
              << std::fixed << std::setprecision(3)
              << (nTotalAee > 0 ? 100.0 * nPassedAee / nTotalAee : 0.0) << "%)\n";
}


void GenAnalyzerAee::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenAnalyzerAee);
