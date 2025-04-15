// -*- C++ -*-
//
// Package:    Gen/GenAnalyzer
// Class:      GenAnalyzer
//
/**\class GenAnalyzer GenAnalyzer.cc Gen/GenAnalyzer/plugins/GenAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Bhim Bam
//         Created:  Tue, 03 Sep 2024 18:24:42 GMT
//
//


#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"

using reco::GenParticle;

int ntotal_event ;
int npassed_event ;

unsigned int runId_;
unsigned int lumiId_;
unsigned long long eventId_;


vector<float> V_att_genHiggs_M_inv_;
vector<float> V_att_genA1_M_inv_;
vector<float> V_att_genA2_M_inv_;
vector<float> V_att_genHiggs_M_;
vector<float> V_att_genA1_M_;
vector<float> V_att_genA2_M_;
vector<float> V_att_dR_A1_A2_;
vector<float> V_att_dR_H_A1_;
vector<float> V_att_dR_H_A2_;
vector<float> V_att_dR_A1_Tau1_;
vector<float> V_att_dR_A1_Tau2_;
vector<float> V_att_dR_A2_Tau3_;
vector<float> V_att_dR_A2_Tau4_;
vector<float> V_att_dR_Tau1_Tau2_;
vector<float> V_att_dR_Tau3_Tau4_;

vector<float> V_att_H_pt_;
vector<float> V_att_A1_pt_;
vector<float> V_att_A2_pt_;
vector<float> V_att_Tau1_pt_;
vector<float> V_att_Tau2_pt_;
vector<float> V_att_Tau3_pt_;
vector<float> V_att_Tau4_pt_;
vector<float> V_att_H_eta_;
vector<float> V_att_A1_eta_;
vector<float> V_att_A2_eta_;
vector<float> V_att_Tau1_eta_;
vector<float> V_att_Tau2_eta_;
vector<float> V_att_Tau3_eta_;
vector<float> V_att_Tau4_eta_;
vector<float> V_att_H_phi_;
vector<float> V_att_A1_phi_;
vector<float> V_att_A2_phi_;
vector<float> V_att_Tau1_phi_;
vector<float> V_att_Tau2_phi_;
vector<float> V_att_Tau3_phi_;
vector<float> V_att_Tau4_phi_;

vector<float> V_att_Tau1_Tau2_deta_;
vector<float> V_att_Tau1_Tau2_dphi_;
vector<float> V_att_Tau3_Tau4_deta_;
vector<float> V_att_Tau3_Tau4_dphi_;


TLorentzVector SetTaus(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
// :
// tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
  isDebug  = iConfig.getParameter<bool>("isDebug");
  print_trigger  = iConfig.getParameter<bool>("print_trigger");

  //now do what ever initialization is needed
  RHTree = fs->make<TTree>("RHTree","Gen info Tree");
  
  RHTree->Branch("Event",  &eventId_);
  RHTree->Branch("Run",  &runId_);
  RHTree->Branch("LumiSection",  &lumiId_);

  RHTree->Branch("GenHiggs_inv",  &V_att_genHiggs_M_inv_);
  RHTree->Branch("GenA1_inv",  &V_att_genA1_M_inv_);
  RHTree->Branch("GenA2_inv",  &V_att_genA2_M_inv_);
  RHTree->Branch("GenHiggs",  &V_att_genHiggs_M_);
  RHTree->Branch("GenA1",  &V_att_genA1_M_);
  RHTree->Branch("GenA2",  &V_att_genA2_M_);
  RHTree->Branch("dR_A1_A2",  &V_att_dR_A1_A2_);
  RHTree->Branch("dR_H_A1",  &V_att_dR_H_A1_);
  RHTree->Branch("dR_H_A2",  &V_att_dR_H_A2_);
  RHTree->Branch("dR_A1_Tau1",  &V_att_dR_A1_Tau1_);
  RHTree->Branch("dR_A1_Tau2",  &V_att_dR_A1_Tau2_);
  RHTree->Branch("dR_A2_Tau3",  &V_att_dR_A2_Tau3_);
  RHTree->Branch("dR_A2_Tau4",  &V_att_dR_A2_Tau4_);
  RHTree->Branch("dR_Tau1_Tau2",  &V_att_dR_Tau1_Tau2_);
  RHTree->Branch("dR_Tau3_Tau4",  &V_att_dR_Tau3_Tau4_);

  RHTree->Branch("H_pt",  &V_att_H_pt_);
  RHTree->Branch("A1_pt",  &V_att_A1_pt_);
  RHTree->Branch("A2_pt",  &V_att_A2_pt_);
  RHTree->Branch("Tau1_pt",  &V_att_Tau1_pt_);
  RHTree->Branch("Tau2_pt",  &V_att_Tau2_pt_);
  RHTree->Branch("Tau3_pt",  &V_att_Tau3_pt_);
  RHTree->Branch("Tau4_pt",  &V_att_Tau4_pt_);
  RHTree->Branch("H_eta",  &V_att_H_eta_);
  RHTree->Branch("A1_eta",  &V_att_A1_eta_);
  RHTree->Branch("A2_eta",  &V_att_A2_eta_);
  RHTree->Branch("Tau1_eta",  &V_att_Tau1_eta_);
  RHTree->Branch("Tau2_eta",  &V_att_Tau2_eta_);
  RHTree->Branch("Tau3_eta",  &V_att_Tau3_eta_);
  RHTree->Branch("Tau4_eta",  &V_att_Tau4_eta_);
  RHTree->Branch("H_phi",  &V_att_H_phi_);
  RHTree->Branch("A1_phi",  &V_att_A1_phi_);
  RHTree->Branch("A2_phi",  &V_att_A2_phi_);
  RHTree->Branch("Tau1_phi",  &V_att_Tau1_phi_);
  RHTree->Branch("Tau2_phi",  &V_att_Tau2_phi_);
  RHTree->Branch("Tau3_phi",  &V_att_Tau3_phi_);
  RHTree->Branch("Tau4_phi",  &V_att_Tau4_phi_);

  RHTree->Branch("Tau1_Tau2_deta",  &V_att_Tau1_Tau2_deta_);
  RHTree->Branch("Tau1_Tau2_dphi",  &V_att_Tau1_Tau2_dphi_);
  RHTree->Branch("Tau3_Tau4_deta",  &V_att_Tau3_Tau4_deta_);
  RHTree->Branch("Tau3_Tau4_dphi",  &V_att_Tau3_Tau4_dphi_);
   
  genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));

}


GenAnalyzer::~GenAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();

  V_att_genHiggs_M_inv_.clear();
  V_att_genA1_M_inv_.clear();
  V_att_genA2_M_inv_.clear();
  V_att_genHiggs_M_.clear();
  V_att_genA1_M_.clear();
  V_att_genA2_M_.clear();

  V_att_dR_A1_A2_.clear();
  V_att_dR_H_A1_.clear();
  V_att_dR_H_A2_.clear();
  V_att_dR_A1_Tau1_.clear();
  V_att_dR_A1_Tau2_.clear();
  V_att_dR_A2_Tau3_.clear();
  V_att_dR_A2_Tau4_.clear();
  V_att_dR_Tau1_Tau2_.clear();
  V_att_dR_Tau3_Tau4_.clear();

  V_att_H_pt_.clear();
  V_att_A1_pt_.clear();
  V_att_A2_pt_.clear();
  V_att_Tau1_pt_.clear();
  V_att_Tau2_pt_.clear();
  V_att_Tau3_pt_.clear();
  V_att_Tau4_pt_.clear();
  V_att_H_eta_.clear();
  V_att_A1_eta_.clear();
  V_att_A2_eta_.clear();
  V_att_Tau1_eta_.clear();
  V_att_Tau2_eta_.clear();
  V_att_Tau3_eta_.clear();
  V_att_Tau4_eta_.clear();
  V_att_H_phi_.clear();
  V_att_A1_phi_.clear();
  V_att_A2_phi_.clear();
  V_att_Tau1_phi_.clear();
  V_att_Tau2_phi_.clear();
  V_att_Tau3_phi_.clear();
  V_att_Tau4_phi_.clear();

  V_att_Tau1_Tau2_deta_.clear();
  V_att_Tau1_Tau2_dphi_.clear();
  V_att_Tau3_Tau4_deta_.clear();
  V_att_Tau3_Tau4_dphi_.clear();



  float genHiggs_mass_inv = -1111.1111;
  float genA1_mass_inv = -1111.1111;
  float genA2_mass_inv = -1111.1111;
  float genHiggs_mass = -1111.1111;
  float genA1_mass = -1111.1111;
  float genA2_mass = -1111.1111;
  float A1_A2_dR = -1111.1111;
  float H_A1_dR = -1111.1111;
  float H_A2_dR = -1111.1111;
  float A1_Tau1_dR = -1111.1111;
  float A1_Tau2_dR = -1111.1111;
  float A2_Tau3_dR = -1111.1111;
  float A2_Tau4_dR = -1111.1111;
  float Tau1_Tau2_dR = -1111.1111;
  float Tau3_Tau4_dR = -1111.1111;

  float H_pt = -1111.1111;
  float A1_pt = -1111.1111;
  float A2_pt = -1111.1111;
  float Tau1_pt = -1111.1111;
  float Tau2_pt = -1111.1111;
  float Tau3_pt = -1111.1111;
  float Tau4_pt = -1111.1111;
  float H_eta = -1111.1111;
  float A1_eta = -1111.1111;
  float A2_eta = -1111.1111;
  float Tau1_eta = -1111.1111;
  float Tau2_eta = -1111.1111;
  float Tau3_eta = -1111.1111;
  float Tau4_eta = -1111.1111;
  float H_phi = -1111.1111;
  float A1_phi = -1111.1111;
  float A2_phi = -1111.1111;
  float Tau1_phi = -1111.1111;
  float Tau2_phi = -1111.1111;
  float Tau3_phi = -1111.1111;
  float Tau4_phi = -1111.1111;

  float Tau1_Tau2_deta = -1111.1111;
  float Tau1_Tau2_dphi = -1111.1111;
  float Tau3_Tau4_deta = -1111.1111;
  float Tau3_Tau4_dphi = -1111.1111;

  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,   genParticles);


  // unsigned int NAs = 0;
  // unsigned int NTau_fromA = 0;
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( abs(iGen->pdgId()) != 15 || abs(iGen->mother()->pdgId()) != 25) continue;
  //   NTau_fromA++;
  // }
  // std::cout << "  >>>>>> Number Tau from  A <<<<<"<<"    "  <<  NTau_fromA << std::endl;
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( abs(iGen->pdgId()) != 25 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15) continue;
  //   NAs++;
  // }
  // std::cout << "  >>>>>> Number of A giving Tau <<<<<"<<"    "  <<  NAs << std::endl;


  bool pass = false;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {

    if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->status()) != 2 || abs(iGen->daughter(0)->daughter(1)->status()) != 2 || abs(iGen->daughter(1)->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->daughter(1)->status()) != 2 ) continue;
    pass = true;

    TLorentzVector GenTau1  = SetTaus(iGen->daughter(0)->daughter(0)->pt(), iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(0)->mass());
    TLorentzVector GenTau2  = SetTaus(iGen->daughter(0)->daughter(1)->pt(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->daughter(1)->mass());
    TLorentzVector GenA1 = GenTau1 + GenTau2;
    TLorentzVector GenTau3  = SetTaus(iGen->daughter(1)->daughter(0)->pt(), iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(0)->mass());
    TLorentzVector GenTau4  = SetTaus(iGen->daughter(1)->daughter(1)->pt(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->daughter(1)->mass());
    TLorentzVector GenA2 = GenTau3 + GenTau4;
    TLorentzVector GenHiggs = GenA1 + GenA2;

    genHiggs_mass_inv = GenHiggs.M();
    genA1_mass_inv = GenA1.M();
    genA2_mass_inv = GenA2.M();
    genHiggs_mass = iGen->mass();
    genA1_mass = iGen->daughter(0)->mass();
    genA2_mass = iGen->daughter(1)->mass();

    V_att_genHiggs_M_inv_.push_back( genHiggs_mass_inv );
    V_att_genA1_M_inv_.push_back( genA1_mass_inv );
    V_att_genA2_M_inv_.push_back( genA2_mass_inv );
    V_att_genHiggs_M_.push_back( genHiggs_mass );
    V_att_genA1_M_.push_back( genA1_mass );
    V_att_genA2_M_.push_back( genA2_mass );


    float dR_A1_A2 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_H_A1 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_H_A2 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_A1_Tau1 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_A1_Tau2 = reco::deltaR( iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_A2_Tau3 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_A2_Tau4 = reco::deltaR( iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_Tau1_Tau2 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi());
    float dR_Tau3_Tau4 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi());

    A1_A2_dR = dR_A1_A2;
    H_A1_dR = dR_H_A1;
    H_A2_dR = dR_H_A2;
    A1_Tau1_dR = dR_A1_Tau1;
    A1_Tau2_dR = dR_A1_Tau2;
    A2_Tau3_dR = dR_A2_Tau3;
    A2_Tau4_dR = dR_A2_Tau4;
    Tau1_Tau2_dR = dR_Tau1_Tau2;
    Tau3_Tau4_dR = dR_Tau3_Tau4;

    V_att_dR_A1_A2_.push_back( A1_A2_dR );
    V_att_dR_H_A1_.push_back( H_A1_dR );
    V_att_dR_H_A2_.push_back( H_A2_dR );
    V_att_dR_A1_Tau1_.push_back( A1_Tau1_dR );
    V_att_dR_A1_Tau2_.push_back( A1_Tau2_dR );
    V_att_dR_A2_Tau3_.push_back( A2_Tau3_dR );
    V_att_dR_A2_Tau4_.push_back( A2_Tau4_dR );
    V_att_dR_Tau1_Tau2_.push_back( Tau1_Tau2_dR );
    V_att_dR_Tau3_Tau4_.push_back( Tau3_Tau4_dR );

    H_pt = iGen->pt();
    A1_pt = iGen->daughter(0)->pt();
    A2_pt = iGen->daughter(1)->pt();
    Tau1_pt = iGen->daughter(0)->daughter(0)->pt();
    Tau2_pt = iGen->daughter(0)->daughter(1)->pt();
    Tau3_pt = iGen->daughter(1)->daughter(0)->pt();
    Tau4_pt = iGen->daughter(1)->daughter(1)->pt();
    H_eta = iGen->eta();
    A1_eta = iGen->daughter(0)->eta();
    A2_eta = iGen->daughter(1)->eta();
    Tau1_eta = iGen->daughter(0)->daughter(0)->eta();
    Tau2_eta = iGen->daughter(0)->daughter(1)->eta();
    Tau3_eta = iGen->daughter(1)->daughter(0)->eta();
    Tau4_eta = iGen->daughter(1)->daughter(1)->eta();
    H_phi = iGen->phi();
    A1_phi = iGen->daughter(0)->phi();
    A2_phi = iGen->daughter(1)->phi();
    Tau1_phi = iGen->daughter(0)->daughter(0)->phi();
    Tau2_phi = iGen->daughter(0)->daughter(1)->phi();
    Tau3_phi = iGen->daughter(1)->daughter(0)->phi();
    Tau4_phi = iGen->daughter(1)->daughter(1)->phi();

    V_att_H_pt_.push_back(H_pt );
    V_att_A1_pt_.push_back(A1_pt );
    V_att_A2_pt_.push_back(A2_pt );
    V_att_Tau1_pt_.push_back(Tau1_pt );
    V_att_Tau2_pt_.push_back(Tau2_pt );
    V_att_Tau3_pt_.push_back(Tau3_pt );
    V_att_Tau4_pt_.push_back(Tau4_pt );
    V_att_H_eta_.push_back(H_eta );
    V_att_A1_eta_.push_back(A1_eta );
    V_att_A2_eta_.push_back(A2_eta );
    V_att_Tau1_eta_.push_back(Tau1_eta );
    V_att_Tau2_eta_.push_back(Tau2_eta );
    V_att_Tau3_eta_.push_back(Tau3_eta );
    V_att_Tau4_eta_.push_back(Tau4_eta );
    V_att_H_phi_.push_back(H_phi );
    V_att_A1_phi_.push_back(A1_phi );
    V_att_A2_phi_.push_back(A2_phi );
    V_att_Tau1_phi_.push_back(Tau1_phi );
    V_att_Tau2_phi_.push_back(Tau2_phi );
    V_att_Tau3_phi_.push_back(Tau3_phi );
    V_att_Tau4_phi_.push_back(Tau4_phi );

    Tau1_Tau2_deta = abs(Tau1_eta-Tau2_eta);
    Tau3_Tau4_deta = abs(Tau3_eta-Tau4_eta);
  
    Tau1_Tau2_dphi = abs(reco::deltaPhi( Tau1_phi, Tau2_phi )); 
    Tau3_Tau4_dphi = abs(reco::deltaPhi( Tau3_phi, Tau4_phi )); 

    V_att_Tau1_Tau2_deta_.push_back( Tau1_Tau2_deta );
    V_att_Tau3_Tau4_deta_.push_back( Tau3_Tau4_deta );
    V_att_Tau1_Tau2_dphi_.push_back( Tau1_Tau2_dphi );
    V_att_Tau3_Tau4_dphi_.push_back( Tau3_Tau4_dphi );


    // std::cout << "  >>>>>> Higgs gen (35) <<<<<"<<"<<< status: "<<iGen->status()<<"<<<pt:  "<<iGen->pt()<<"<<<eta:  "<<iGen->eta()<<"<<<phi:  "<<iGen->phi()<<"<<<mass:  "<<iGen->mass() << std::endl;
    // std::cout << "  >>>>>> Higgs LorentzVector (35) <<<<<"<<"<<<mass:  "<<GenHiggs.M() << std::endl;
    // std::cout << "  >>>>>> A1 gen (25) <<<<<"<<"<<< status: "<<iGen->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->mass() << std::endl;
    // std::cout << "  >>>>>> A1 LorentzVector (25) <<<<<"<<"<<<mass:  "<< GenA1.M() << std::endl;
    // std::cout << "  >>>>>> A2 gen (25) <<<<<"<<"<<< status: "<<iGen->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->mass() << std::endl;
    // std::cout << "  >>>>>> A2 LorentzVector (25) <<<<<"<<"<<<mass:  "<< GenA2.M() << std::endl;
    // std::cout << "  >>>>>> Tau1 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(0)->mass() << std::endl;
    // std::cout << "  >>>>>> Tau2 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(1)->mass() << std::endl;
    // std::cout << "  >>>>>> Tau3 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(0)->mass() << std::endl;
    // std::cout << "  >>>>>> Tau4 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(1)->mass() << std::endl;
    //
    // std::cout << "  >>>>>> dR_A1_A2:  "<<dR_A1_A2<< std::endl;
    // std::cout << "  >>>>>> dR_H_A1:  "<<dR_H_A1<< std::endl;
    // std::cout << "  >>>>>> dR_H_A2:  "<<dR_H_A2<< std::endl;
    // std::cout << "  >>>>>> dR_A1_Tau1:  "<<dR_A1_Tau1<< std::endl;
    // std::cout << "  >>>>>> dR_A1_Tau2:  "<<dR_A1_Tau2<< std::endl;
    // std::cout << "  >>>>>> dR_A2_Tau3:  "<<dR_A2_Tau3<< std::endl;
    // std::cout << "  >>>>>> dR_A2_Tau4:  "<<dR_A2_Tau4<< std::endl;
    // std::cout << "  >>>>>> dR_Tau1_Tau2:  "<<dR_Tau1_Tau2<< std::endl;
    // std::cout << "  >>>>>> dR_Tau3_Tau4:  "<<dR_Tau3_Tau4<< std::endl;

  } // gen particle collection, looping over Higgs. 
  ntotal_event++;
  if (pass) {
    npassed_event++;
    //fillTrigger( iEvent, iSetup );
    RHTree->Fill();

  }



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
GenAnalyzer::beginJob()
{
  ntotal_event = 0;
  npassed_event = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer::endJob()
{
  std::cout << "  >>>>>> Total events selected events <<<<<  "<<npassed_event<<"/"<<ntotal_event<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
