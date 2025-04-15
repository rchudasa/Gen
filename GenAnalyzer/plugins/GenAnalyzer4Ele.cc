// -*- C++ -*-
//
// Package:    Gen/GenAnalyzer4Ele
// Class:      GenAnalyzer4Ele
//
/**\class GenAnalyzer4Ele GenAnalyzer4Ele.cc Gen/GenAnalyzer/plugins/GenAnalyzer4Ele.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ruchi Chudasama
//         Created:  Wed, 04 Dec 2024 15:38:50 GMT
//
//


#include "Gen/GenAnalyzer/interface/GenAnalyzer4Ele.h"

using reco::GenParticle;

int nTotal ;
int nPassed ;

unsigned int runID_;
unsigned int lumiID_;
unsigned long long eventID_;


vector<float> V_aee_genHiggs_M_inv_;
vector<float> V_aee_genA1_M_inv_;
vector<float> V_aee_genA2_M_inv_;
vector<float> V_aee_genHiggs_M_;
vector<float> V_aee_genA1_M_;
vector<float> V_aee_genA2_M_;
vector<float> V_aee_dR_A1_A2_;
vector<float> V_aee_dR_H_A1_;
vector<float> V_aee_dR_H_A2_;
vector<float> V_aee_dR_A1_Ele1_;
vector<float> V_aee_dR_A1_Ele2_;
vector<float> V_aee_dR_A2_Ele3_;
vector<float> V_aee_dR_A2_Ele4_;
vector<float> V_aee_dR_Ele1_Ele2_;
vector<float> V_aee_dR_Ele3_Ele4_;

vector<float> V_aee_H_pt_;
vector<float> V_aee_A1_pt_;
vector<float> V_aee_A2_pt_;
vector<float> V_aee_Ele1_pt_;
vector<float> V_aee_Ele2_pt_;
vector<float> V_aee_Ele3_pt_;
vector<float> V_aee_Ele4_pt_;
vector<float> V_aee_H_eta_;
vector<float> V_aee_A1_eta_;
vector<float> V_aee_A2_eta_;
vector<float> V_aee_Ele1_eta_;
vector<float> V_aee_Ele2_eta_;
vector<float> V_aee_Ele3_eta_;
vector<float> V_aee_Ele4_eta_;
vector<float> V_aee_H_phi_;
vector<float> V_aee_A1_phi_;
vector<float> V_aee_A2_phi_;
vector<float> V_aee_Ele1_phi_;
vector<float> V_aee_Ele2_phi_;
vector<float> V_aee_Ele3_phi_;
vector<float> V_aee_Ele4_phi_;

vector<float> V_aee_Ele1_Ele2_deta_;
vector<float> V_aee_Ele1_Ele2_dphi_;
vector<float> V_aee_Ele3_Ele4_deta_;
vector<float> V_aee_Ele3_Ele4_dphi_;


TLorentzVector SetEles(Float_t Ele_pt, Float_t Ele_eta, Float_t Ele_phi, Float_t Ele_mass){
  TLorentzVector Ele_Candidate;
  Ele_Candidate.SetPtEtaPhiM(Ele_pt, Ele_eta, Ele_phi, Ele_mass);
  return Ele_Candidate;
}

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer4Ele::GenAnalyzer4Ele(const edm::ParameterSet& iConfig)
// :
// tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
  isDebug  = iConfig.getParameter<bool>("isDebug");
  print_trigger  = iConfig.getParameter<bool>("print_trigger");

  //now do what ever initialization is needed
  RHTree = fs->make<TTree>("RHTree","Gen info Tree");
  
  RHTree->Branch("Event",  &eventID_);
  RHTree->Branch("Run",  &runID_);
  RHTree->Branch("LumiSection",  &lumiID_);

  RHTree->Branch("GenHiggs_inv",  &V_aee_genHiggs_M_inv_);
  RHTree->Branch("GenA1_inv",  &V_aee_genA1_M_inv_);
  RHTree->Branch("GenA2_inv",  &V_aee_genA2_M_inv_);
  RHTree->Branch("GenHiggs",  &V_aee_genHiggs_M_);
  RHTree->Branch("GenA1",  &V_aee_genA1_M_);
  RHTree->Branch("GenA2",  &V_aee_genA2_M_);
  RHTree->Branch("dR_A1_A2",  &V_aee_dR_A1_A2_);
  RHTree->Branch("dR_H_A1",  &V_aee_dR_H_A1_);
  RHTree->Branch("dR_H_A2",  &V_aee_dR_H_A2_);
  RHTree->Branch("dR_A1_Ele1",  &V_aee_dR_A1_Ele1_);
  RHTree->Branch("dR_A1_Ele2",  &V_aee_dR_A1_Ele2_);
  RHTree->Branch("dR_A2_Ele3",  &V_aee_dR_A2_Ele3_);
  RHTree->Branch("dR_A2_Ele4",  &V_aee_dR_A2_Ele4_);
  RHTree->Branch("dR_Ele1_Ele2",  &V_aee_dR_Ele1_Ele2_);
  RHTree->Branch("dR_Ele3_Ele4",  &V_aee_dR_Ele3_Ele4_);

  RHTree->Branch("H_pt",  &V_aee_H_pt_);
  RHTree->Branch("A1_pt",  &V_aee_A1_pt_);
  RHTree->Branch("A2_pt",  &V_aee_A2_pt_);
  RHTree->Branch("Ele1_pt",  &V_aee_Ele1_pt_);
  RHTree->Branch("Ele2_pt",  &V_aee_Ele2_pt_);
  RHTree->Branch("Ele3_pt",  &V_aee_Ele3_pt_);
  RHTree->Branch("Ele4_pt",  &V_aee_Ele4_pt_);
  RHTree->Branch("H_eta",  &V_aee_H_eta_);
  RHTree->Branch("A1_eta",  &V_aee_A1_eta_);
  RHTree->Branch("A2_eta",  &V_aee_A2_eta_);
  RHTree->Branch("Ele1_eta",  &V_aee_Ele1_eta_);
  RHTree->Branch("Ele2_eta",  &V_aee_Ele2_eta_);
  RHTree->Branch("Ele3_eta",  &V_aee_Ele3_eta_);
  RHTree->Branch("Ele4_eta",  &V_aee_Ele4_eta_);
  RHTree->Branch("H_phi",  &V_aee_H_phi_);
  RHTree->Branch("A1_phi",  &V_aee_A1_phi_);
  RHTree->Branch("A2_phi",  &V_aee_A2_phi_);
  RHTree->Branch("Ele1_phi",  &V_aee_Ele1_phi_);
  RHTree->Branch("Ele2_phi",  &V_aee_Ele2_phi_);
  RHTree->Branch("Ele3_phi",  &V_aee_Ele3_phi_);
  RHTree->Branch("Ele4_phi",  &V_aee_Ele4_phi_);

  RHTree->Branch("Ele1_Ele2_deta",  &V_aee_Ele1_Ele2_deta_);
  RHTree->Branch("Ele1_Ele2_dphi",  &V_aee_Ele1_Ele2_dphi_);
  RHTree->Branch("Ele3_Ele4_deta",  &V_aee_Ele3_Ele4_deta_);
  RHTree->Branch("Ele3_Ele4_dphi",  &V_aee_Ele3_Ele4_dphi_);
   
  genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));

}


GenAnalyzer4Ele::~GenAnalyzer4Ele()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer4Ele::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  eventID_ = iEvent.id().event();
  runID_ = iEvent.id().run();
  lumiID_ = iEvent.id().luminosityBlock();

  V_aee_genHiggs_M_inv_.clear();
  V_aee_genA1_M_inv_.clear();
  V_aee_genA2_M_inv_.clear();
  V_aee_genHiggs_M_.clear();
  V_aee_genA1_M_.clear();
  V_aee_genA2_M_.clear();

  V_aee_dR_A1_A2_.clear();
  V_aee_dR_H_A1_.clear();
  V_aee_dR_H_A2_.clear();
  V_aee_dR_A1_Ele1_.clear();
  V_aee_dR_A1_Ele2_.clear();
  V_aee_dR_A2_Ele3_.clear();
  V_aee_dR_A2_Ele4_.clear();
  V_aee_dR_Ele1_Ele2_.clear();
  V_aee_dR_Ele3_Ele4_.clear();

  V_aee_H_pt_.clear();
  V_aee_A1_pt_.clear();
  V_aee_A2_pt_.clear();
  V_aee_Ele1_pt_.clear();
  V_aee_Ele2_pt_.clear();
  V_aee_Ele3_pt_.clear();
  V_aee_Ele4_pt_.clear();
  V_aee_H_eta_.clear();
  V_aee_A1_eta_.clear();
  V_aee_A2_eta_.clear();
  V_aee_Ele1_eta_.clear();
  V_aee_Ele2_eta_.clear();
  V_aee_Ele3_eta_.clear();
  V_aee_Ele4_eta_.clear();
  V_aee_H_phi_.clear();
  V_aee_A1_phi_.clear();
  V_aee_A2_phi_.clear();
  V_aee_Ele1_phi_.clear();
  V_aee_Ele2_phi_.clear();
  V_aee_Ele3_phi_.clear();
  V_aee_Ele4_phi_.clear();

  V_aee_Ele1_Ele2_deta_.clear();
  V_aee_Ele1_Ele2_dphi_.clear();
  V_aee_Ele3_Ele4_deta_.clear();
  V_aee_Ele3_Ele4_dphi_.clear();



  float genHiggs_mass_inv = -1111.1111;
  float genA1_mass_inv = -1111.1111;
  float genA2_mass_inv = -1111.1111;
  float genHiggs_mass = -1111.1111;
  float genA1_mass = -1111.1111;
  float genA2_mass = -1111.1111;
  float A1_A2_dR = -1111.1111;
  float H_A1_dR = -1111.1111;
  float H_A2_dR = -1111.1111;
  float A1_Ele1_dR = -1111.1111;
  float A1_Ele2_dR = -1111.1111;
  float A2_Ele3_dR = -1111.1111;
  float A2_Ele4_dR = -1111.1111;
  float Ele1_Ele2_dR = -1111.1111;
  float Ele3_Ele4_dR = -1111.1111;

  float H_pt = -1111.1111;
  float A1_pt = -1111.1111;
  float A2_pt = -1111.1111;
  float Ele1_pt = -1111.1111;
  float Ele2_pt = -1111.1111;
  float Ele3_pt = -1111.1111;
  float Ele4_pt = -1111.1111;
  float H_eta = -1111.1111;
  float A1_eta = -1111.1111;
  float A2_eta = -1111.1111;
  float Ele1_eta = -1111.1111;
  float Ele2_eta = -1111.1111;
  float Ele3_eta = -1111.1111;
  float Ele4_eta = -1111.1111;
  float H_phi = -1111.1111;
  float A1_phi = -1111.1111;
  float A2_phi = -1111.1111;
  float Ele1_phi = -1111.1111;
  float Ele2_phi = -1111.1111;
  float Ele3_phi = -1111.1111;
  float Ele4_phi = -1111.1111;

  float Ele1_Ele2_deta = -1111.1111;
  float Ele1_Ele2_dphi = -1111.1111;
  float Ele3_Ele4_deta = -1111.1111;
  float Ele3_Ele4_dphi = -1111.1111;

  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,   genParticles);



  bool pass = false;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {

    if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 ) continue;
    //if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 11 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 11 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 11 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 11 ) continue;
    //if ( abs(iGen->daughter(0)->daughter(0)->status()) != 1 || abs(iGen->daughter(0)->daughter(1)->status()) != 1 || abs(iGen->daughter(1)->daughter(0)->status()) != 1 || abs(iGen->daughter(1)->daughter(1)->status()) != 1 ) continue;

    std::cout << "Daughter pdgID:"<< iGen->daughter(0)->pdgId() << " status:"<< iGen->daughter(0)->status() ;
    std::cout << " grand Daughter pdgID:"<< iGen->daughter(0)->daughter(0)->pdgId() << " status:"<< iGen->daughter(0)->daughter(0)->status() ;
    std::cout << " grand Daughter 2 pdgID:"<< iGen->daughter(0)->daughter(1)->pdgId() << " status:"<< iGen->daughter(0)->daughter(1)->status()<<std::endl;
    std::cout << " Daughter 2 pdgID:"<< iGen->daughter(1)->pdgId() << " status:"<< iGen->daughter(1)->status();
    std::cout << " grand Daughter pdgID:"<< iGen->daughter(1)->daughter(0)->pdgId() << " status:"<< iGen->daughter(1)->daughter(0)->status() ;
    std::cout << " grand Daughter 2 pdgID:"<< iGen->daughter(1)->daughter(1)->pdgId() << " status:"<< iGen->daughter(1)->daughter(1)->status()<<std::endl<<std::endl;

    //if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 36 || iGen->daughter(1)->pdgId() != 36 ) continue;
    //if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 11 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 11 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 11 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 11 ) continue;
    //if ( abs(iGen->daughter(0)->daughter(0)->status()) != 2 || abs(iGen->daughter(0)->daughter(1)->status()) != 2 || abs(iGen->daughter(1)->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->daughter(1)->status()) != 2 ) continue;
    pass = true;

    TLorentzVector GenEle1  = SetEles(iGen->daughter(0)->daughter(0)->pt(), iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(0)->mass());
    TLorentzVector GenEle2  = SetEles(iGen->daughter(0)->daughter(1)->pt(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->daughter(1)->mass());
    TLorentzVector GenA1 = GenEle1 + GenEle2;
    TLorentzVector GenEle3  = SetEles(iGen->daughter(1)->daughter(0)->pt(), iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(0)->mass());
    TLorentzVector GenEle4  = SetEles(iGen->daughter(1)->daughter(1)->pt(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->daughter(1)->mass());
    TLorentzVector GenA2 = GenEle3 + GenEle4;
    TLorentzVector GenHiggs = GenA1 + GenA2;

    genHiggs_mass_inv = GenHiggs.M();
    genA1_mass_inv = GenA1.M();
    genA2_mass_inv = GenA2.M();
    genHiggs_mass = iGen->mass();
    genA1_mass = iGen->daughter(0)->mass();
    genA2_mass = iGen->daughter(1)->mass();

    V_aee_genHiggs_M_inv_.push_back( genHiggs_mass_inv );
    V_aee_genA1_M_inv_.push_back( genA1_mass_inv );
    V_aee_genA2_M_inv_.push_back( genA2_mass_inv );
    V_aee_genHiggs_M_.push_back( genHiggs_mass );
    V_aee_genA1_M_.push_back( genA1_mass );
    V_aee_genA2_M_.push_back( genA2_mass );


    float dR_A1_A2 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_H_A1 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_H_A2 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_A1_Ele1 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_A1_Ele2 = reco::deltaR( iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
    float dR_A2_Ele3 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_A2_Ele4 = reco::deltaR( iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
    float dR_Ele1_Ele2 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi());
    float dR_Ele3_Ele4 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi());

    A1_A2_dR = dR_A1_A2;
    H_A1_dR = dR_H_A1;
    H_A2_dR = dR_H_A2;
    A1_Ele1_dR = dR_A1_Ele1;
    A1_Ele2_dR = dR_A1_Ele2;
    A2_Ele3_dR = dR_A2_Ele3;
    A2_Ele4_dR = dR_A2_Ele4;
    Ele1_Ele2_dR = dR_Ele1_Ele2;
    Ele3_Ele4_dR = dR_Ele3_Ele4;

    V_aee_dR_A1_A2_.push_back( A1_A2_dR );
    V_aee_dR_H_A1_.push_back( H_A1_dR );
    V_aee_dR_H_A2_.push_back( H_A2_dR );
    V_aee_dR_A1_Ele1_.push_back( A1_Ele1_dR );
    V_aee_dR_A1_Ele2_.push_back( A1_Ele2_dR );
    V_aee_dR_A2_Ele3_.push_back( A2_Ele3_dR );
    V_aee_dR_A2_Ele4_.push_back( A2_Ele4_dR );
    V_aee_dR_Ele1_Ele2_.push_back( Ele1_Ele2_dR );
    V_aee_dR_Ele3_Ele4_.push_back( Ele3_Ele4_dR );

    H_pt = iGen->pt();
    A1_pt = iGen->daughter(0)->pt();
    A2_pt = iGen->daughter(1)->pt();
    Ele1_pt = iGen->daughter(0)->daughter(0)->pt();
    Ele2_pt = iGen->daughter(0)->daughter(1)->pt();
    Ele3_pt = iGen->daughter(1)->daughter(0)->pt();
    Ele4_pt = iGen->daughter(1)->daughter(1)->pt();
    H_eta = iGen->eta();
    A1_eta = iGen->daughter(0)->eta();
    A2_eta = iGen->daughter(1)->eta();
    Ele1_eta = iGen->daughter(0)->daughter(0)->eta();
    Ele2_eta = iGen->daughter(0)->daughter(1)->eta();
    Ele3_eta = iGen->daughter(1)->daughter(0)->eta();
    Ele4_eta = iGen->daughter(1)->daughter(1)->eta();
    H_phi = iGen->phi();
    A1_phi = iGen->daughter(0)->phi();
    A2_phi = iGen->daughter(1)->phi();
    Ele1_phi = iGen->daughter(0)->daughter(0)->phi();
    Ele2_phi = iGen->daughter(0)->daughter(1)->phi();
    Ele3_phi = iGen->daughter(1)->daughter(0)->phi();
    Ele4_phi = iGen->daughter(1)->daughter(1)->phi();

    V_aee_H_pt_.push_back(H_pt );
    V_aee_A1_pt_.push_back(A1_pt );
    V_aee_A2_pt_.push_back(A2_pt );
    V_aee_Ele1_pt_.push_back(Ele1_pt );
    V_aee_Ele2_pt_.push_back(Ele2_pt );
    V_aee_Ele3_pt_.push_back(Ele3_pt );
    V_aee_Ele4_pt_.push_back(Ele4_pt );
    V_aee_H_eta_.push_back(H_eta );
    V_aee_A1_eta_.push_back(A1_eta );
    V_aee_A2_eta_.push_back(A2_eta );
    V_aee_Ele1_eta_.push_back(Ele1_eta );
    V_aee_Ele2_eta_.push_back(Ele2_eta );
    V_aee_Ele3_eta_.push_back(Ele3_eta );
    V_aee_Ele4_eta_.push_back(Ele4_eta );
    V_aee_H_phi_.push_back(H_phi );
    V_aee_A1_phi_.push_back(A1_phi );
    V_aee_A2_phi_.push_back(A2_phi );
    V_aee_Ele1_phi_.push_back(Ele1_phi );
    V_aee_Ele2_phi_.push_back(Ele2_phi );
    V_aee_Ele3_phi_.push_back(Ele3_phi );
    V_aee_Ele4_phi_.push_back(Ele4_phi );

    Ele1_Ele2_deta = abs(Ele1_eta-Ele2_eta);
    Ele3_Ele4_deta = abs(Ele3_eta-Ele4_eta);
  
    Ele1_Ele2_dphi = abs(reco::deltaPhi( Ele1_phi, Ele2_phi )); 
    Ele3_Ele4_dphi = abs(reco::deltaPhi( Ele3_phi, Ele4_phi )); 

    V_aee_Ele1_Ele2_deta_.push_back( Ele1_Ele2_deta );
    V_aee_Ele3_Ele4_deta_.push_back( Ele3_Ele4_deta );
    V_aee_Ele1_Ele2_dphi_.push_back( Ele1_Ele2_dphi );
    V_aee_Ele3_Ele4_dphi_.push_back( Ele3_Ele4_dphi );


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
  nTotal++;
  if (pass) {
    nPassed++;
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
GenAnalyzer4Ele::beginJob()
{
  nTotal = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer4Ele::endJob()
{
  std::cout << "  >>>>>> Total events selected events <<<<<  "<<nPassed<<"/"<<nTotal<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer4Ele::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(GenAnalyzer4Ele);
