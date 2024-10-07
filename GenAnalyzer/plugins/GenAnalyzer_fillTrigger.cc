#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"
// vector<std::string> V_accept_trigger_name_;


TH1D *H_accept_trigger;
float V_accept_trigger;
vector<float> V_accept_trigger_;

//-----------------------now do what ever initialization is needed

void GenAnalyzer::branchesTrigger(TTree* tree, edm::Service<TFileService> &fs)
{
  H_accept_trigger     = fs->make<TH1D>("h_accept_trigger"   , "accept_trigger;accept_trigger;Events"                 ,  100,  0, 500);
  tree->Branch("Num_accepted_triggers",  &V_accept_trigger_);
  // RHTree->Branch("accepted_trigger_name",  &V_accept_trigger_name_);
}

// ---------------------- Fill tree with trigger info  ------------
void GenAnalyzer::fillTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  V_accept_trigger = -1111.1111;
  std::cout << " >>>>>> Checking TriggerResults" << std::endl;
  float pass_trigger = 0;
  // Study of trigger bit
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(triggerResultsToken_,hltresults);
  if (!hltresults.isValid()) {
    if(print_trigger) std::cout << "!!! Error in getting TriggerResults product from Event !!!" << std::endl;
  }
  int ntrigs = hltresults->size();
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);

  for (int itrig = 0; itrig != ntrigs; ++itrig)
  {

    std::string trigName = triggerNames.triggerName(itrig);
    bool accept = hltresults->accept(itrig);
    if (!(accept))
    {
      if(print_trigger) std::cout << " Rejected Triggers !!!!!!!!!!> " << trigName << std::endl;
      continue;
    }
    else
    {
      pass_trigger = pass_trigger+1;
      if(print_trigger) std::cout << " Accept Triggers---------> " << trigName << std::endl;

    }

    // char* cstr = new char[trigName.size() + 1];
    // std::strcpy(cstr, trigName.c_str());
    // V_accept_trigger_name_.push_back(cstr);

  }
  std::cout << " Number of Accept Triggers in this events ---------> " << pass_trigger << std::endl;



  V_accept_trigger    = pass_trigger;
  V_accept_trigger_.push_back(V_accept_trigger);
  H_accept_trigger->Fill(V_accept_trigger);
  RHTree->Fill();

  V_accept_trigger_.clear();
  // V_accept_trigger_name_.clear();


}
