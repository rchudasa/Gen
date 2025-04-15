import os
cfg='GenAnalyzer/python/conFig4Ele_cfg.py'
#cfg='GenAnalyzer/python/conFig_cfg.py'
# inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGeneration/gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/crab_gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/230510_053835/0000/GEN_HToAAToTauTau_M10_2018UL_90.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/GEN_HToAAToTauTau_M13_2018UL.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/analysis/MCGeneration/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/Unboosted_GEN_HToAAToTauTau_M8_2018UL.root'
# inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGeneration/signal_withTrigger/GEN_PreMix_HToAATo4Tau_M_3p7_pythia8_2018UL/sim_HToAATo4Tau_M_3p7/240407_130328/0000/AODSIM_HToAATo4Tau_777.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/crab_GEN_SIM_HToAATo4Tau_hadronic_tauDecay_M3p7/240927_094611/0000/GEN_SIM_HToAATo4Tau_M3p7_40.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Ele_M0p2_Run3_2023/crab_GEN_SIM_HToAATo4Ele_M0p2/241007_073344/0000/GEN_SIM_HToAATo4Ele_0p2GeV_52.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Ele_M4_Run3_2023/crab_GEN_SIM_HToAATo4Ele_M4/241007_073556/0000/GEN_SIM_HToAATo4Ele_4GeV_75.root'
#inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/CMSSW_13_0_17/src/MCProduction/E2E-HToAATo4Ele/GEN_SIM_HToAATo4Ele_pdgID25_0p2GeV.root'
inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/check/CMSSW_13_0_14/src/MCProduction/E2E-HToAATo4Ele/HToAATo4Ele_1p2GeV.root'
maxEvents_=-1
skipEvents_=0#
outputFile_='GenInfo_only_H2AA4Ele_1p2GeV.root'
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(cmd)
os.system(cmd)
