import os
cfg='GenAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m3p6To18_pt30To300_v2/AOD_ATo4Tau_Hadronic_m3p6To18/241103_222918/0000/AOD_ATo2Tau_extra_collection_1.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M14/250108_222401/0000/GEN_SIM_HToAATo4Tau_M14_2.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/GEN_HToAAToTauTau_M13_2018UL.root'
maxEvents_=10
skipEvents_=0#
outputFile_='GenInfo_only_ATo2Tau_fromAOD.root'
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(cmd)
os.system(cmd)
