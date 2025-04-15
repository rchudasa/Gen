import os
cfg='GenAnalyzer/python/conFig_cfg.py'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M3p7/250108_164018/0000/GEN_SIM_HToAATo4Tau_M3p7_1.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M4_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M4/250108_164058/0000/GEN_SIM_HToAATo4Tau_M4_97.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M5_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M5/250108_164117/0000/GEN_SIM_HToAATo4Tau_M5_1.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M6_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M6/250108_164136/0000/GEN_SIM_HToAATo4Tau_M6_8.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M8_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M8/250108_221209/0000/GEN_SIM_HToAATo4Tau_M8_4.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M10_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M10/250108_221723/0000/GEN_SIM_HToAATo4Tau_M10_761.root'
#inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M12_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M12/250108_222012/0000/GEN_SIM_HToAATo4Tau_M12_976.root'
inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/crab_GEN_SIM_HToAATo4Tau_tauDecay_M14/250108_222401/0000/GEN_SIM_HToAATo4Tau_M14_2.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/GEN_HToAAToTauTau_M13_2018UL.root'
maxEvents_=-1
skipEvents_=0#
outputFile_='GenInfo_only_H2AA4Tau_M14GeV.root'
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(cmd)
os.system(cmd)
