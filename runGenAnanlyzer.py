import os
cfg='GenAnalyzer/python/ConfFile_cfg.py'
cfg='GenAnalyzer/python/conFig_cfg.py'
inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/CMSSW_13_0_17/src/MCProduction/E2E-HToAATo4Tau/GEN_SIM_HToAATo4Tau_M3p7_withGenFilter.root'
maxEvents_=-1
skipEvents_=0#
outputFile_='GenInfo_H2AA4Tau_M3p7GeV_withGenFilter.root'
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(cmd)
os.system(cmd)
