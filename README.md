This repo is only to access the Gen information and Trigger

### 1. Set up CMS release
```
export SCRAM_ARCH=el8_amd64_gcc11

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsrel CMSSW_13_0_13

cd CMSSW_13_0_13/src

cmsenv

```
### 2. Clone GenAnalyzer

` git clone git@github.com:bhbam/Gen.git `
