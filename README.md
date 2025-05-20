# VertexCompositeAnalysis

Example of setting up and running V0 reconstruction inside jets in pp

cmssw-el7
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_10_6_4_patch1
cd CMSSW_10_6_4_patch1/src
cmsenv
git clone -b master git@github.com:prdas634/VertexCompositeAnalysis.git
cd VertexCompositeAnalysis
scram b -j8
voms-proxy-init --voms cms
cd VertexCompositeProducer/test/

## V0 reconstruction inside jets in pp using MINIAOD

For data:
cmsRun TestPD_Data_ppRun2UL_MiniAOD.py 

For MC:
cmsRun TestPD_MC_ppRun2UL_MiniAOD.py