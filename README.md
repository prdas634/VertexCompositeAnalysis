# VertexCompositeAnalysis

Example of setting up and running V0 reconstruction inside jets in pp

cmsrel CMSSW_10_6_4_patch1

cd CMSSW_10_6_4_patch1/src

cmssw-el7

cmsenv

git clone -b master https://github.com/prdas634/VertexCompositeAnalysis.git

cd VertexCompositeAnalysis

scram b -j8

voms-proxy-init --voms cms

cd VertexCompositeProducer/test

## V0 reconstruction inside jets in pp using MINIAOD
cmsRun ppRun2UL_V0Both_MiniAOD_cfg_testPD.py
