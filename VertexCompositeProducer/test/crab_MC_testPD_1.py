from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'TestPD_MC_ppRun2UL_MiniAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['v0ana_mc_wJet_test.root']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4000
#config.JobType.inputFiles = ['HeavyIonRPVRcd_PbPb2018_offline.db']
config.section_('Data')
#config.Data.inputDataset = '/JetHT/Run2018D-12Nov2019_UL2018-v4/MINIAOD'
#config.Data.inputDataset = '/JetMET0/Run2023B-19Dec2023-v1/MINIAOD'
#config.Data.inputDataset = '/JetMET1/Run2023D-22Sep2023_v2-v1/MINIAOD'
#config.Data.inputDataset = '/JetMET/Run2022G-22Sep2023-v2/MINIAOD'
config.Data.inputDataset = '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/user/x/xiaoyul/Cert_Collisions2023_366442_370790_Golden_JSON.txt'
#.txt file convert from raw data of /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json
#config.Data.lumiMask = '/afs/cern.ch/user/x/xiaoyul/Cert_Collisions2022_355100_362760_Golden_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/work/p/prdas/private/TestV0Reco/JSON_RunFiles/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'


config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 150
config.Data.publication = False

#config.Data.outLFNDirBase = '/store/user/pgardner/MINIAOD_2018_UL_D_ak8_new'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/Run3_jet_trees/'

config.Data.outLFNDirBase = '/store/group/phys_heavyions/prdas/Run2_v1_mc_wjet_trees_testMC/'

#config.Data.totalUnits        = 1 # root file or lumi section for test only

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
