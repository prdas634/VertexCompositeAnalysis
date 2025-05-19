from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'TestPD_Data_ppRun2UL_MiniAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['v0ana_wJet_test.root']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4000
#config.JobType.inputFiles = ['HeavyIonRPVRcd_PbPb2018_offline.db']
config.section_('Data')

config.Data.inputDataset = '/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt

config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 150
config.Data.publication = False

config.Data.outLFNDirBase = '/store/group/phys_heavyions/prdas/Run2_v0_wjet_trees_DataTest/'

#config.Data.totalUnits        = 1 # root file or lumi section for test only

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
