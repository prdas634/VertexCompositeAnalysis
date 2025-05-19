import FWCore.ParameterSet.Config as cms
#import HLTrigger.HLTfilters.hltHighLevel_cfi
process = cms.Process("ANALYSIS")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#                                     TryToContinue = cms.untracked.vstring('ProductNotFound'))

# Define the input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('/store/data/Run2018C/JetHT/MINIAOD/UL2018_MiniAODv2-v1/00000/9F030B27-DBB0-EE46-A8FB-64FE0C417EE9.root')
#/store/data/Run2018D/JetHT/MINIAOD/UL2018_MiniAODv2-v2/2560000/0B3138BB-D2C9-AB49-B79F-32711AE032F4.root')
#/store/data/Run2018D/JetHT/MINIAOD/UL2018_MiniAODv2-v2/2560000/0337EBE4-98DB-AA47-A434-604D39F0CF02.root')
#/store/data/Run2018D/JetHT/MINIAOD/UL2018_MiniAODv2-v2/25530000/2E63C1BC-8EE1-D146-A831-1193C4E90C14.root')
#/store/data/Run2018C/JetHT/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/015CDDF7-DBCB-094D-861E-54F73012741A.root')
#        fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/flowcorr/JetHT/Ak8Jet500Skim_JetHT_Run2018D/210123_050442/0000/ppRun2UL_MINIAOD_17.root'),
#        fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/flowcorr/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/Ak8Jet500Skim_QCDPt470_Pythia8_UL18/210128_140537/0000/ppRun2UL_MINIAOD_10.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('106X_dataRun2_v32')

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.lambdaana = particleAna.clone(
  recoParticles = cms.InputTag("generalLambdaCandidatesNew"),
  triggerInfo = cms.untracked.VPSet([
    cms.PSet(path = cms.string('HLT_AK8PFJet500_v*')),
  ]),
  selectEvents = cms.string("eventFilter_step"),
)

process.antilambdaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalAntiLambdaCandidatesNew")
)

process.kshortana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalKshortCandidatesNew")
)

process.xiana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalXiCandidatesNew")
)

process.antixiana = process.xiana.clone(
  recoParticles = cms.InputTag("generalAntiXiCandidatesNew")
)

process.omegaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalOmegaCandidatesNew")
)

process.antiomegaana = process.omegaana.clone(
  recoParticles = cms.InputTag("generalAntiOmegaCandidatesNew")
)

process.generalanaNewSeq = cms.Sequence(process.lambdaana * process.antilambdaana * process.kshortana * process.xiana * process.antixiana * process.omegaana * process.antiomegaana)
process.generalana_step = cms.EndPath( process.generalanaNewSeq )


# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.HLTPaths = [
    'HLT_AK8PFJet500*'
    #'HLT_AK8PFJet500_v*',
    ]
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hlt500= process.hltFilter.clone()

process.eventFilterHLT500 = cms.Sequence(process.hlt500)

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter
)
process.eventFilter_step = cms.Path( process.eventFilter )

# Define the analysis steps
process.v0rereco_step = cms.Path(process.eventFilter
                                 * process.generalLambdaCandidatesNew
                                 * process.generalAntiLambdaCandidatesNew
                                 * process.generalKshortCandidatesNew
                                 * process.generalXiCandidatesNew
                                 * process.generalAntiXiCandidatesNew
                                 * process.generalOmegaCandidatesNew
                                 * process.generalAntiOmegaCandidatesNew
                                 )

#process.generalLambdaCandidatesNew.fitAlgo = cms.vuint32([3])
#process.generalKshortCandidatesNew.fitAlgo = cms.vuint32([3])
#process.generalXiCandidatesNew.fitAlgo = cms.vuint32([3])
#process.generalOmegaCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalLambdaCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalAntiLambdaCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalKshortCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalXiCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalAntiXiCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalOmegaCandidatesNew.fitAlgo = cms.vuint32([0])
process.generalAntiOmegaCandidatesNew.fitAlgo = cms.vuint32([0])


process.analyzerOffline = cms.EDAnalyzer('TrackAnalyzer',
                                         vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                         packedCandSrc = cms.InputTag("packedPFCandidates"),
                                         lostTracksSrc = cms.InputTag("lostTracks"),
                                         jets2 = cms.InputTag('slimmedJetsAK8'),
                                         #pfjetH = cms.InputTag('hltAK8PFJetsCorrectedMatchedToCaloJets200')
                                         )


# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('v0ana_data_wJet_testAll.root'))

#process.output = cms.OutputModule("PoolOutputModule",
 #                                 fileName = cms.untracked.string('ppRun2UL_SKIM_MINIAOD_wJet.root'),
  #                                SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
  #dataset = cms.untracked.PSet(
   #                                   dataTier = cms.untracked.string('MINIAOD')
    #                              )
     #                             )
#process.output.outputCommands = cms.untracked.vstring('drop *',
 #                                                     'keep *_*_*_ANASKIM1'
  #                                                    )

#process.output_step = cms.EndPath(process.output)

#main forest sequence
process.runAnalyzer500 = cms.Path(
    #process.analyzerOnline *
    process.eventFilterHLT500 *
    process.analyzerOffline
)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
 #   process.eventFilterHLT500 *,
#    process.analyzerOffline,
    process.runAnalyzer500,
    process.v0rereco_step,
    process.generalana_step,
  #  process.output_step
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
changeToMiniAOD(process)
