import FWCore.ParameterSet.Config as cms
process = cms.Process("ANALYSIS")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer20UL18MiniAODv2/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/100000/021E9C3E-CEB9-4343-9B3A-12790E92E8D9.root'))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v16')

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.lambdaana = particleAna_mc.clone(
    recoParticles = cms.InputTag("generalLambdaCandidatesNew"),
    triggerInfo = cms.untracked.VPSet([
        cms.PSet(path = cms.string('HLT_AK8PFJet500_v*')),
    ]),
    selectEvents = cms.string("eventFilter_step"),
    genParticles = cms.untracked.InputTag("mergedGenParticles"),
    #  addSource    = cms.untracked.bool(False),
    autoFillPdgId = cms.untracked.bool(False),
    genPdgId     = cms.untracked.vuint32([3122]),
    #  saveTree = cms.untracked.bool(False)
)

process.antilambdaana = process.lambdaana.clone(
    recoParticles = cms.InputTag("generalAntiLambdaCandidatesNew")
)

process.kshortana = process.lambdaana.clone(
    recoParticles = cms.InputTag("generalKshortCandidatesNew"),
    genPdgId = cms.untracked.vuint32([310])
)

#process.kshortana.genPdgId = cms.untracked.vuint32([310])

process.xiana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalXiCandidatesNew"),
    genPdgId = cms.untracked.vuint32([3312])
)

#process.xiana.genPdgId = cms.untracked.vuint32([3312])

process.antixiana = process.xiana.clone(
  recoParticles = cms.InputTag("generalAntiXiCandidatesNew")
)

process.omegaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalOmegaCandidatesNew"),
    genPdgId = cms.untracked.vuint32([3334])
)

#process.omegaana.genPdgId = cms.untracked.vuint32([3334])

process.antiomegaana = process.omegaana.clone(
  recoParticles = cms.InputTag("generalAntiOmegaCandidatesNew")
)

process.generalanaNewSeq = cms.Sequence(process.lambdaana * process.antilambdaana * process.kshortana * process.xiana * process.antixiana * process.omegaana * process.antiomegaana)
process.generalana_step = cms.EndPath( process.generalanaNewSeq )

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_AK8PFJet500*'
    #'HLT_AK8PFJet500_v*',
]

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

process.generalLambdaCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalAntiLambdaCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalKshortCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalXiCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalAntiXiCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalOmegaCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalAntiOmegaCandidatesNew.fitAlgo = cms.vuint32([3])

process.analyzerOffline = cms.EDAnalyzer('TrackAnalyzer',
                                         vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                         packedCandSrc = cms.InputTag("packedPFCandidates"),
                                         lostTracksSrc = cms.InputTag("lostTracks"),
                                         jets2 = cms.InputTag('slimmedJetsAK8'),
                                         #pfjetH = cms.InputTag('hltAK8PFJetsCorrectedMatchedToCaloJets200')
                                         doGen = cms.untracked.bool(True),
                                         genEvtInfo = cms.InputTag("generator"),
                                         packedGen = cms.InputTag("packedGenParticles"),
                                         genJets = cms.InputTag("slimmedGenJetsAK8"),
                                         puSummaryInfo = cms.InputTag("slimmedAddPileupInfo")
                                     )

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('v0ana_mc_wJet_testAll.root'))

#process.output = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string('ppRun2UL_SKIM_MINIAOD.root'),
#  SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
# dataset = cms.untracked.PSet(
#  dataTier = cms.untracked.string('MINIAOD')
#)
#)
#process.output.outputCommands = cms.untracked.vstring('drop *',
#      'keep *_*_*_ANASKIM1'
#)

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
    process.runAnalyzer500,
    process.v0rereco_step,
    process.generalana_step,
    #    process.output_step
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAODMC
changeToMiniAODMC(process)
