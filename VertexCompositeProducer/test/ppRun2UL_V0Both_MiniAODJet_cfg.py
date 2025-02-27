import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM1")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
#        fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it///store/data/Run2018C/JetHT/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/015CDDF7-DBCB-094D-861E-54F73012741A.root'),
        fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/flowcorr/JetHT/Ak8Jet500Skim_JetHT_Run2018D/210123_050442/0000/ppRun2UL_MINIAOD_17.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('106X_dataRun2_v32')

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.lambdaana = particleAna.clone(
  recoParticles = cms.InputTag("generalLambdaCandidatesNew"),
  triggerInfo = cms.untracked.VPSet([
    cms.PSet(path = cms.string('HLT_AK8PFJet500_v*')), 
  ]),
  selectEvents = cms.string("eventFilter_step"),
)

process.kshortana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalKshortCandidatesNew")
)

process.xiana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalXiCandidatesNew")
)

process.omegaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalOmegaCandidatesNew")
)

process.generalanaNewSeq = cms.Sequence(process.lambdaana * process.kshortana * process.xiana * process.omegaana)
process.generalana_step = cms.EndPath( process.generalanaNewSeq )


# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_AK8PFJet500_v*', 
    ]

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter 
)
process.eventFilter_step = cms.Path( process.eventFilter )

# Define the analysis steps
process.v0rereco_step = cms.Path(process.eventFilter 
                               * process.generalLambdaCandidatesNew
                               * process.generalKshortCandidatesNew
                               * process.generalXiCandidatesNew 
                               * process.generalOmegaCandidatesNew
                               )

process.generalLambdaCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalKshortCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalXiCandidatesNew.fitAlgo = cms.vuint32([3])
process.generalOmegaCandidatesNew.fitAlgo = cms.vuint32([3])

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('v0ana.root'))

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('ppRun2UL_SKIM_MINIAOD.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('MINIAOD')
    )
)
process.output.outputCommands = cms.untracked.vstring('drop *',
      'keep *_*_*_ANASKIM1'
)

process.output_step = cms.EndPath(process.output)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.v0rereco_step,
    process.generalana_step,
    process.output_step
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAODJet
changeToMiniAODJet(process)
