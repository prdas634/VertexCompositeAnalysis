import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('RhoPbPb2018SKIM',eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.EventContent.EventContent_cff')


# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# DEBUG
#process.MessageLogger.cerr.threshold = "DEBUG"
#process.MessageLogger.debugModules = ["ntrkFilterRho"]

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
)

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
###'/store/hidata/HIRun2018A/HIMinimumBias1/AOD/PromptReco-v1/000/326/577/00000/532A6440-350D-2446-89FB-3CA9E8335E4F.root'
#'/store/hidata/HIRun2018A/HIMinimumBias6/AOD/04Apr2019-v1/70012/FEB2298F-D3AA-5A41-A999-02AEC7648569.root'
'file:output_numEvent100.root'
   ),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_v6')


# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","reRECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)


# Add PbPb event plane
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.hiEvtPlane.trackTag = cms.InputTag("generalTracks")
process.hiEvtPlane.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlane.useNtrk = cms.untracked.bool(False)
process.hiEvtPlane.caloCentRef = cms.double(-1)
process.hiEvtPlane.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRef = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2018_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
process.evtplane_seq = cms.Sequence(process.hiEvtPlane * process.hiEvtPlaneFlat)


# Rho candidate rereco
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")
process.generalRhoCandidatesNew = process.generalParticles.clone(
    pdgId = cms.int32(113),
    mass = cms.double(0.775),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(1.3),

    preSelection = cms.string(""
       "charge==0"
       "&& userFloat('dauPtSum') >= 1.6 && userFloat('dauEtaDiff') <= 1.0"
       ),
    pocaSelection = cms.string(""
       "pt >= 1.0"
       "&& userFloat('dca') >= 0 && userFloat('dca') <= 9999."
       ),
    postSelection = cms.string(""
       "userFloat('vertexProb') >= 0.02"
       "&& userFloat('normChi2') <= 9999.0"
       ),
    finalSelection = cms.string(""
       "userFloat('rVtxMag') >= 0.0 && userFloat('rVtxSig') >= 2.0"
       "&& userFloat('lVtxMag') >= 0.0 && userFloat('lVtxSig') >= 3.0"
       "&& cos(userFloat('angle3D')) >= -2.0 && cos(userFloat('angle2D')) >= -2.0"
       "&& abs(userFloat('angle3D')) <= 0.2 && abs(userFloat('angle2D')) <= 0.2"
       "&& abs(rapidity) < 2.0"
       ),
#
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
           selection = cms.string(
              "pt>1.0 && abs(eta)<2.4"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& (normalizedChi2/hitPattern.trackerLayersWithMeasurement)<0.18"
              "&& numberOfValidHits >=11"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) < 3.0 && abs(userFloat("dxySig")) < 3.0'
              '&& (track.algo!=6 || userFloat("mva")>=0.98)'
              )
           ),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
           selection = cms.string(
              "pt>1.0 && abs(eta)<2.4"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& (normalizedChi2/hitPattern.trackerLayersWithMeasurement)<0.18"
              "&& numberOfValidHits >=11"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) < 3.0 && abs(userFloat("dxySig")) < 3.0'
              '&& (track.algo!=6 || userFloat("mva")>=0.98)'
              )
           )
    ])
  )
process.generalRhoCandidatesNew.mva = cms.InputTag("generalTracks","MVAValues") ###cesar:to change iter6 tracking mva cut

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.generalanaNew = particleAna.clone(
    saveTree = False,
  recoParticles = cms.InputTag("generalRhoCandidatesNew"),
  triggerInfo = cms.untracked.VPSet([
    cms.PSet(path = cms.string('HLT_HIMinimumBias_*')), # Minimum bias
  ]),
  selectEvents = cms.string("ntrkFilterRho_step"),
)

process.generalanaNewSeq = cms.Sequence(process.generalanaNew)
process.generalana_step = cms.EndPath( process.generalanaNewSeq )

# Add PbPb collision event selection
process.load("VertexCompositeAnalysis.VertexCompositeProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')

process.load('VertexCompositeAnalysis.VertexCompositeProducer.ntrkUtils_cff')

# Define the event selection sequence

process.ntrkFilterRho = process.ntrkFilter.clone(
      useCent = cms.untracked.bool(True),
      nMultMin = cms.untracked.int32(0), # >=
      nMultMax = cms.untracked.int32(200), # <
      centBinMin = cms.untracked.int32(100), # >=
      centBinMax = cms.untracked.int32(200) # <
    )

process.ntrkFilterRho_seq = cms.Sequence(process.ntrkFilterRho)

process.eventFilter_HM = cms.Sequence(
    process.offlinePrimaryVerticesRecovery *
    process.hfCoincFilter2Th4 *
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )
process.ntrkFilterRho_step = cms.Path(process.eventFilter_HM * process.ntrkFilterRho_seq)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq * process.evtplane_seq * process.ntrkFilterRho_seq)
process.rhorereco_step = cms.Path(process.eventFilter_HM * process.ntrkFilterRho_seq* process.generalRhoCandidatesNew)

# Define the output
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('rhoana_ntuple_PbPb2018.root')
                                  )

process.output_HM = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PbPb2018_SKIM_AOD.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step', 'ntrkFilterRho_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)
process.output_HM.outputCommands = cms.untracked.vstring('drop *',
      'keep *_*_*_RhoPbPb2018SKIM',
      'keep *_offlinePrimaryVerticesRecovery_*_*',
      'keep *_hiEvtPlane_*_*',
      'keep *_centralityBin_*_*',
      'keep *_hiCentrality_*_RhoPbPb2018SKIM',
)


process.output_HM_step = cms.EndPath(process.output_HM)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.ntrkFilterRho_step,
    process.rhorereco_step,
    process.output_HM_step,
    process.generalana_step
)

# Add the event selection filters
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter)
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.clusterCompatibilityFilter)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter , process.Flag_clusterCompatibilityFilter ]
for P in eventFilterPaths:
      process.schedule.insert(0, P)

# peripheral pv recovery
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
