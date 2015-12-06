import FWCore.ParameterSet.Config as cms

process = cms.Process("RecHitAnalysis")

# Include the RandomNumberGeneratorService definition
#process.load("IOMC.RandomEngine.IOMC_cff")

# Famos Common inputs 
#process.load("FastSimulation.Configuration.CommonInputs_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoTracker.TrackProducer.TrackRefitter_cfi") 
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')



process.RecHitAnalysis = cms.EDAnalyzer("RecHitAnalayzer",
    track_label = cms.InputTag( "TrackRefitter"),
    simhit_label = cms.VInputTag(
        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),     
        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTECHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"),     
        cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
    ),
    verbose = cms.int32(0),
    isFastSimOnly = cms.bool(False),
    outfile = cms.string("rechitanalyzer_SimHitCollections_Full.root"),
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
#    debugFlag = cms.untracked.bool(True),
#    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring(
	#'file:/afs/cern.ch/work/a/ajafari/FastSim/385751AF-783D-E511-918C-0025905B8598.root'
	#'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-RECO/75X_mcRun2_asymptotic_v3-v1/00000/3AEBD3AF-783D-E511-8D7D-0025905A48F2.root',
	'file:./FullSimInput/step3_RAW2DIGI_L1Reco_RECO_EI_PAT_VALIDATION_DQM.root',
    )
)
process.seq = cms.Sequence(process.MeasurementTrackerEvent * process.TrackRefitter * process.RecHitAnalysis)
if process.RecHitAnalysis.isFastSimOnly:
    print "I am fast simulation"
    process.load('FastSimulation.Configuration.Reconstruction_BefMix_cff')
    process.RecHitAnalysis.simhit_label =  cms.VInputTag(cms.InputTag( "famosSimHits","TrackerHits"))
    process.RecHitAnalysis.track_label = cms.InputTag( "generalTracks"),
    process.RecHitAnalysis.outfile = cms.string("rechitanalyzer_SimHitCollections_Fast.root")
    process.source.fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/75X_mcRun2_asymptotic_v3_FastSim-v1/00000/0879F22D-763D-E511-BF50-0025905938B4.root')
    process.seq  = cms.Sequence(process.RecHitAnalysis)
process.Path = cms.Path(process.seq)


