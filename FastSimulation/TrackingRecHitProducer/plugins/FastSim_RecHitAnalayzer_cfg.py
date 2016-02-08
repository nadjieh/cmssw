import FWCore.ParameterSet.Config as cms

process = cms.Process("RecHitAnalysis")


# Famos Common inputs 
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
process.load('FastSimulation.Configuration.Reconstruction_BefMix_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')



process.RecHitAnalysis = cms.EDAnalyzer("RecHitAnalayzer",
    track_label = cms.InputTag( "generalTracks"),
    simhit_label =  cms.VInputTag(cms.InputTag( "famosSimHits","TrackerHits")),
    verbose = cms.int32(0),
    isFastSimOnly = cms.bool(True),
    outfile = cms.string("rechitanalyzer_SimHitCollections_Fast.root")
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(9000)
)
process.source = cms.Source("PoolSource",
#    debugFlag = cms.untracked.bool(True),
#    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/75X_mcRun2_asymptotic_v3_FastSim-v1/00000/0879F22D-763D-E511-BF50-0025905938B4.root'
        '/store/relval/CMSSW_8_0_0_pre2/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/76X_mcRun2_asymptotic_v12_FastSim-v1/00000/10174562-EE8D-E511-ADE6-0025905B85D8.root',
        '/store/relval/CMSSW_8_0_0_pre2/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/76X_mcRun2_asymptotic_v12_FastSim-v1/00000/74274DA5-F28D-E511-865F-002618943982.root',
        #'/store/relval/CMSSW_8_0_0_pre2/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/76X_mcRun2_asymptotic_v12_FastSim-v1/00000/B6CBB93E-F28D-E511-83BC-0CC47A4D769C.root',
        #'/store/relval/CMSSW_8_0_0_pre2/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/76X_mcRun2_asymptotic_v12_FastSim-v1/00000/C2CC30B6-F18D-E511-BE59-002618943982.root',
        #'/store/relval/CMSSW_8_0_0_pre2/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/76X_mcRun2_asymptotic_v12_FastSim-v1/00000/C61EC6B6-F18D-E511-AC4E-0CC47A4C8E66.root',
    )
)

process.seq  = cms.Sequence(process.RecHitAnalysis)
process.Path = cms.Path(process.seq)


