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
process.load('FastSimulation.Configuration.Reconstruction_BefMix_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "75X_mcRun2_asymptotic_v3::All"
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#process.load("FastSimulation.Configuration.FamosSequences_cff")

# Magnetic Field (new mapping, 3.8 and 4.0T)
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

#process.load("FastSimulation.Configuration.mixNoPU_cfi")
#process.mix.playback = cms.untracked.bool(True)
# RecHit Analysis ###
#process.load("FastSimulation.TrackingRecHitProducer.FamosRecHitAnalysis_cfi")

process.RecHitAnalysis = cms.EDAnalyzer("RecHitAnalayzer",
    track_label = cms.InputTag( "generalTracks"),
    simhit_label = cms.InputTag( "famosSimHits","TrackerHits"),
    verbose = cms.int32(0),
    isFastSimOnly = cms.bool(True),
    outfile = cms.string("rechitanalyzer_fastLayer.root"),
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
#    debugFlag = cms.untracked.bool(True),
#    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring(
#	'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-RECO/75X_mcRun2_asymptotic_v3-v1/00000/385751AF-783D-E511-918C-0025905B8598.root'
#FASTSIM------#
	'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/75X_mcRun2_asymptotic_v3_FastSim-v1/00000/0879F22D-763D-E511-BF50-0025905938B4.root'
#FULLSIM------#
	#'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-RECO/75X_mcRun2_asymptotic_v3-v1/00000/385751AF-783D-E511-918C-0025905B8598.root',
	#'file:/afs/cern.ch/work/a/ajafari/FastSim/385751AF-783D-E511-918C-0025905B8598.root'
	#'/store/relval/CMSSW_7_5_1/RelValSingleMuPt100_UP15/GEN-SIM-RECO/75X_mcRun2_asymptotic_v3-v1/00000/3AEBD3AF-783D-E511-8D7D-0025905A48F2.root',
    )
)

process.Path = cms.Path(process.RecHitAnalysis)


