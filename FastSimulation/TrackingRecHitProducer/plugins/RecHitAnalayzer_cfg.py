import FWCore.ParameterSet.Config as cms

process = cms.Process("RecHitAnalysis")

# Include the RandomNumberGeneratorService definition
#process.load("IOMC.RandomEngine.IOMC_cff")

# Famos Common inputs 
#process.load("FastSimulation.Configuration.CommonInputs_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MCRUN2_74_V9::All"
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
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
#    debugFlag = cms.untracked.bool(True),
#    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_3/RelValSingleMuPt100_UP15/GEN-SIM-DIGI-RECO/MCRUN2_74_V9_FastSim-v11/00000/220B5FCA-B606-E511-B79F-0025905A60B6.root')
)

process.Path = cms.Path(process.RecHitAnalysis)


