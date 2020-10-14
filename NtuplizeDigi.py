# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms
# from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_cff import Run3

import os
import sys
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('REPR',Run3)
# process = cms.Process('REPR',eras.Phase2C9)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)



options = VarParsing.VarParsing ('standard')
options.register('ifile', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "input file name")
options.parseArguments()
ifile = options.ifile
print("process number: ", ifile)
# resublist = [52,83,84,85,125,136,159,158,160,122,161,194,157,193,198,197,250,199,200,202,203,206,212,213,214]
# if ifile >=len(resublist):
    # quit()
# ifile = resublist[ifile]
# inputFile = ""
# with open("/afs/cern.ch/user/s/siluo/Work/Muon/filenames/DYToLL.txt") as filenames:
#     for i, line in enumerate(filenames):
#         if i == options.ifile:
#             inputFile = line
#             break
# inputFile = inputFile.strip('\n')

inputFile = 'file:step2bis'+str(ifile)+'.root';

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
# process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')

process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.Timing = cms.Service("Timing",
          #summaryOnly = cms.untracked.bool(False),
          #useJobReport = cms.untracked.bool(True)
#)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root'),
    #fileNames = cms.untracked.vstring('root://eoscms/eos/cms/store/relval/CMSSW_10_6_0_pre3/RelValMinBias_14TeV/GEN-SIM-DIGI-RAW/105X_upgrade2023_realistic_v5_2023D41noPU-v1/10000/E6CBA1C6-7A2E-A540-97B3-DE2C30AB70C8.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_pre3/RelValMuGunPt2To100/GEN-SIM-DIGI-RAW/105X_upgrade2023_realistic_v5_2023D41noPU-v2/10000/602E0B41-B698-6340-AC68-517578FEC457.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_106X_upgrade2023_realistic_v2_2023D41PU200-v1/10000/FEA5D564-937A-8D4B-9C9A-696EFC05AB58.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_patch2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/BC7B5A96-E3D2-ED48-81FC-35EF57134127.root'),
    # fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/DarkSUSY_mH_125_mGammaD_20_cT_1_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/20000/567B5C3B-984E-F94E-99B6-E63DE07E1C60.root'),
    # fileNames = cms.untracked.vstring('/store/mc/PhaseIITDRSpring19DR/DarkSUSY_mH_125_mGammaD_20_cT_1_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/20000/8124D852-95B1-C344-891C-20163A2E5A95.root'),
    # fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v1/60000/E0D5C6A5-B855-D14F-9124-0B2C9B28D0EA.root'),
    fileNames = cms.untracked.vstring(inputFile),
    # fileNames = cms.untracked.vstring('file:step2bis.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Ntuplize nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
#     outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

process.TFileService = cms.Service("TFileService", fileName = cms.string('out/out_'+str(ifile)+'.root'), closeFileFast = cms.untracked.bool(True))

from GEMCode.GEMValidation.simTrackMatching_cfi import simTrackPSet
process.NtupleMaker = cms.EDAnalyzer('NtupleMaker',
                                       simTrackPSet,
                                       verbose = cms.untracked.int32(1),
                                       cscStations = cms.vstring("CSC_ALL","CSC_ME11", "CSC_ME1a", "CSC_ME1b", "CSC_ME12", "CSC_ME13",
                                                                 "CSC_ME21", "CSC_ME22", "CSC_ME31", "CSC_ME32", "CSC_ME41", "CSC_ME42"),
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30000.0),      # only save TPs with |z0| < X cm
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       )
ana = process.NtupleMaker
ana.simTrack.minEta = 1.2
ana.simTrack.maxEta = 2.4
ana.simTrack.minPt = 3
ana.gemSimHit.verbose = 0
ana.gemStripDigi.verbose = 0
ana.gemStripDigi.matchDeltaStrip = 2
ana.gemPadDigi.verbose = 0
ana.gemCoPadDigi.verbose = 0
ana.gemPadCluster.verbose = 0
ana.cscComparatorDigi.verbose = 0
ana.cscWireDigi.verbose = 0
ana.cscALCT.verbose = 0
ana.cscALCT.minBX = 2
ana.cscALCT.maxBX = 4
ana.cscCLCT.verbose = 0
ana.cscCLCT.minBX = 6
ana.cscCLCT.maxBX = 8
ana.cscLCT.verbose = 0
ana.muon.inputTag = cms.InputTag("gmtStage2Digis","Muon")
ana.gemStripDigi.inputTag = cms.InputTag("muonGEMDigis")
process.ana = cms.Path(ana)

process.endjob_step = cms.EndPath(process.endOfProcess)
# process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth,process.ana)
# process.schedule = cms.Schedule(process.L1simulation_step,process.TTTracksEmulationWithTruth,process.ana)
# process.schedule = cms.Schedule(process.L1simulation_step,process.TTTracksEmulationWithTruth,process.CSCRecHit, process.CSCStub, process.GEMRecHit, process.ana)
# process.schedule = cms.Schedule(process.L1simulation_step,process.ana)
process.schedule = cms.Schedule(process.ana)


# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
# process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
