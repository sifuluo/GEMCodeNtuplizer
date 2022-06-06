# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step2bis.py --filein file:step3.root --fileout file:step2bis.root --mc --eventcontent FEVTDEBUG --datatier GEN-SIM-DIGI-L1 --conditions auto:phase1_2021_realistic --step L1 --geometry DB:Extended --era Run3 --python_filename step2bis_L1.py --no_exec -n 10


import os
import sys
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('ifile', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "input file index")
options.register('dataset', "CMSSW12/Run3/RVSMPt10", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "input dataset name")
options.register('outputtag', "ME1/RVSMPt1000noPU", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "output folder name")
options.parseArguments()
ifile = options.ifile
datatag = options.dataset
outputtag = options.outputtag
IsRun4 = True
TestNEvent = -1 # Event number for local test
IsLocal = False
if not IsLocal: TestNEvent = -1 #precaution for batch job only run incompletely
if ifile < 0: ifile = 0 #precaution for testing local private sample

ChangeAlgo = False
B_matchWithHS = True
B_assignGEMCSCBending = False
B_migitageSlopeByCosi = True
unpack = False
inputFile = ""
outputname = ""
if IsLocal:
  inputFile = 'file:step2_' + str(ifile) + '.root'
  outputname = 'out/RelVal_' + str(ifile) + '.root'
  # inputFile = 'file:step2pu_0.root'
  # outputname = 'out/out_step2.root'

else:
  with open("/afs/cern.ch/user/s/siluo/Work/Muon/filenames/"+datatag+".txt") as filenames:
    for i, line in enumerate(filenames):
      if i == options.ifile:
        inputFile = line
        inputFile = inputFile.strip('\n')
        break
  outputname = '/eos/user/s/siluo/'+outputtag+'/InProgress/out_'+str(ifile)+'.root'

if IsLocal: print("Testing {} events of file {}.".format(TestNEvent, inputFile))
else: print("Processing {}th file of dataset {}.".format(ifile, datatag))

import FWCore.ParameterSet.Config as cms

# if IsRun4:
#   from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
#   process = cms.Process('ReL1',Phase2C11I13M9)
# else:
#   from Configuration.Eras.Era_Run3_cff import Run3
#   process = cms.Process('ReL1',Run3)

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('ReL1',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
# if IsRun4:
#   process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
#   # process.load('Configuration.Geometry.GeometryExtended2026D76_cff')
# else:
#   process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if not IsLocal:
  process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(
  # input = cms.untracked.int32(-1),
  input = cms.untracked.int32(TestNEvent),
  output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(inputFile),
  secondaryFileNames = cms.untracked.vstring(),
  # dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
  # inputCommands = cms.untracked.vstring(
  #       'keep *',
  #       'drop *_simMuonGEMDigis_*_*'
  #   )
)

process.options = cms.untracked.PSet(
  FailPath = cms.untracked.vstring(),
  IgnoreCompletely = cms.untracked.vstring(),
  Rethrow = cms.untracked.vstring(),
  SkipEvent = cms.untracked.vstring(),
  allowUnscheduled = cms.obsolete.untracked.bool,
  canDeleteEarly = cms.untracked.vstring(),
  emptyRunLumiMode = cms.obsolete.untracked.string,
  eventSetup = cms.untracked.PSet(
    forceNumberOfConcurrentIOVs = cms.untracked.PSet(

    ),
    numberOfConcurrentIOVs = cms.untracked.uint32(1)
  ),
  fileMode = cms.untracked.string('FULLMERGE'),
  forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
  makeTriggerResults = cms.obsolete.untracked.bool,
  numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
  numberOfConcurrentRuns = cms.untracked.uint32(1),
  numberOfStreams = cms.untracked.uint32(0),
  numberOfThreads = cms.untracked.uint32(1),
  printDependencies = cms.untracked.bool(False),
  sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
  throwIfIllegalParameter = cms.untracked.bool(True),
  wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
  annotation = cms.untracked.string('step2bis.py nevts:10'),
  name = cms.untracked.string('Applications'),
  version = cms.untracked.string('$Revision: 1.19 $')
)

# # Output definition
#
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-L1'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:out/out.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if IsRun4:
  # process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
  process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun3_2021_realistic_v14', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

if ChangeAlgo:
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE11.matchWithHS = cms.bool(B_matchWithHS)
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE21.matchWithHS = cms.bool(B_matchWithHS)
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE11.assignGEMCSCBending = cms.bool(B_assignGEMCSCBending)
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE21.assignGEMCSCBending = cms.bool(B_assignGEMCSCBending)
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE11.mitigateSlopeByCosi = cms.bool(B_migitageSlopeByCosi)
  process.simCscTriggerPrimitiveDigis.tmbPhase2GE21.mitigateSlopeByCosi = cms.bool(B_migitageSlopeByCosi)

# process.simCscTriggerPrimitiveDigis.tmbPhase2GE21.enableMatchGEMandME1a = cms.bool(False)
# process.simCscTriggerPrimitiveDigis.tmbPhase2GE21.enableMatchGEMandME1b = cms.bool(False)
# process.simCscTriggerPrimitiveDigis.commonParam.runME21ILT = cms.bool(False)

## Run-3 patterns with CCLUT, with ILT
process.simCscTriggerPrimitiveDigisReL1 = process.simCscTriggerPrimitiveDigis.clone()
# process.simCscTriggerPrimitiveDigisReL1.commonParam.run3 = cms.bool(True)
# process.simCscTriggerPrimitiveDigisReL1.commonParam.runCCLUT_OTMB = cms.bool(True)
# process.simCscTriggerPrimitiveDigisReL1.commonParam.runCCLUT_TMB = cms.bool(True)
# process.simCscTriggerPrimitiveDigisReL1.commonParam.runME11ILT = True
# process.simCscTriggerPrimitiveDigisReL1.commonParam.GEMPadDigiClusterProducer = cms.InputTag("simMuonGEMPadDigiClusters")
process.simEmtfDigisReL1 = process.simEmtfDigis.clone()
process.simEmtfDigisReL1.CSCInput = cms.InputTag('simCscTriggerPrimitiveDigisReL1','MPCSORTED',process._Process__name)
#
# ## redefine the L1-step
process.SimL1Emulator = cms.Sequence(
    process.simMuonGEMPadDigis *
    process.simMuonGEMPadDigiClusters *
    process.simCscTriggerPrimitiveDigis *
    process.simCscTriggerPrimitiveDigisReL1 *
    process.simEmtfDigis *
    process.simEmtfDigisReL1
)
#
# if not IsRun4:
#   process.muonGEMDigis.useDBEMap = True
#   process.simMuonGEMPadDigis.InputCollection = "muonGEMDigis"

# if IsRun4:
#   process.muonGEMDigis.useDBEMap = False
#   process.simMuonGEMPadDigis.InputCollection = "simMuonGEMDigis"
#   # process.muonGEMDigis.readMultiBX = True
# process.simCscTriggerPrimitiveDigis.commonParam.runCCLUT = cms.bool(True)

# NtupleMaker
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputname), closeFileFast = cms.untracked.bool(True))

from GEMCode.GEMValidation.simTrackMatching_cfi import simTrackPSet
from L1Trigger.CSCTriggerPrimitives.params.gemcscParams import gemcscPSets
process.NtupleMaker = cms.EDAnalyzer('NtupleMaker',
                                       simTrackPSet,
                                       gemcscPSets = gemcscPSets,
                                       verbose = cms.untracked.int32(0),
                                       cscStations = cms.vstring("CSC_ALL","CSC_ME11", "CSC_ME1a", "CSC_ME1b", "CSC_ME12", "CSC_ME13",
                                                                 "CSC_ME21", "CSC_ME22", "CSC_ME31", "CSC_ME32", "CSC_ME41", "CSC_ME42"),
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30000.0),      # only save TPs with |z0| < X cm
                                       IsRun4 = cms.bool(IsRun4),
                                       Print_matchCscStubs = cms.bool(False),
                                       Print_allCscStubs = cms.bool(False),
                                       Print_all = cms.bool(False),
                                       Print_ALCT = cms.bool(False),
                                       Print_CLCT = cms.bool(False),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       useGEMs = cms.bool(True),
                                       )
ana = process.NtupleMaker
# ana.simTrack.minEta = 1.2
# ana.simTrack.maxEta = 2.4
# ana.simTrack.minPt = 3
# ana.gemSimHit.verbose = 0
# # ana.gemStripDigi.verbose = 0
# ana.gemStripDigi.matchDeltaStrip = 2
# ana.gemPadDigi.verbose = 0
# ana.gemCoPadDigi.verbose = 0
# ana.gemPadCluster.verbose = 0
# ana.cscComparatorDigi.verbose = 0
# ana.cscWireDigi.verbose = 0
# ana.cscALCT.verbose = 0
# ana.cscALCT.minBX = 2
# ana.cscALCT.maxBX = 4
# ana.cscCLCT.verbose = 0
# ana.cscCLCT.minBX = 6
# ana.cscCLCT.maxBX = 8
# ana.cscLCT.verbose = 0
# ana.cscLCT.addGhostLCTs = cms.bool(True)
# ana.cscLCT.addGhosts = cms.bool(True)
ana.cscLCT.inputTag   = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.cscCLCT.inputTag  = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.cscALCT.inputTag  = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.cscMPLCT.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisReL1","MPCSORTED","ReL1")
ana.cscShower.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.cscPreCLCT.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.cscCLCTPreTrigger.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
ana.gemCoPadDigi.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisReL1","","ReL1")
# ana.gemStripDigi = cms.PSet(
#   verbose = cms.int32(0),
#   inputTag = cms.InputTag("simMuonGEMDigis" if not unpack else "muonGEMDigis"),
#   minBX = cms.int32(-1),
#   maxBX = cms.int32(1),
#   matchDeltaStrip = cms.int32(1),
#   matchToSimLink = cms.bool(False)
# )
process.ana = cms.Path(ana)

# process.tmbPhase2GE11.matchWithHS = cms.bool(False)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.muonGEMDigis_step = cms.Path(process.muonGEMDigis)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
# process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

if unpack:
  process.schedule = cms.Schedule(process.muonGEMDigis_step,process.L1simulation_step,process.ana,process.endjob_step)
else:
  process.schedule = cms.Schedule(process.L1simulation_step,process.ana,process.endjob_step)
  # process.schedule = cms.Schedule(process.ana,process.endjob_step)
# Schedule definition
# if IsRun4:
#   process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step,process.ana,process.endjob_step)
# else:
#   process.schedule = cms.Schedule(process.muonGEMDigis_step,process.L1simulation_step,process.ana,process.endjob_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
