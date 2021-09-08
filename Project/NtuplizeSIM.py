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
TestNEvent = 100
IsPrivate = False
IsLocal = False
IsFullRun = True
if ifile < 0:
  IsLocal = True
if ifile == -1:
  IsFullRun = False
if IsPrivate:
  IsRun4 = True
# ifile: >=0 number of file to process. -1: process 100 event local file. -2: process all local file.
# print("process number: ", ifile)

inputFile = ""
outputname = ""
if IsLocal:
  inputFile = 'file:' + ('Private' if IsPrivate else 'step1') + ('Run4' if IsRun4 else 'Run3') + '.root'
  outputname = 'out/' + ('Private_' if IsPrivate else 'out_') + ('Run4' if IsRun4 else 'Run3') + '.root'
else:
  with open("/afs/cern.ch/user/s/siluo/Work/Muon/filenames/"+datatag+".txt") as filenames:
    for i, line in enumerate(filenames):
      if i == options.ifile:
        inputFile = line
        inputFile = inputFile.strip('\n')
        break
  outputname = '/eos/user/s/siluo/Muon/'+outputtag+'/InProgress/out_'+str(ifile)+'.root'

if ifile >=  0: print("Processing {}th file of dataset {}.".format(ifile,datatag))
if ifile == -1: print("Testing " + str(TestNEvent) + " events of " + inputFile)
if ifile == -2: print("Testing all events of " + inputFile)

import FWCore.ParameterSet.Config as cms

if IsRun4:
  from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
  process = cms.Process('ReL1',Phase2C11I13M9)
else:
  from Configuration.Eras.Era_Run3_cff import Run3
  process = cms.Process('ReL1',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
if IsRun4:
  process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
  process.load('Configuration.Geometry.GeometryExtended2026D76_cff')
else:
  process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if not IsLocal:
  process.MessageLogger.cerr.FwkReport.reportEvery = 100

# process.MessageLogger.suppressWarning = cms.untracked.vstring('GEMRawToDigiModule')
# process.MessageLogger.suppressWarning = cms.untracked.vstring("muonGEMDigis","simEmtfDigis")
# process.MessageLogger.suppressError = cms.untracked.vstring("simCscTriggerPrimitiveDigisRun3CCLUTILT")
# process.MessageLogger.suppressError = cms.untracked.vstring("simCscTriggerPrimitiveDigisRun3CCLUT")
# process.MessageLogger.suppressError = cms.untracked.vstring("simCscTriggerPrimitiveDigis")
# process.MessageLogger.cerr.threshold = cms.untracked.string('ERROR')

process.maxEvents = cms.untracked.PSet(
  # input = cms.untracked.int32(-1),
  input = cms.untracked.int32(TestNEvent if not IsFullRun else -1),
  output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(inputFile),
  # fileNames = cms.untracked.vstring(('file:step1' + ('Run4' if IsRun4 else 'Run3') + '.root') if IsLocal else inputFile),
  secondaryFileNames = cms.untracked.vstring()
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
  # process.GlobalTag = GlobalTag(process.GlobalTag, '113X_mcRun4_realistic_v7', '')
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')


from GEMCode.GEMValidation.cscTriggerCustoms import addCSCTriggerRun3
process = addCSCTriggerRun3(process)

from GEMCode.GEMValidation.cscTriggerCustoms import runOn120XMC
process = runOn120XMC(process)
if IsRun4:
  process.muonGEMDigis.useDBEMap = False
  # process.muonGEMDigis.readMultiBX = True
process.simCscTriggerPrimitiveDigis.commonParam.runCCLUT = cms.bool(True)

# NtupleMaker
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputname), closeFileFast = cms.untracked.bool(True))

from GEMCode.GEMValidation.simTrackMatching_cfi import simTrackPSet
process.NtupleMaker = cms.EDAnalyzer('NtupleMaker',
                                       simTrackPSet,
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
ana.simTrack.minEta = 1.2
ana.simTrack.maxEta = 2.4
ana.simTrack.minPt = 3
ana.gemSimHit.verbose = 0
# ana.gemStripDigi.verbose = 0
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
# ana.cscLCT.addGhostLCTs = cms.bool(True)
# ana.cscLCT.addGhosts = cms.bool(True)
ana.cscLCT.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisILT","","ReL1")
ana.cscCLCT.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisILT","","ReL1")
ana.cscALCT.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisILT","","ReL1")
# ana.muon.inputTag = cms.InputTag("gmtStage2Digis","Muon")
ana.gemStripDigi = cms.PSet(
  verbose = cms.int32(0),
  # inputTag = cms.InputTag("muonGEMDigis"),
  inputTag = cms.InputTag("simMuonGEMDigis" if IsRun4 else "muonGEMDigis"),
  minBX = cms.int32(-1),
  maxBX = cms.int32(1),
  matchDeltaStrip = cms.int32(1),
  matchToSimLink = cms.bool(False)
)
# ana.gemCoPadDigi.inputTag = cms.InputTag("simCscTriggerPrimitiveDigisILT","")
process.ana = cms.Path(ana)


# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.muonGEMDigis_step = cms.Path(process.muonGEMDigis)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
# process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.muonGEMDigis_step,process.L1simulation_step,process.ana,process.endjob_step)
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
