import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from RecoTauTag.RecoTau.tools import runTauIdMVA

options = VarParsing ('python')

#$$
#options.register('pileup', 200,
#                 VarParsing.multiplicity.singleton,
#                 VarParsing.varType.int,
#                 "Specify the pileup in the sample (used for choosing B tag MVA thresholds)"
#)
#$$

options.register('debug', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "For Debug purposes"
)

options.register('rerunBtag', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Rerun the B tagging algorithms using new training"
)


options.register('GlobalTag', '110X_mcRun4_realistic_v3',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Specify the global tag for the release "
)


options.register('Analyzr', 'Validator',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Specify which ntuplzer you want to run "
)

options.parseArguments()

process = cms.Process("MyAna", Phase2C9)

# Geometry, GT, and other standard sequences
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.GlobalTag.globaltag = "111X_mcRun4_realistic_T15_v1"


# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
)


# Input

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) ) 

options.inputFiles = ['file:input.root']

#options.secondaryInputFiles = '/store/mc/PhaseIITDRSpring19DR/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/AODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/130000/773F8021-12B1-7E4B-9B0B-F0C188EC63B8.root'



#options.inputFiles = ['/store/relval/CMSSW_11_0_0_pre13/RelValTTbar_14TeV/MINIAODSIM/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v2/20000/5E63BB51-0E53-104E-9ED4-7B2D73B5C930.root']


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                options.inputFiles 
        ),
                            
)

process.source.inputCommands = cms.untracked.vstring("keep *")

# HGCAL EGamma ID
#process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")
#process.load("RecoEgamma.Phase2InterimID.phase2EgammaRECO_cff")
# analysis
moduleName = options.Analyzr  
process.myana = cms.EDAnalyzer(moduleName)
process.load("TreeMaker.Ntuplzr."+moduleName+"_cfi")

process.myana.debug = options.debug
process.myana.extendFormat = True

#$$
postfix='WithNewTraining'
patJetSource = 'selectedUpdatedPatJets'+postfix
#$$
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
# The updateJetCollection function uncorrect the jets from MiniAOD and 
# then recorrect them using the curren set of JEC in the event setup, recalculates
# btag discriminators from new training
# And the new name of the updated jet collection becomes selectedUpdatedPatJets+postfix
if options.rerunBtag:
    updateJetCollection(
        process,
        jetSource      = cms.InputTag('slimmedJetsPuppi'),
        pvSource       = cms.InputTag('offlineSlimmedPrimaryVertices'),
        svSource       = cms.InputTag('slimmedSecondaryVertices'),
        pfCandidates	= cms.InputTag('packedPFCandidates'),
        jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = [
            'pfDeepCSVJetTags:probb', 
            'pfDeepCSVJetTags:probbb',
            'pfDeepFlavourJetTags:probb',
            'pfDeepFlavourJetTags:probbb',
            'pfDeepFlavourJetTags:problepb'
        ],
        postfix = postfix
    )
    process.myana.jets = cms.InputTag(patJetSource)
    #print patJetSource

    postfix='CHSWithNewTraining'
    patJetSource = 'selectedUpdatedPatJets'+postfix
    updateJetCollection(
	process,
	jetSource      = cms.InputTag('slimmedJets'),
	pvSource       = cms.InputTag('offlineSlimmedPrimaryVertices'),
	svSource       = cms.InputTag('slimmedSecondaryVertices'),
	pfCandidates	= cms.InputTag('packedPFCandidates'),
	jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
	btagDiscriminators = [
	    'pfDeepCSVJetTags:probb', 
	    'pfDeepCSVJetTags:probbb',
	    'pfDeepFlavourJetTags:probb',
	    'pfDeepFlavourJetTags:probbb',
	    'pfDeepFlavourJetTags:problepb'
	],
	postfix = postfix
    )
    #print patJetSource
    process.myana.jetschs = cms.InputTag(patJetSource)


tauIdEmbedder = runTauIdMVA.TauIDEmbedder(
    process, cms, updatedTauName = "slimmedTausNewID",
    toKeep = ["2017v2", "newDM2017v2", "newDMPhase2v1", "deepTau2017v2p1",  "againstEle2018", "againstElePhase2v1"]
)
tauIdEmbedder.runTauID()
tauSrc_InputTag = cms.InputTag('slimmedTausNewID')# to be taken for any n-tuplizer

for key in options._register.keys():
    print "{:<20} : {}".format(key, getattr(options, key))


process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string("file.root")
                                   )

#$$
#Trick to make it work with rerunBtag
process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)
#$$

# run
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("output.root"),
#                               outputCommands = cms.untracked.vstring("drop *", "keep *_slimmedTausNewID_*_*"))    
#process.check  = cms.OutputModule("PoolOutputModule",
#                               compressionAlgorithm = cms.untracked.string('LZMA'),
#                               compressionLevel = cms.untracked.int32(4),
#                               eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#                               dataset = cms.untracked.PSet(                                                                                                                                                #dataTier = cms.untracked.string('AODSIM'),                                                                                                                      #filterName = cms.untracked.string('')                                                                                                                                  #),
#                               fileName = cms.untracked.string("out.root"),
#                               SelectEvents = cms.untracked.PSet(                                                                                              
#                                                 SelectEvents = cms.vstring("p")                                                                                #                                            ) 
#                              )
#


#$$

process.p = cms.Path(process.rerunMvaIsolationSequence*process.slimmedTausNewID*process.myana, process.tsk)
#process.e = cms.EndPath(process.out)#process.check)

#open('ntupleFileDump.py','w').write(process.dumpPython())


