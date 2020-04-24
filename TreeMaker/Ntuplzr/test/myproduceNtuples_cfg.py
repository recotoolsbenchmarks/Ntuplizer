import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
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
#$$ process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')

process.load('Configuration.Geometry.GeometryExtended2026D52Reco_cff')

#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, options.GlobalTag, '')

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


#options.inputFiles = ['/store/relval/CMSSW_11_1_0_pre6/RelValTTbar_14TeV/MINIAODSIM/111X_mcRun3_2021_realistic_v3-v1/20000/DAF5417B-BAE7-B843-A47F-8C728774661F.root']
options.inputFiles = ['/store/relval/CMSSW_11_1_0_pre1/RelValTTbar_14TeV/MINIAODSIM/PU25ns_110X_mcRun4_realistic_v2_2026D52PU200_ext1-v1/20000/DDF640CE-0514-3143-A1EC-27DA32BA2047.root']
#options.inputFiles = '/store/mc/PhaseIITDRSpring19MiniAOD/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/MINIAODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/130000/B155F8E9-1F3C-A741-8E86-5BEABD3AFC13.root'
#options.inputFiles = ['/store/mc/PhaseIITDRSpring19MiniAOD/WpWpJJ_EWK_TuneCP5_14TeV-madgraph-pythia8/MINIAODSIM/PU200_106X_upgrade2023_realistic_v3-v2/110000/863B8ABC-8FBA-0742-8092-71B024307821.root',
#                      '/store/mc/PhaseIITDRSpring19MiniAOD/WpWpJJ_EWK_TuneCP5_14TeV-madgraph-pythia8/MINIAODSIM/PU200_106X_upgrade2023_realistic_v3-v2/110000/91644A10-3DBF-9249-8D98-70061BC075E9.root',
#                      '/store/mc/PhaseIITDRSpring19MiniAOD/WpWpJJ_EWK_TuneCP5_14TeV-madgraph-pythia8/MINIAODSIM/PU200_106X_upgrade2023_realistic_v3-v2/110000/A885BB63-20F8-5E41-A0E5-3BE38BA5B0BF.root']

#options.inputFiles = ['/store/mc/PhaseIITDRSpring19MiniAOD/WpWpJJ_QCD_TuneCP5_14TeV-madgraph-pythia8/MINIAODSIM/PU200_106X_upgrade2023_realistic_v3-v2/250000/1BDEBB5D-01C8-6645-8F72-082549F9DD9B.root']
#options.inputFiles = ['/store/mc/PhaseIITDRSpring19MiniAOD/PhotonFlatPt8To150/MINIAODSIM/NoPU_106X_upgrade2023_realistic_v3-v1/230000/E8D9CEB8-64AC-C641-9E59-C488B6BE706F.root']
#options.inputFiles = ['/store/relval/CMSSW_10_6_0_patch2/RelValTTbar_14TeV/MINIAODSIM/PU25ns_106X_upgrade2023_realistic_v3_2023D41PU200-v1/10000/FFF6B0DE-6A3A-D04C-9A77-C195F92F8577.root']
#options.secondaryInputFiles = '/store/mc/PhaseIITDRSpring19DR/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/AODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/130000/773F8021-12B1-7E4B-9B0B-F0C188EC63B8.root'



#options.inputFiles = ['/store/relval/CMSSW_11_0_0_pre13/RelValTTbar_14TeV/MINIAODSIM/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v2/20000/5E63BB51-0E53-104E-9ED4-7B2D73B5C930.root']


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                options.inputFiles 
        ),
                            
)

process.source.inputCommands = cms.untracked.vstring("keep *")

# HGCAL EGamma ID
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")
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
#process.check  = cms.OutputModule("PoolOutputModule",
#                               compressionAlgorithm = cms.untracked.string('LZMA'),
#                               compressionLevel = cms.untracked.int32(4),
#                               eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#                               #dataset = cms.untracked.PSet(                                                                                                              #                                  #dataTier = cms.untracked.string('AODSIM'),                                                                                              #                                  #filterName = cms.untracked.string('')                                                                                                   #                               #),                                                                                                                                         #                               fileName = cms.untracked.string("out.root"),
#                               #SelectEvents = cms.untracked.PSet(                                                                                                         #                               #               SelectEvents = cms.vstring("p")                                                                                             #                               #               )                                                                                                                           #                               )
#


#$$
# process.p = cms.Path(process.phase2Egamma*process.myana)
process.p = cms.Path(process.phase2Egamma*process.myana, process.tsk)
#process.p = cms.Path(process.myana, process.tsk)
#process.e = cms.EndPath(process.check)
#process.p = cms.Path(process.myana)
#$$

#open('ntupleFileDump.py','w').write(process.dumpPython())


