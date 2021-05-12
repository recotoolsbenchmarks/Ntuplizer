import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('Validator',
                       debug              = cms.bool(False),
                       extendFormat       = cms.bool(False),
                       applyjec           = cms.bool(False),
                       vertices           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                       vertices4D         = cms.InputTag("offlineSlimmedPrimaryVertices4D"),
                       pfCandid           = cms.InputTag("packedPFCandidates"),
                       pileUp             = cms.InputTag("slimmedAddPileupInfo"),
                       genParts           = cms.InputTag("prunedGenParticles"),
                       genJets            = cms.InputTag("slimmedGenJets"),
                       genMet             = cms.InputTag("genMetTrue"),
                       photonsEB          = cms.InputTag("slimmedPhotons"),
                       photonsEE          = cms.InputTag("photonsFromMultiCl"),
                       photons_EBmva      = cms.InputTag("photonMVAIDProducerEB"),
                       photons_EEmva      = cms.InputTag("photonMVAIDProducerHGCal"),
                       electronsEB        = cms.InputTag("slimmedElectrons"),
                       electronsEE        = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
                       electrons_EBLoose  = cms.InputTag("electronMVAIDProducerEB","LooseWP"),
                       electrons_EBMedium = cms.InputTag("electronMVAIDProducerEB","MediumWP"),
                       electrons_EBTight  = cms.InputTag("electronMVAIDProducerEB","TightWP"),
                       electrons_EEmva    = cms.InputTag("electronMVAIDProducerHGCal"),
                       muons              = cms.InputTag("slimmedMuons" ),
                       taus               = cms.InputTag("slimmedTausNewID"),
                       jets               = cms.InputTag("slimmedJetsPuppi","","ReRECO"),
                       jetschs            = cms.InputTag("slimmedJets"),
                       fatjets            = cms.InputTag("slimmedJetsAK8"),
                       met                = cms.InputTag("slimmedMETsPuppi","","ReRECO"),
                       metpf              = cms.InputTag("slimmedMETs"),
                       
)
