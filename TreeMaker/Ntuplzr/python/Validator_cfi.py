import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('Validator',
                       debug          = cms.bool(True),               
                       vertices       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                       pileUp         = cms.InputTag("slimmedAddPileupInfo"),
                       genParts       = cms.InputTag("prunedGenParticles"),
                       genJets        = cms.InputTag("slimmedGenJets"),
                       genMet         = cms.InputTag("genMetTrue"),
                       #photons        = cms.InputTag("slimmedPhotons"),
                       #electrons      = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
                       photons        = cms.InputTag("phase2Photons"),
                       electrons      = cms.InputTag("phase2Electrons"),
                       muons          = cms.InputTag("slimmedMuons" ),
                       taus           = cms.InputTag("slimmedTaus"),
                       jets           = cms.InputTag("slimmedJetsPuppi"),
                       met            = cms.InputTag("slimmedMETsPuppi"),
)
