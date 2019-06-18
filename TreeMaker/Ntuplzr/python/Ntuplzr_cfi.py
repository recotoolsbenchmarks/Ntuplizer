import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('Ntuplzr',
                       debug          = cms.bool(True),               
                       genJets        = cms.InputTag("slimmedGenJets"),
                       genParts       = cms.InputTag("prunedGenParticles"),
                       #photons        = cms.InputTag("slimmedPhotons"),
                       #electrons      = cms.InputTag("slimmedElectrons"),
                       photons        = cms.InputTag("phase2Photons"),
                       ecalRecHits    = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                       electrons      = cms.InputTag("phase2Electrons"),
                       muons          = cms.InputTag("slimmedMuons" ),
                       jets           = cms.InputTag("slimmedJetsPuppi"),
                       met            = cms.InputTag("slimmedMETsPuppi"),
                       taus           = cms.InputTag("slimmedTaus"),
                       vertices       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                       beamspot       = cms.InputTag("offlineBeamSpot"),
                       conversions    = cms.InputTag("reducedEgamma", "reducedConversions","PAT")
)
