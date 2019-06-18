 // -*- C++ -*-
 //
 // Package:    FullSimAnalyser/Ntuplzr
 // Class:      Ntuplzr
 // 
 /**\class Ntuplzr Ntuplzr.cc FullSimAnalyser/Ntuplzr/plugins/Ntuplzr.cc

  Description: [one line class summary]

  Implementation:
      [Notes on implementation]
 */
 //
 // Original Author:  Sandhya Jain
 //         Created:  Wed, 29 Mar 2017 05:29:37 GMT
 //
 //


 // system include files
 #include <memory>
 #include <iostream>
 // user include files
 #include "FWCore/Framework/interface/Frameworkfwd.h"
 #include "FWCore/Framework/interface/one/EDAnalyzer.h"
 #include "FWCore/Framework/interface/Event.h"
 #include "FWCore/Framework/interface/MakerMacros.h"
 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/MessageLogger/interface/MessageLogger.h"//
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "CommonTools/UtilAlgos/interface/TFileService.h"

 #include "DataFormats/PatCandidates/interface/Muon.h"
 #include "DataFormats/MuonReco/interface/MuonSelectors.h"
 #include "DataFormats/PatCandidates/interface/Electron.h"
 #include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
 #include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
 #include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

 #include "DataFormats/PatCandidates/interface/Jet.h"
 #include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
 #include "DataFormats/PatCandidates/interface/MET.h"
 #include "DataFormats/PatCandidates/interface/PackedCandidate.h"
 #include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
 #include "DataFormats/JetReco/interface/GenJet.h"
 #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
 #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
 #include "DataFormats/PatCandidates/interface/Photon.h"
 #include "DataFormats/PatCandidates/interface/Tau.h"

 #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
 #include "TrackingTools/Records/interface/TransientTrackRecord.h"
 #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
 #include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
 #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
 #include "SimTracker/Records/interface/TrackAssociatorRecord.h"
 #include "DataFormats/VertexReco/interface/Vertex.h"
 #include "DataFormats/VertexReco/interface/VertexFwd.h"
 #include "Geometry/Records/interface/MuonGeometryRecord.h"
 #include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
 #include "Geometry/GEMGeometry/interface/ME0Geometry.h"

 #include "FWCore/Framework/interface/ESHandle.h"

 #include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
 #include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
 #include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
 #include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
 #include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
 #include "CommonTools/CandUtils/interface/AddFourMomenta.h"
 #include "DataFormats/Math/interface/deltaR.h"
 #include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
 //#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

 #include "TFile.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TTree.h"
 #include "TLorentzVector.h"
 #include "Math/GenVector/VectorUtil.h"

 using namespace std;
 using namespace reco;
 //
 // class declaration
 //

 // If the analyzer does not use TFileService, please remove
 // the template argument to the base class so the class inherits
 // from  edm::one::EDAnalyzer<> and also remove the line from
 // constructor "usesResource("TFileService");"
 // This will improve performance in multithreaded jobs.
 const Int_t kMaxVertices = 300;
 const Int_t kMaxWeights = 1500;
 const Int_t kMaxParticle = 10000;
 const Int_t kMaxGenJet = 800;
 const Int_t kMaxPhoton = 800;
 const Int_t kMaxElectron = 200;
 const Int_t kMaxMuonLoose = 200;
 const Int_t kMaxJet = 800;
 const Int_t kMaxTau = 200;
 const Int_t kMaxMissingET = 1;


 class Ntuplzr : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
 public:
   explicit Ntuplzr(const edm::ParameterSet&);
   ~Ntuplzr();

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
   virtual void beginJob() override;
   virtual void beginRun(edm::Run const&, edm::EventSetup const&) ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endRun(edm::Run const&, edm::EventSetup const&) ;
   virtual void endJob() override;

   //#bool isLooseElec    (const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot); 
   //#bool isMediumElec   (const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot); 
   //#bool isTightElec    (const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot); 
   bool isME0MuonSel   (reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
   bool isME0MuonSelNew(reco::Muon, double, double, double, edm::EventSetup const& );
   bool debug_;
   edm::Service<TFileService> fs_;
   edm::EDGetTokenT<std::vector<reco::GenJet>>      genJetsToken_      ;
   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_     ;
   edm::EDGetTokenT<GenEventInfoProduct>            generatorToken_    ;
   edm::EDGetTokenT<LHEEventProduct>                generatorlheToken_ ;
   edm::EDGetTokenT<std::vector<pat::Photon>>       photnsToken_       ; 
   edm::EDGetTokenT<EcalRecHitCollection>           ecalRecHitsToken_  ;
   edm::EDGetTokenT<std::vector<pat::Electron>>     elecsToken_        ;
   edm::EDGetTokenT<std::vector<pat::Muon>>         muonsToken_        ;
   edm::EDGetTokenT<std::vector<pat::Jet>>          jetsToken_         ;
   edm::EDGetTokenT<std::vector<pat::MET>>          metToken_          ;
   edm::EDGetTokenT<std::vector<pat::Tau>>          tausToken_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>      verticesToken_     ;
   edm::EDGetTokenT<reco::BeamSpot>                 bsToken_           ;
   edm::EDGetTokenT<std::vector<reco::Conversion>>  convToken_         ;


   PFJetIDSelectionFunctor jetIDLoose_ ;
   PFJetIDSelectionFunctor jetIDTight_ ;
   const ME0Geometry*      ME0Geometry_;

   TTree* mytree;
   int Ngenjet;
   float Ptgenjet[kMaxGenJet], Etagenjet[kMaxGenJet], Phigenjet[kMaxGenJet], Massgenjet[kMaxGenJet];

   //int Nstable, DecayID;
   //int ne, nve, nm, nvm, nt, nvt;
   int Nevt;
   float genWeight[kMaxWeights];

   int Nvtx;
   float Vtx_pt2[kMaxVertices];

   int Ngenlepton,Ngenphoton, Ngenparticle;
   int Statusgenphoton[kMaxParticle],PIdgenlepton[kMaxParticle],Chargegenlepton[kMaxParticle],Statusgenlepton[kMaxParticle];
   float Pgenphoton[kMaxParticle],Pxgenphoton[kMaxParticle],Pygenphoton[kMaxParticle],Pzgenphoton[kMaxParticle],Egenphoton[kMaxParticle],Ptgenphoton[kMaxParticle],Etagenphoton[kMaxParticle],Phigenphoton[kMaxParticle];
   float Pgenlepton[kMaxParticle],Pxgenlepton[kMaxParticle],Pygenlepton[kMaxParticle],Pzgenlepton[kMaxParticle],Egenlepton[kMaxParticle],Ptgenlepton[kMaxParticle],Etagenlepton[kMaxParticle],Phigenlepton[kMaxParticle],Massgenlepton[kMaxParticle],IsolationVargenlepton[kMaxParticle];


   float Ptgenparticle[kMaxParticle],Etagenparticle[kMaxParticle],Phigenparticle[kMaxParticle],Massgenparticle[kMaxParticle];
   int PIdgenparticle[kMaxParticle], Statusgenparticle[kMaxParticle], M1genparticle[kMaxParticle], M2genparticle[kMaxParticle], D1genparticle[kMaxParticle], D2genparticle[kMaxParticle] ;


   int Npho;
   float MVApho[kMaxPhoton],Ptpho[kMaxPhoton], Etapho[kMaxPhoton], Phipho[kMaxPhoton],Epho[kMaxPhoton],isLpho[kMaxPhoton],isTpho[kMaxPhoton];
   //TH2F * photonEcorr_{nullptr};

   int Nelec, Chargeelec[kMaxElectron], isLE[kMaxElectron], isME[kMaxElectron],isTE[kMaxElectron];
   float Ptelec[kMaxElectron], Etaelec[kMaxElectron], Phielec[kMaxElectron], IsolationVarelec[kMaxElectron], Pxelec[kMaxElectron], Pyelec[kMaxElectron], Pzelec[kMaxElectron], Eelec[kMaxElectron], MVAelec[kMaxElectron], Masselec[kMaxElectron];
   float OldIsolationVarelec[kMaxElectron];
   int MissingHitselec[kMaxElectron];
   int ChargeConsistentelec[kMaxElectron], ChargeConsistent1elec[kMaxElectron], ChargeConsistent2elec[kMaxElectron], ConvVetoElec[kMaxElectron];
   int isEBelec[kMaxElectron];
   float kfhitselec[kMaxElectron];
   float gsfhitselec[kMaxElectron];
   float kfchi2elec[kMaxElectron];
   float gsfchi2elec[kMaxElectron];
   float fbremelec[kMaxElectron];
   float eelepoutelec[kMaxElectron], deltaetaelec[kMaxElectron], deltaphielec[kMaxElectron];
   float oldsigmaietaietaelec[kMaxElectron];
   float oldsigmaiphiiphielec[kMaxElectron];
   float oldcircularityelec[kMaxElectron];
   float oldr9elec[kMaxElectron];
   float scletawidthelec[kMaxElectron];
   float sclphiwidthelec[kMaxElectron];
   float heelec[kMaxElectron];
   float epelec[kMaxElectron];
   float eseedpoutelec[kMaxElectron];
   float deltaetaseedelec[kMaxElectron];
   float deltaphiseedelec[kMaxElectron];




   int Nmuon, Chargemuon[kMaxMuonLoose];
   float Ptmuon[kMaxMuonLoose], Etamuon[kMaxMuonLoose], Phimuon[kMaxMuonLoose], IsolationVarmuon[kMaxMuonLoose],Pxmuon[kMaxMuonLoose], Pymuon[kMaxMuonLoose], Pzmuon[kMaxMuonLoose], Emuon[kMaxMuonLoose], Massmuon[kMaxMuonLoose];
   int isLM[kMaxMuonLoose], isTM[kMaxMuonLoose];
   bool CutBasedIdLoosemuon[kMaxMuonLoose], CutBasedIdMediummuon[kMaxMuonLoose], CutBasedIdTightmuon[kMaxMuonLoose];
   bool PFIsoLoosemuon[kMaxMuonLoose], PFIsoMediummuon[kMaxMuonLoose], PFIsoTightmuon[kMaxMuonLoose];

   int Njet, NoDjet[kMaxJet]; 
   int Lj[kMaxJet], Tj[kMaxJet] ;
   float Ptjet[kMaxJet], Etajet[kMaxJet], Phijet[kMaxJet], Massjet[kMaxJet];
   float NHEFjet[kMaxJet], NEEFjet[kMaxJet];
   float CHEFjet[kMaxJet], CEMEFjet[kMaxJet];
   float JECFjet[kMaxJet], bMVAjet[kMaxJet];
   float DeepCSVjet[kMaxJet];

   int Nmet;
   float Met[kMaxMissingET],Phimet[kMaxMissingET];


   int Ntau, Chargetau[kMaxTau];
   float DMtau[kMaxTau], ChargedIsotau[kMaxTau] ;
   float Pttau[kMaxTau], Etatau[kMaxTau], Phitau[kMaxTau], Masstau[kMaxTau];



 };

 //
 // constants, enums and typedefs
 //

 //
 // static data member definitions
 //

 //
 // constructors and destructor
 //
 Ntuplzr::Ntuplzr(const edm::ParameterSet& iConfig):
   //pileup_(iConfig.getParameter<unsigned int>("pileup")),
   debug_(iConfig.getParameter<bool>("debug")),
   genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
   genPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
   photnsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
   ecalRecHitsToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalRecHits"))),
   elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
   muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
   jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
   metToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("met"))),
   tausToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
   verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
   convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
   jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
   jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT)
 {
    //now do what ever initialization is needed


   if(debug_)  std::cout<<"Here I am : in constructor "<<std::endl;

   //if ( iConfig.existsAs<edm::FileInPath>("photonEcorr") ) {
   //  auto photonEcorrFile = iConfig.getParameter<edm::FileInPath>("photonEcorr");
   //  TFile * photonEcorr = TFile::Open(photonEcorrFile.fullPath().c_str());
   //  photonEcorr_ = (TH2F*) photonEcorr->Get("combinedECorrection");
   //  photonEcorr_->SetDirectory(0);  // don't delete
   //  delete photonEcorr;
   //}

   ME0Geometry_=0;
   usesResource("TFileService");
   mytree   = fs_->make<TTree>("mytree","TestTree");
   // mytree->Branch("Nstable",&Nstable, "Nstable/I");
   //mytree->Branch("DecayID",&DecayID, "DecayID/I");
   mytree->Branch("Nvtx",&Nvtx, "Nvtx/I");
   mytree->Branch("Vtx_pt2",Vtx_pt2, "Vtx_pt2[Nvtx]/F");

   mytree->Branch("Nevt",&Nevt, "Nevt/I");
   mytree->Branch("genWeight",genWeight, "genWeight[Nevt]/F");


   mytree->Branch("Ngenparticle",&Ngenparticle, "Ngenparticle/I");
   mytree->Branch("PIdgenparticle",PIdgenparticle, "PIdgenparticle[Ngenparticle]/I");
   mytree->Branch("Statusgenparticle",Statusgenparticle, "Statusgenparticle[Ngenparticle]/I");
   mytree->Branch("Ptgenparticle",Ptgenparticle, "Ptgenparticle[Ngenparticle]/F");
   mytree->Branch("Etagenparticle",Etagenparticle, "Etagenparticle[Ngenparticle]/F");
   mytree->Branch("Phigenparticle",Phigenparticle, "Phigenparticle[Ngenparticle]/F");
   mytree->Branch("Massgenparticle",Massgenparticle, "Massgenparticle[Ngenparticle]/F");
   mytree->Branch("M1genparticle",M1genparticle, "M1genparticle[Ngenparticle]/I");
   mytree->Branch("M2genparticle",M2genparticle, "M2genparticle[Ngenparticle]/I");
   mytree->Branch("D1genparticle",D1genparticle, "D1genparticle[Ngenparticle]/I");
   mytree->Branch("D2genparticle",D2genparticle, "D2genparticle[Ngenparticle]/I");


   mytree->Branch("Ngenlepton",&Ngenlepton, "Ngenlepton/I");
   mytree->Branch("PIdgenlepton",PIdgenlepton, "PIdgenlepton[Ngenlepton]/I");
   mytree->Branch("Chargegenlepton",Chargegenlepton, "Chargegenlepton[Ngenlepton]/I");
   mytree->Branch("Statusgenlepton",Statusgenlepton, "Statusgenlepton[Ngenlepton]/I");
   mytree->Branch("Pgenlepton",Pgenlepton, "Pgenlepton[Ngenlepton]/F");
   mytree->Branch("Pxgenlepton",Pxgenlepton, "Pxgenlepton[Ngenlepton]/F");
   mytree->Branch("Pygenlepton",Pygenlepton, "Pygenlepton[Ngenlepton]/F");
   mytree->Branch("Pzgenlepton",Pzgenlepton, "Pzgenlepton[Ngenlepton]/F");
   mytree->Branch("Egenlepton",Egenlepton, "Egenlepton[Ngenlepton]/F");
   mytree->Branch("Ptgenlepton",Ptgenlepton, "Ptgenlepton[Ngenlepton]/F");
   mytree->Branch("Etagenlepton",Etagenlepton, "Etagenlepton[Ngenlepton]/F");
   mytree->Branch("Phigenlepton",Phigenlepton, "Phigenlepton[Ngenlepton]/F");
   mytree->Branch("Massgenlepton",Massgenlepton, "Massgenlepton[Ngenlepton]/F");
   mytree->Branch("IsolationVargenlepton",IsolationVargenlepton, "IsolationVargenlepton[Ngenlepton]/F");




   mytree->Branch("Ngenphoton",&Ngenphoton, "Ngenphoton/I");
   mytree->Branch("Statusgenphoton",Statusgenphoton, "Statusgenphoton[Ngenphoton]/I");
   mytree->Branch("Pgenphoton",Pgenphoton, "Pgenphoton[Ngenphoton]/F");
   mytree->Branch("Pxgenphoton",Pxgenphoton, "Pxgenphoton[Ngenphoton]/F");
   mytree->Branch("Pygenphoton",Pygenphoton, "Pygenphoton[Ngenphoton]/F");
   mytree->Branch("Pzgenphoton",Pzgenphoton, "Pzgenphoton[Ngenphoton]/F");
   mytree->Branch("Egenphoton",Egenphoton, "Egenphoton[Ngenphoton]/F");
   mytree->Branch("Ptgenphoton",Ptgenphoton, "Ptgenphoton[Ngenphoton]/F");
   mytree->Branch("Etagenphoton",Etagenphoton, "Etagenphoton[Ngenphoton]/F");
   mytree->Branch("Phigenphoton",Phigenphoton, "Phigenphoton[Ngenphoton]/F");

   mytree->Branch("Ngenjet",&Ngenjet, "Ngenjet/I");
   mytree->Branch("Ptgenjet",Ptgenjet, "Ptgenjet[Ngenjet]/F");
   mytree->Branch("Etagenjet",Etagenjet, "Etagenjet[Ngenjet]/F");
   mytree->Branch("Phigenjet",Phigenjet, "Phigenjet[Ngenjet]/F");
   mytree->Branch("Massgenjet",Massgenjet, "Massgenjet[Ngenjet]/F");

   mytree->Branch("Npho",&Npho, "Npho/I");
   mytree->Branch("Ptpho",Ptpho, "Ptpho[Npho]/F");
   mytree->Branch("Etapho",Etapho, "Etapho[Npho]/F");
   mytree->Branch("Phipho",Phipho, "Phipho[Npho]/F");
   mytree->Branch("Epho",Epho, "Epho[Npho]/F");
   mytree->Branch("MVApho",MVApho, "MVApho[Npho]/F");
   mytree->Branch("isLpho",isLpho, "isLpho[Npho]/I");
   mytree->Branch("isTpho",isTpho, "isTpho[Npho]/I");

   mytree->Branch("Nelec",&Nelec, "Nelec/I");
   mytree->Branch("Ptelec",Ptelec, "Ptelec[Nelec]/F");
   mytree->Branch("Etaelec",Etaelec, "Etaelec[Nelec]/F");
   mytree->Branch("Phielec",Phielec, "Phielec[Nelec]/F");
   mytree->Branch("Chargeelec",Chargeelec, "Chargeelec[Nelec]/I");
   mytree->Branch("IsolationVarelec",IsolationVarelec, "IsolationVarelec[Nelec]/F");
   mytree->Branch("OldIsolationVarelec",OldIsolationVarelec, "OldIsolationVarelec[Nelec]/F");
   mytree->Branch("MissingHitselec",MissingHitselec, "MissingHitselec[Nelec]/I");
   mytree->Branch("ConvVetoElec",ConvVetoElec, "ConvVetoElec[Nelec]/I");
   mytree->Branch("ChargeConsistentelec",ChargeConsistentelec, "ChargeConsistentelec[Nelec]/I");
   mytree->Branch("ChargeConsistent1elec",ChargeConsistent1elec, "ChargeConsistent1elec[Nelec]/I");
   mytree->Branch("ChargeConsistent2elec",ChargeConsistent2elec, "ChargeConsistent2elec[Nelec]/I");
   mytree->Branch("Pxelec",Pxelec, "Pxelec[Nelec]/F");
   mytree->Branch("Pyelec",Pyelec, "Pyelec[Nelec]/F");
   mytree->Branch("Pzelec",Pzelec, "Pzelec[Nelec]/F");
   mytree->Branch("Eelec",Eelec, "Eelec[Nelec]/F");
   mytree->Branch("Masselec",Masselec, "Masselec[Nelec]/F");
   mytree->Branch("isLE", isLE, "isLE[Nelec]/I");
   mytree->Branch("isME", isME, "isME[Nelec]/I");
   mytree->Branch("isTE", isTE, "isTE[Nelec]/I");
   mytree->Branch("MVAelec",MVAelec, "MVAelec[Nelec]/F");
   mytree->Branch("isEBelec",isEBelec,"isEBelec[Nelec]/I");
   mytree->Branch("fbremelec",fbremelec,"fbremelec[Nelec]/F");
   mytree->Branch("kfhitselec",kfhitselec,"kfhitselec[Nelec]/F");
   mytree->Branch("kfchi2elec",kfchi2elec,"kfchi2elec[Nelec]/F");
   mytree->Branch("gsfhitselec",gsfhitselec,"gsfhitselec[Nelec]/F");
   mytree->Branch("gsfchi2elec",gsfchi2elec,"gsfchi2elec[Nelec]/F");
   mytree->Branch("eelepoutelec",eelepoutelec,"eelepoutelec[Nelec]/F");
   mytree->Branch("deltaetaelec",deltaetaelec,"deltaetaelec[Nelec]/F");
   mytree->Branch("deltaphielec",deltaphielec,"deltaphielec[Nelec]/F");
   mytree->Branch("oldsigmaietaietaelec",oldsigmaietaietaelec,"oldsigmaietaietaelec[Nelec]/F");
   mytree->Branch("oldsigmaiphiiphielec",oldsigmaiphiiphielec,"oldsigmaiphiiphielec[Nelec]/F");
   mytree->Branch("oldcircularityelec",oldcircularityelec,"oldcircularityelec[Nelec]/F");
   mytree->Branch("oldr9elec",oldr9elec,"oldr9elec[Nelec]/F");
   mytree->Branch("scletawidthelec",scletawidthelec,"scletawidthelec[Nelec]/F");
   mytree->Branch("sclphiwidthelec",sclphiwidthelec,"sclphiwidthelec[Nelec]/F");
   mytree->Branch("heelec",heelec,"heelec[Nelec]/F");
   mytree->Branch("epelec",epelec,"epelec[Nelec]/F");
   mytree->Branch("eseedpoutelec",eseedpoutelec,"eseedpoutelec[Nelec]/F");
   mytree->Branch("deltaetaseedelec",deltaetaseedelec,"deltaetaseedelec[Nelec]/F");
   mytree->Branch("deltaphiseedelec",deltaphiseedelec,"deltaphiseedelec[Nelec]/F");

   mytree->Branch("Nmuon",&Nmuon, "Nmuon/I");
   mytree->Branch("Ptmuon",Ptmuon, "Ptmuon[Nmuon]/F");
   mytree->Branch("Etamuon",Etamuon, "Etamuon[Nmuon]/F");
   mytree->Branch("Phimuon",Phimuon, "Phimuon[Nmuon]/F");
   mytree->Branch("Chargemuon",Chargemuon, "Chargemuon[Nmuon]/I");
   mytree->Branch("IsolationVarmuon",IsolationVarmuon, "IsolationVarmuon[Nmuon]/F");
   mytree->Branch("CutBasedIdLoosemuon",CutBasedIdLoosemuon, "CutBasedIdLoosemuon[Nmuon]/F");
   mytree->Branch("CutBasedIdMediummuon",CutBasedIdMediummuon, "CutBasedIdMediummuon[Nmuon]/F");
   // mytree->Branch("CutBasedIdMediumPromptmuon",CutBasedIdMediumPromptmuon, "CutBasedIdMediumPromptmuon[Nmuon]/F");
   mytree->Branch("CutBasedIdTightmuon",CutBasedIdTightmuon, "CutBasedIdTightmuon[Nmuon]/F");
   //mytree->Branch("CutBasedIdGlobalHighPtmuon",CutBasedIdGlobalHighPtmuon, "CutBasedIdGlobalHighPtmuon[Nmuon]/F");
   //mytree->Branch("CutBasedIdTrkHighPtmuon",CutBasedIdTrkHighPtmuon, "CutBasedIdTrkHighPtmuon[Nmuon]/F");
   mytree->Branch("PFIsoLoosemuon",PFIsoLoosemuon, "PFIsoLoosemuon[Nmuon]/F");
   mytree->Branch("PFIsoMediummuon",PFIsoMediummuon, "PFIsoMediummuon[Nmuon]/F");
   mytree->Branch("PFIsoTightmuon",PFIsoTightmuon, "PFIsoTightmuon[Nmuon]/F");
   //mytree->Branch("PFIsoTightmuon",PFIsoTightmuon, "PFIsoTightmuon[Nmuon]/F");
   mytree->Branch("Emuon",Emuon, "Emuon[Nmuon]/F");
   mytree->Branch("Massmuon",Massmuon, "Massmuon[Nmuon]/F");
   mytree->Branch("isLM", isLM , "isLM[Nmuon]/I");
   mytree->Branch("isTM", isTM, "isTM[Nmuon]/I");

   mytree->Branch("Njet",&Njet, "Njet/I");
   mytree->Branch("Ptjet",Ptjet, "Ptjet[Njet]/F");
   mytree->Branch("Etajet",Etajet, "Etajet[Njet]/F");
   mytree->Branch("Phijet",Phijet, "Phijet[Njet]/F");
   mytree->Branch("Lj", Lj, "Lj[Njet]/I");
   mytree->Branch("Tj", Tj, "Tj[Njet]/I");
   mytree->Branch("bMVAjet",bMVAjet, "bMVAjet[Njet]/F");
   mytree->Branch("DeepCSVjet",DeepCSVjet,"DeepCSVjet[Njet]/F");
   mytree->Branch("Massjet",Massjet, "Massjet[Njet]/F");
   mytree->Branch("NoDjet", NoDjet, "NoDjet[Njet]/I" );
   mytree->Branch("NHEFjet", NHEFjet, "NHEFjet[Njet]/F" );
   mytree->Branch("NEEFjet", NEEFjet, "NEEFjet[Njet]/F" );
   mytree->Branch("CHEFjet", CHEFjet, "CHEFjet[Njet]/F" );
   mytree->Branch("CEMEFjet", CEMEFjet, "CEMEFjet[Njet]/F" );
   mytree->Branch("JECFjet", JECFjet, "JECFjet[Njet]/F" );

   mytree->Branch("Nmet",&Nmet, "Nmet/I");
   mytree->Branch("Met", Met, "Met[Nmet]/F");
   mytree->Branch("Phimet",Phimet, "Phimet[Nmet]/F");



  mytree->Branch("Ntau",&Ntau, "Ntau/I");
  mytree->Branch("Pttau",Pttau, "Pttau[Ntau]/F");
  mytree->Branch("Etatau",Etatau, "Etatau[Ntau]/F");
  mytree->Branch("Phitau",Phitau, "Phitau[Ntau]/F");
  mytree->Branch("Chargetau",Chargetau, "Chargetau[Ntau]/I");
  mytree->Branch("Masstau",Masstau, "Masstau[Ntau]/F");
  mytree->Branch("DMtau",DMtau, "DMtau[Ntau]/F");
  mytree->Branch("ChargedIsotau",ChargedIsotau, "ChargedIsotau[Ntau]/F");
  //mytree->Branch("isLM", isLM , "isLM[Ntau]/I");
  //mytree->Branch("isTM", isTM, "isTM[Ntau]/I");
  //

  if(debug_)      std::cout<<"Here I am : ending constructor "<<std::endl;


}


Ntuplzr::~Ntuplzr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntuplzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   if(debug_) std::cout<<"here starts the event:"<<std::endl;
   Handle<std::vector<reco::GenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   
   Handle<std::vector<reco::GenJet>> genJets;
   iEvent.getByToken(genJetsToken_, genJets);

   Handle<std::vector<pat::Photon>> photns;
   iEvent.getByToken(photnsToken_, photns);
   
   Handle<EcalRecHitCollection> ecalRecHits;
   iEvent.getByToken(ecalRecHitsToken_, ecalRecHits);

   Handle<std::vector<pat::Electron>> elecs;
   iEvent.getByToken(elecsToken_, elecs);
   
   Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonsToken_, muons);
   
   Handle<std::vector<pat::Jet>> jets;
   iEvent.getByToken(jetsToken_, jets);
   
   Handle<std::vector<pat::MET>> met;
   iEvent.getByToken(metToken_, met);


   Handle<std::vector<pat::Tau>> taus;
   iEvent.getByToken(tausToken_, taus);

   Handle<std::vector<reco::Vertex>> vertices;
   iEvent.getByToken(verticesToken_, vertices);

   Handle<reco::ConversionCollection> conversions;
   iEvent.getByToken(convToken_, conversions);
   //const reco::Conversion &conversion = *conversions.product();  

   Handle<reco::BeamSpot> bsHandle;
   iEvent.getByToken(bsToken_, bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();  

   edm::Handle<GenEventInfoProduct> evt;
   iEvent.getByToken( generatorToken_,evt);

   edm::Handle<LHEEventProduct> evet;
   iEvent.getByToken(generatorlheToken_, evet);

   if(debug_) std::cout<<"Here I am : got handles right "<<std::endl;  
   //Nstable = 0;
   // ne      = 0;
   // nve     = 0;
   // nm      = 0;
   // nvm     = 0;
   // nt      = 0;
   // nvt     = 0;
   //DecayID = 0;
   Nevt       = 0;
   Nvtx       = 0;
   Ngenjet    = 0;
   Ngenlepton = 0;
   Ngenphoton = 0;
   Nelec      = 0;
   Njet       = 0;
   Nmuon      = 0;
   Nmet       = 0;
   Npho       = 0;
   Ntau       = 0;

   int prVtx = -1;
   for (size_t i = 0; i < vertices->size(); i++) {
     if (vertices->at(i).isFake()) continue;
     if (vertices->at(i).ndof() <= 4) continue;
     if (prVtx < 0) prVtx = i;
     Vtx_pt2[Nvtx] = vertices->at(i).p4().pt();
     Nvtx++;
   }
   if (prVtx < 0) return;
   if(debug_)  std::cout<<"Here I am : got vertex infor right "<<std::endl;

   // Generator weights 
   genWeight[0] = 1.0;
   if(evt.isValid()) {
     genWeight[0] = evt->weight();
     Nevt++;
     for (unsigned int i = 1; i < evt->weights().size(); i++) {
       genWeight[Nevt]=evt->weights()[i];
       Nevt++;
       if (Nevt > kMaxWeights) break;
     }
   }


   if(evet.isValid()) {
    double asdd=evet->originalXWGTUP();
    double sum=0.;
    for(unsigned int i=0  ; i<evet->weights().size();i++) {
      double asdde=evet->weights()[i].wgt;
      genWeight[Nevt]=genWeight[0]*asdde/asdd;
      sum+=genWeight[Nevt];
      Nevt++;
    }
   }  

   //GenJet information
   std::vector<size_t> jGenJets;                                            
   for(size_t ij= 0 ; ij < genJets->size(); ij++)
     {
       if (genJets->at(ij).pt() < 20.) continue;
       if (fabs(genJets->at(ij).eta()) > 5) continue;

       bool overlaps = false;
       for (size_t j = 0; j < genParts->size(); j++) {
	 if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
	 if (fabs(genJets->at(ij).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(ij).p4()) < 0.01) {
	   overlaps = true;
	   break;
	 }
       }
       if (overlaps) continue;
       jGenJets.push_back(ij);

       Ptgenjet[Ngenjet] = genJets->at(ij).pt();
       Etagenjet[Ngenjet]= genJets->at(ij).eta();
       Phigenjet[Ngenjet]= genJets->at(ij).phi();
       Massgenjet[Ngenjet]= genJets->at(ij).mass();
       Ngenjet++;

       if(Ngenjet>kMaxGenJet) break;       
     }
   if(debug_)   std::cout<<"Here I am : got genjet infor right "<<std::endl;



   //vector<reco::GenParticle>::const_iterator itParticle;
   vector<const reco::Candidate *> vectorCandidate;
   vector<const reco::Candidate *>::iterator itCandidate;

   for( vector<reco::GenParticle>::const_iterator itPackedParticle = genParts->begin(); itPackedParticle != genParts->end(); ++itPackedParticle)
     {
       vectorCandidate.push_back(&*itPackedParticle);
     }

   //GenParticle information
   for (size_t i = 0; i < genParts->size(); i++) {

     PIdgenparticle[Ngenparticle]       = genParts->at(i).pdgId();
     Statusgenparticle[Ngenparticle]    = genParts->at(i).status();
     Ptgenparticle[Ngenparticle]        = genParts->at(i).pt();
     Phigenparticle[Ngenparticle]       = genParts->at(i).phi();
     Etagenparticle[Ngenparticle]       = genParts->at(i).eta();
     Massgenparticle[Ngenparticle]      = genParts->at(i).mass();
     if(genParts->at(i).mother())
       {
	 std::cout<<"directly:"<<genParts->at(i).mother()->pdgId()<<std::endl;
     	 itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), genParts->at(i).mother());
	 if(itCandidate != vectorCandidate.end())
	   {
	     std::cout<<"index:"<< distance(vectorCandidate.begin(), itCandidate)<<std::endl;
	     M1genparticle[Ngenparticle] = distance(vectorCandidate.begin(), itCandidate);
	     M2genparticle[Ngenparticle] = distance(vectorCandidate.begin(), itCandidate);
	   }

       }
     else
       {
	 M1genparticle[Ngenparticle] =-99;
	 M2genparticle[Ngenparticle] =-99;
       }


     itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), genParts->at(i).daughter(0));
     if(itCandidate != vectorCandidate.end()) D1genparticle[Ngenparticle] = distance(vectorCandidate.begin(), itCandidate);

     itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), genParts->at(i).daughter(genParts->at(i).numberOfDaughters() - 1));
     if(itCandidate != vectorCandidate.end()) D2genparticle[Ngenparticle] = distance(vectorCandidate.begin(), itCandidate);
//     M1genparticle[Ngenparticle]   = genParts->at(i).mother1();
//     M2genparticle[Ngenparticle]   = genParts->at(i).mother2();
//     D1genparticle[Ngenparticle]   = genParts->at(i).daughter1();
//     D2genparticle[Ngenparticle]   = genParts->at(i).daughter2();
//     
     
     Ngenparticle++;
     if(Ngenparticle>kMaxParticle) break;       
   }
   if(debug_)   std::cout<<"Here I am : got genlep infor right "<<std::endl;



   //GenLep information
   for (size_t i = 0; i < genParts->size(); i++) {
     if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
     if (genParts->at(i).pt() < 10.) continue;
     if (fabs(genParts->at(i).eta()) > 3.) continue;
     double genIso = 0.;
   
     for (size_t j = 0; j < jGenJets.size(); j++) {
       if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
       std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
       for (size_t k = 0; k < jconst.size(); k++) {
	 double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
	 if (deltaR < 0.01 || deltaR > 0.4) continue;
	 genIso = genIso + jconst[k]->pt();
       }
     }   
     
     genIso = genIso / genParts->at(i).pt();
     PIdgenlepton[Ngenlepton]       = genParts->at(i).pdgId();
     Chargegenlepton[Ngenlepton]    = genParts->at(i).charge();
     Statusgenlepton[Ngenlepton]    = genParts->at(i).status();
     Pgenlepton[Ngenlepton]         = genParts->at(i).p();
     Pxgenlepton[Ngenlepton]        = genParts->at(i).px();
     Pygenlepton[Ngenlepton]        = genParts->at(i).py();
     Pzgenlepton[Ngenlepton]        = genParts->at(i).pz();
     Egenlepton[Ngenlepton]         = genParts->at(i).energy();
     Ptgenlepton[Ngenlepton]        = genParts->at(i).pt();
     Phigenlepton[Ngenlepton]       = genParts->at(i).phi();
     Etagenlepton[Ngenlepton]       = genParts->at(i).eta();
     Massgenlepton[Ngenlepton]      = genParts->at(i).mass();
     IsolationVargenlepton[Ngenlepton] = genIso; 
     Ngenlepton++;
     if(Ngenlepton>kMaxParticle) break;       
   }
   if(debug_)   std::cout<<"Here I am : got genlep infor right "<<std::endl;

   //GenPhoton information
   for (size_t i = 0; i < genParts->size(); i++) {
     if (abs(genParts->at(i).pdgId()) != 22) continue;
     if (genParts->at(i).pt() < 10.) continue;
     if (fabs(genParts->at(i).eta()) > 3.) continue;

     Statusgenphoton[Ngenphoton]= genParts->at(i).status();
     Pgenphoton[Ngenphoton]     = genParts->at(i).p();
     Pxgenphoton[Ngenphoton]    = genParts->at(i).px();
     Pygenphoton[Ngenphoton]    = genParts->at(i).py();
     Pzgenphoton[Ngenphoton]    = genParts->at(i).pz();
     Egenphoton[Ngenphoton]     = genParts->at(i).energy();
     Ptgenphoton[Ngenphoton]    = genParts->at(i).pt();
     Phigenphoton[Ngenphoton]   = genParts->at(i).phi();
     Etagenphoton[Ngenphoton]   = genParts->at(i).eta();
     Ngenphoton++;
     if(Ngenphoton>kMaxParticle) break;       
   }
   if(debug_)   std::cout<<"Here I am : got genpho infor right "<<std::endl;
   
   //Photon information                                                                                                         
   for(size_t ip= 0 ; ip < photns->size(); ip++)
     {
       if(photns->at(ip).pt() < 10.) continue;
       if(fabs(photns->at(ip).eta()) > 3.) continue;
       float mvaValue = photns->at(ip).userFloat("mvaValue");
       bool isEB = photns->at(ip).isEB();
       bool isLoose = 0;
       bool isTight = 0;
       if( isEB )
	 {
	   isLoose = (mvaValue > 0.00);
	   isTight = (mvaValue > 0.56);
	 }     
       else
	 {
	   isLoose = (mvaValue > 0.20);
	   isTight = (mvaValue > 0.68);
	 }          
       
       MVApho[Npho]= mvaValue;
       Ptpho[Npho] = photns->at(ip).pt();
       Etapho[Npho]= photns->at(ip).eta();
       Phipho[Npho]= photns->at(ip).phi();
       Epho[Npho]  = photns->at(ip).energy();
       isLpho[Npho]= (int) isLoose;
       isTpho[Npho]= (int) isTight;
       Npho++;
       if(Npho>kMaxPhoton) break;
     }
   if(debug_)   std::cout<<"Here I am : got pho infor right "<<std::endl;

   //Electron information                                                                                                       
   for(size_t ie= 0 ; ie < elecs->size(); ie++)  {
     if(elecs->at(ie).pt() < 10.) continue;
     if (fabs(elecs->at(ie).eta()) > 3.) continue;
     float mvaValue = elecs->at(ie).userFloat("mvaValue");
     //float mvaValue = 0.25 ;
     bool isEB = elecs->at(ie).isEB();
     
     
     bool isLoose = 0;
     bool isMedium = 0;
     bool isTight = 0;

     if( isEB ) {
       if (elecs->at(ie).pt() < 20.) {
	 isLoose = (mvaValue > -0.661);
	 isMedium = (mvaValue > 0.885);
	 isTight = (mvaValue > 0.986);
	 

       }
       else {
	 isLoose = (mvaValue > -0.797);
	 isMedium = (mvaValue > 0.723);
	 isTight = (mvaValue > 0.988);
       }
     }
     else {
       if (not (elecs->at(ie).userFloat("hgcElectronID:ecEnergy") > 0)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:sigmaUU") > 0)) continue;
       if (not (elecs->at(ie).fbrem() > -1)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:measuredDepth") < 40)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:nLayers") > 20)) continue;
       if (elecs->at(ie).pt() < 20.) {
	 isLoose = (mvaValue > -0.320);
	 isMedium = (mvaValue > 0.777);
	 isTight = (mvaValue > 0.969);
       }
       else {
	 isLoose = (mvaValue > -0.919);
	 isMedium = (mvaValue > 0.591);
	 isTight = (mvaValue > 0.983);
       }
     }
     //isLE[Nelec] = isLooseElec(elecs->at(ie),conversion,beamspot);
     //isTE[Nelec] = isTightElec(elecs->at(ie),conversion,beamspot);

     if(isEB) 
       IsolationVarelec[Nelec] = (elecs->at(ie).puppiNoLeptonsChargedHadronIso() + elecs->at(ie).puppiNoLeptonsNeutralHadronIso() + elecs->at(ie).puppiNoLeptonsPhotonIso()) / elecs->at(ie).pt();
     else 
       IsolationVarelec[Nelec] = (elecs->at(ie).userFloat("hgcElectronID:caloIsoRing1") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing2") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing3") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing4")) / elecs->at(ie).energy();


     isLE[Nelec]                 = (int) isLoose;
     isME[Nelec]                 = (int) isMedium;
     isTE[Nelec]                 = (int) isTight;
     MVAelec[Nelec]              = mvaValue;
     Ptelec[Nelec]               = elecs->at(ie).pt();
     Etaelec[Nelec]              = elecs->at(ie).eta();
     Phielec[Nelec]              = elecs->at(ie).phi();
     Chargeelec[Nelec]           = elecs->at(ie).charge();
     Pxelec[Nelec]               = elecs->at(ie).pz();
     Pyelec[Nelec]               = elecs->at(ie).py();
     Pzelec[Nelec]               = elecs->at(ie).pz();
     Masselec[Nelec]             = elecs->at(ie).mass();
     Eelec[Nelec]                = elecs->at(ie).energy();
     OldIsolationVarelec[Nelec]  = elecs->at(ie).dr03TkSumPt()/elecs->at(ie).pt(); 
     MissingHitselec[Nelec]      = elecs->at(ie).gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); 
     ChargeConsistentelec[Nelec] = (int) elecs->at(ie).isGsfCtfScPixChargeConsistent();
     ChargeConsistent1elec[Nelec] = (int) elecs->at(ie).isGsfScPixChargeConsistent();
     ChargeConsistent2elec[Nelec] = (int) elecs->at(ie).isGsfCtfChargeConsistent();
     ConvVetoElec[Nelec]         = (int) elecs->at(ie).passConversionVeto();
     isEBelec[Nelec]             = (int) isEB;
     fbremelec[Nelec]            = elecs->at(ie).fbrem();
     reco::TrackRef myTrackRef   = elecs->at(ie).closestCtfTrackRef();
     bool validKF                = myTrackRef.isAvailable() && myTrackRef.isNonnull();
     kfhitselec[Nelec]           = (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.;
     kfchi2elec[Nelec]           = (validKF) ? myTrackRef->normalizedChi2() : -1.;
     gsfhitselec[Nelec]          = elecs->at(ie).gsfTrack()->hitPattern().trackerLayersWithMeasurement();
     gsfchi2elec[Nelec]          = elecs->at(ie).gsfTrack()->normalizedChi2();
     eelepoutelec[Nelec]         = elecs->at(ie).eEleClusterOverPout();
     deltaetaelec[Nelec]         = elecs->at(ie).deltaEtaEleClusterTrackAtCalo();
     deltaphielec[Nelec]         = elecs->at(ie).deltaPhiEleClusterTrackAtCalo();
     oldsigmaietaietaelec[Nelec] = elecs->at(ie).full5x5_sigmaIetaIeta();
     oldsigmaiphiiphielec[Nelec] = elecs->at(ie).full5x5_sigmaIphiIphi();
     oldcircularityelec[Nelec]   = ( elecs->at(ie).showerShape().e5x5!=0.) ? 1.-elecs->at(ie).showerShape().e1x5/elecs->at(ie).showerShape().e5x5 : -1. ;
     oldr9elec[Nelec]            = elecs->at(ie).full5x5_r9();
     scletawidthelec[Nelec]      = elecs->at(ie).superCluster()->etaWidth();
     sclphiwidthelec[Nelec]      = elecs->at(ie).superCluster()->phiWidth();
     heelec[Nelec]               = elecs->at(ie).hcalOverEcal();
     epelec[Nelec]               = elecs->at(ie).eSuperClusterOverP();
     eseedpoutelec[Nelec]        = elecs->at(ie).eEleClusterOverPout();
     deltaetaseedelec[Nelec]     = elecs->at(ie).deltaEtaSeedClusterTrackAtCalo();
     deltaphiseedelec[Nelec]     = elecs->at(ie).deltaPhiSeedClusterTrackAtCalo();
     
     //std::cout<<"elec Mass:"<<elecs->at(ie).mass()<<std::endl;
     
     Nelec++;
    if(Nelec>kMaxElectron) break;
   }
   //std::cout<<Nelec<<std::endl;
   if(debug_)   std::cout<<"Here I am : got elec infor right "<<std::endl;
   // Muon information                                                                                                     
   for(size_t im= 0 ; im < muons->size(); im++)
     {
    if (muons->at(im).pt() < 2.) continue;
    if (fabs(muons->at(im).eta()) > 2.8) continue;

	   Ptmuon[Nmuon] = muons->at(im).pt();
	   Etamuon[Nmuon]= muons->at(im).eta();
	   Phimuon[Nmuon]= muons->at(im).phi();
	   Chargemuon[Nmuon]= muons->at(im).charge();
	   Pxmuon[Nmuon] = muons->at(im).px();
	   Pymuon[Nmuon] = muons->at(im).py();
	   Pzmuon[Nmuon] = muons->at(im).pz();
	   Emuon[Nmuon]= muons->at(im).energy();
	   Massmuon[Nmuon]= muons->at(im).mass();
	   // IsolationVarmuon[Nmuon] = (muons->at(im).puppiNoLeptonsChargedHadronIso() + muons->at(im).puppiNoLeptonsNeutralHadronIso() + muons->at(im).puppiNoLeptonsPhotonIso()) / muons->at(im).pt();
	   IsolationVarmuon[Nmuon] = muons->at(im).trackIso()/muons->at(im).pt();

	   CutBasedIdLoosemuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdLoose);
	   CutBasedIdMediummuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdMedium);
	   CutBasedIdTightmuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdTight);
	   PFIsoLoosemuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoLoose);
	   PFIsoMediummuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoMedium);

	   PFIsoTightmuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoTight);

	   double dPhiCut = std::min(std::max(1.2/muons->at(im).p(),1.2/100),0.056);
	   double dPhiBendCut = std::min(std::max(0.2/muons->at(im).p(),0.2/100),0.0096);
	   isLM[Nmuon] = (int) (fabs(muons->at(im).eta()) < 2.4 && muon::isLooseMuon(muons->at(im))) || (fabs(muons->at(im).eta()) > 2.4 && isME0MuonSelNew(muons->at(im), 0.077, dPhiCut, dPhiBendCut, iSetup));
	 	   
	   bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
	   if (muons->at(im).innerTrack().isNonnull()){
	     ipxy = std::abs(muons->at(im).muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
	     ipz = std::abs(muons->at(im).muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
	     validPxlHit = muons->at(im).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
	     highPurity = muons->at(im).innerTrack()->quality(reco::Track::highPurity);
	   }      
	   dPhiCut = std::min(std::max(1.2/muons->at(im).p(),1.2/100),0.032);
	   dPhiBendCut = std::min(std::max(0.2/muons->at(im).p(),0.2/100),0.0041);
	   isTM[Nmuon] = (int) (fabs(muons->at(im).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(im),vertices->at(prVtx))) || (fabs(muons->at(im).eta()) > 2.4 && isME0MuonSelNew(muons->at(im), 0.048, dPhiCut, dPhiBendCut, iSetup) && ipxy && ipz && validPxlHit && highPurity);
	   //isTM[Nmuon] = (int) muon::isTightMuon(muons->at(im),vertices->at(prVtx));
	   //std::cout<<"Muon Mass:"<<muons->at(im).mass()<<std::endl;
	   Nmuon++;
	   if(Nmuon>kMaxMuonLoose) break;
     }
     
   if(debug_)   std::cout<<"Here I am : got muon infor right "<<std::endl;
   //std::cout<<"Muon"<<std::endl;

   //Jet information                                                                                                            
   for(size_t ij= 0 ; ij < jets->size(); ij++)
     {
       //    cout<<"ij = "<<ij<<" "<<jets->size()<<" "<<jets->at(ij).pt()<<" "<<jets->at(ij).eta()<<endl;

       if (jets->at(ij).pt() < 20.) continue;
       if (fabs(jets->at(ij).eta()) > 5) continue;
       
       Ptjet[Njet]  = jets->at(ij).pt();
       Etajet[Njet] = jets->at(ij).eta();
       Phijet[Njet] = jets->at(ij).phi();
       Massjet[Njet]= jets->at(ij).mass();
       
       NoDjet[Njet]   = jets->at(ij).numberOfDaughters();
       NHEFjet[Njet] = jets->at(ij).neutralHadronEnergyFraction();
       NEEFjet[Njet] = jets->at(ij).neutralEmEnergyFraction();
       CHEFjet[Njet] = jets->at(ij).chargedHadronEnergyFraction();
       CEMEFjet[Njet] = jets->at(ij).chargedEmEnergyFraction();
       JECFjet[Njet] = jets->at(ij).pt()/jets->at(ij).correctedJet("Uncorrected").pt();
       if( (jets->at(ij).numberOfDaughters() > 1 ) && ( jets->at(ij).neutralEmEnergyFraction()< 0.99 ) && ( jets->at(ij).neutralHadronEnergyFraction() < 0.99 ))
	 {
	   Lj[Njet] = 1; 
	 }
       else Lj[Njet] = 0;
       //std::cout<<"Tight ID Njet : " << Njet <<std::endl;
       //pat::strbitset retTight = jetIDTight_.getBitTemplate();
       //retTight.set(false);
       //Tj[Njet]=(int) jetIDTight_(jets->at(ij), retTight);
       if( (jets->at(ij).numberOfDaughters() > 1 ) && ( jets->at(ij).neutralEmEnergyFraction()< 0.9 ) && ( jets->at(ij).neutralHadronEnergyFraction() < 0.9 ))
       	 {
	   //  std::cout<< "This jet should pass Tight ID " <<  std::endl;
       	   Tj[Njet]= 1;
       	 }
       else Tj[Njet]= 0;
       //std::cout<<"jet Mass:"<<jets->at(ij).mass()<<std::endl;

       bMVAjet[Njet]         = (float) jets->at(ij).bDiscriminator("pfCombinedMVAV2BJetTags");  ;
       DeepCSVjet[Njet]      = (float) jets->at(ij).bDiscriminator("pfDeepCSVJetTags:probb") + jets->at(ij).bDiscriminator("pfDeepCSVJetTags:probbb");
          
       Njet++;
       if(Njet>kMaxJet) break;
     }
   
   if(debug_)   std::cout<<"Here I am : got jet infor right "<<std::endl;
    //Met information      
    for(size_t imet= 0 ; imet < met->size(); imet++)
      {
        Met[Nmet]    = met->at(imet).pt();
        Phimet[Nmet] = met->at(imet).phi();
        Nmet++;
      }


   if(debug_)   std::cout<<"Here I am : got met infor right "<<std::endl;
    // Taus info

    for (size_t it = 0; it < taus->size(); it++) {
      if (taus->at(it).pt()<15.) continue; 
      if (fabs(taus->at(it).eta()) > 3.0) continue;
      if (taus->at(it).tauID("decayModeFinding")<0) continue;     
      
      Chargetau[Ntau]     = taus->at(it).charge();
      Pttau[Ntau]         = taus->at(it).pt();
      Etatau[Ntau]        = taus->at(it).eta();
      Phitau[Ntau]        = taus->at(it).phi();
      Masstau[Ntau]       = taus->at(it).mass();
      DMtau[Ntau]         = taus->at(it).decayMode();
      ChargedIsotau[Ntau] = taus->at(it).tauID("chargedIsoPtSum");
      Ntau++;
      if(Ntau>kMaxTau) break; 
    }
    if(debug_)   std::cout<<"Here I am : got tau infor right "<<std::endl;
    if(debug_) std::cout<<"beforeend"<<std::endl;
 
   mytree->Fill();
   if(debug_) std::cout<<"end"<<std::endl;
}





// ------------ method check that an e passes loose ID ----------------------------------
/*  bool
  Ntuplzr::isLooseElec(const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
  if (patEl.hcalOverEcal() > 6.741) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 73.76) return false;
  //if (ConversionTools::hasMatchedConversion(patEl, convCol, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
  bool
  Ntuplzr::isMediumElec(const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
  if (patEl.hcalOverEcal() > 7.371) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 22.6) return false;
  //if (ConversionTools::hasMatchedConversion(patEl, convCol, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
  bool
  Ntuplzr::isTightElec(const pat::Electron & patEl, const reco::ConversionCollection &convCol, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
  if (patEl.hcalOverEcal() > 4.492) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 18.26) return false;
  //if (ConversionTools::hasMatchedConversion(patEl, convCol, beamspot.position())) return false;
  return true;
}
*/



// ------------ method to improve ME0 muon ID ----------------
  bool 
Ntuplzr::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = std::abs(chamber->x - segment->x);
          deltaY   = std::abs(chamber->y - segment->y);
          pullX    = std::abs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = std::abs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = std::abs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

bool 
Ntuplzr::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut, edm::EventSetup const& iSetup)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;


    if(!ME0Geometry_)
      //throw std::runtime_error("Ntuplzr::isME0MuonSelNew: muon geometry not loaded");
      {
	edm::ESHandle<ME0Geometry> hGeom;
	iSetup.get<MuonGeometryRecord>().get(hGeom);
	ME0Geometry_ =( &*hGeom);
	if(!ME0Geometry_) return result;
      }

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      if (chamber->detector() == 5){

        for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

          LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
          LocalPoint seg_loc_coord(segment->x, segment->y, 0);
          LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
          LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

          const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);
	  if(!me0chamber)continue;

          GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
          GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);

          //double segDPhi = segment->me0SegmentRef->deltaPhi();
          // need to check if this works
          double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
          double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);

          deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
          deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
          deltaPhiBend = std::abs(segDPhi - trackDPhi);

          if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;

        }
      }
    }

  }

  return result;

}






// ------------ method called once each job just before starting event loop  ------------
void 
Ntuplzr::beginJob()
{

}


void
Ntuplzr::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);


  /*  edm::Handle<LHERunInfoProduct> run; 
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
 
  iRun.getByToken( "externalLHEProducer", abc );
  LHERunInfoProduct myLHERunInfoProduct = *(abc.product());
 
  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
    }*/
}

// ------------ method called when ending the processing of a run  ------------
  void
  Ntuplzr::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }


// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuplzr::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplzr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplzr);

//  LocalWords:  im
