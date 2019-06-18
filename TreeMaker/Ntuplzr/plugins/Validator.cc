// -*- C++ -*-
//
// Package:    TreeMaker/Ntuplzr
// Class:      Validator
// 
/**\class Validator Validator.cc TreeMaker/Validator/plugins/Validator.cc

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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
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


class Validator : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit Validator(const edm::ParameterSet&);
  ~Validator();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) ;
  virtual void endJob() override;

  bool isME0MuonSelNew(reco::Muon, double, double, double, edm::EventSetup const& );
  float calculate_demetraIsolation(const pat::Tau&) const;
  bool debug_;
  edm::Service<TFileService> fs_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      verticesToken_     ;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_     ;
  edm::EDGetTokenT<std::vector<reco::GenJet>>      genJetsToken_      ;
  edm::EDGetTokenT<std::vector<pat::Photon>>       photnsToken_       ; 
  edm::EDGetTokenT<std::vector<pat::Electron>>     elecsToken_        ;
  edm::EDGetTokenT<std::vector<pat::Muon>>         muonsToken_        ;
  edm::EDGetTokenT<std::vector<pat::Tau>>          tausToken_         ;
  edm::EDGetTokenT<std::vector<pat::Jet>>          jetsToken_         ;
  edm::EDGetTokenT<std::vector<pat::MET>>          metToken_          ;
  //edm::EDGetTokenT<std::vector<reco::Conversion>>  convToken_         ;

  const ME0Geometry*      ME0Geometry_;

  TTree* mytree;
  int Nevt;
  
  int Nvtx;
  float Vtx_pt2[kMaxVertices];
  
  int Ngenparticle;
  float Ptgenparticle[kMaxParticle],Etagenparticle[kMaxParticle],Phigenparticle[kMaxParticle],Massgenparticle[kMaxParticle];
  int PIdgenparticle[kMaxParticle], Statusgenparticle[kMaxParticle], M1genparticle[kMaxParticle], M2genparticle[kMaxParticle], D1genparticle[kMaxParticle], D2genparticle[kMaxParticle] ;

  int Ngenjet;
  float Ptgenjet[kMaxGenJet], Etagenjet[kMaxGenJet], Phigenjet[kMaxGenJet], Massgenjet[kMaxGenJet];

  int Npho;
  float MVApho[kMaxPhoton],Ptpho[kMaxPhoton], Etapho[kMaxPhoton], Phipho[kMaxPhoton],Masspho[kMaxPhoton],Isolationpho[kMaxPhoton], IDVarpho[kMaxPhoton];
  uint32_t IsoPasspho[kMaxPhoton], IDPasspho[kMaxPhoton];

  int Nelec, Chargeelec[kMaxElectron];
  float Ptelec[kMaxElectron], Etaelec[kMaxElectron], Phielec[kMaxElectron], Isolationelec[kMaxElectron], MVAelec[kMaxElectron], Masselec[kMaxElectron], IDVarelec[kMaxElectron];
  uint32_t IsoPasselec[kMaxElectron], IDPasselec[kMaxElectron];


  int Nmuon, Chargemuon[kMaxMuonLoose];
  float Ptmuon[kMaxMuonLoose], Etamuon[kMaxMuonLoose], Phimuon[kMaxMuonLoose], Isolationmuon[kMaxMuonLoose], Massmuon[kMaxMuonLoose], IDVarmuon[kMaxMuonLoose];
  uint32_t IsoPassmuon[kMaxMuonLoose], IDPassmuon[kMaxMuonLoose] ;


  int Ntau, Chargetau[kMaxTau];
  float DMtau[kMaxTau], Isolationtau[kMaxTau], Pttau[kMaxTau], Etatau[kMaxTau], Phitau[kMaxTau], Masstau[kMaxTau], IsoFunctiontau[kMaxTau];
  uint32_t IsoPasstau[kMaxTau];

  int Njet; 
  float Ptjet[kMaxJet], Etajet[kMaxJet], Phijet[kMaxJet], Massjet[kMaxJet];
  uint32_t IDPassjet[kMaxJet];
  float bMVAjet[kMaxJet];
  float DeepCSVjetTag1[kMaxJet];
  float DeepCSVjetTag2[kMaxJet];

  int Nmet;
  float Met[kMaxMissingET],Phimet[kMaxMissingET];
};



Validator::Validator(const edm::ParameterSet& iConfig):
  debug_(iConfig.getParameter<bool>("debug")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  photnsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  tausToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("met")))

{

  if(debug_)  std::cout<<"Here I am : in constructor "<<std::endl;

  ME0Geometry_ = 0;
  Nevt         = 0;

  usesResource("TFileService");
  mytree   = fs_->make<TTree>("mytree","TestTree");
  mytree->Branch("Nevt",&Nevt, "Nevt/I");
  mytree->Branch("Nvtx",&Nvtx, "Nvtx/I");
  mytree->Branch("Vtx_pt2",Vtx_pt2, "Vtx_pt2[Nvtx]/F");

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

  mytree->Branch("Ngenjet",&Ngenjet, "Ngenjet/I");
  mytree->Branch("Ptgenjet",Ptgenjet, "Ptgenjet[Ngenjet]/F");
  mytree->Branch("Etagenjet",Etagenjet, "Etagenjet[Ngenjet]/F");
  mytree->Branch("Phigenjet",Phigenjet, "Phigenjet[Ngenjet]/F");
  mytree->Branch("Massgenjet",Massgenjet, "Massgenjet[Ngenjet]/F");

  mytree->Branch("Npho",&Npho, "Npho/I");
  mytree->Branch("Ptpho",Ptpho, "Ptpho[Npho]/F");
  mytree->Branch("Etapho",Etapho, "Etapho[Npho]/F");
  mytree->Branch("Phipho",Phipho, "Phipho[Npho]/F");
  mytree->Branch("Masspho",Masspho, "Masspho[Npho]/F");
  mytree->Branch("IDVarpho", IDVarpho, "IDVarpho[Npho]/F");
  mytree->Branch("Isolationpho",Isolationpho, "Isolationpho[Npho]/F");
  mytree->Branch("IDPasspho", IDPasspho, "IDPasspho[Npho]/i");
  mytree->Branch("IsoPasspho", IsoPasspho, "IsoPasspho[Npho]/i");


  mytree->Branch("Nelec",&Nelec, "Nelec/I");
  mytree->Branch("Ptelec",Ptelec, "Ptelec[Nelec]/F");
  mytree->Branch("Etaelec",Etaelec, "Etaelec[Nelec]/F");
  mytree->Branch("Phielec",Phielec, "Phielec[Nelec]/F");
  mytree->Branch("Masselec",Masselec, "Masselec[Nelec]/F");
  mytree->Branch("Chargeelec",Chargeelec, "Chargeelec[Nelec]/I");
  mytree->Branch("IDVarelec", IDVarelec, "IDVarelec[Nelec]/F");
  mytree->Branch("Isolationelec",Isolationelec, "Isolationelec[Nelec]/F");
  mytree->Branch("IDPasselec", IDPasselec, "IDPasselec[Nelec]/i");
  mytree->Branch("IsoPasselec", IsoPasselec, "IsoPasselec[Nelec]/i");
  
  mytree->Branch("Nmuon",&Nmuon, "Nmuon/I");
  mytree->Branch("Ptmuon",Ptmuon, "Ptmuon[Nmuon]/F");
  mytree->Branch("Etamuon",Etamuon, "Etamuon[Nmuon]/F");
  mytree->Branch("Phimuon",Phimuon, "Phimuon[Nmuon]/F");
  mytree->Branch("Massmuon",Massmuon, "Massmuon[Nmuon]/F");
  mytree->Branch("Chargemuon",Chargemuon, "Chargemuon[Nmuon]/I");
  mytree->Branch("IDVarmuon", IDVarmuon, "IDVarmuon[Nmuon]/F");
  mytree->Branch("Isolationmuon",Isolationmuon, "Isolationmuon[Nmuon]/F");
  mytree->Branch("IDPassmuon", IDPassmuon, "IDPassmuon[Nmuon]/i");
  mytree->Branch("IsoPassmuon", IsoPassmuon, "IsoPassmuon[Nmuon]/i");
  
  mytree->Branch("Ntau",&Ntau, "Ntau/I");
  mytree->Branch("Pttau",Pttau, "Pttau[Ntau]/F");
  mytree->Branch("Etatau",Etatau, "Etatau[Ntau]/F");
  mytree->Branch("Phitau",Phitau, "Phitau[Ntau]/F");
  mytree->Branch("Masstau",Masstau, "Masstau[Ntau]/F");
  mytree->Branch("Chargetau",Chargetau, "Chargetau[Ntau]/I");
  mytree->Branch("DMtau",DMtau, "DMtau[Ntau]/F");
  mytree->Branch("Isolationtau",Isolationtau, "Isolationtau[Ntau]/F");
  mytree->Branch("IsoFunctiontau",IsoFunctiontau, "IsoFunctiontau[Ntau]/F");
  mytree->Branch("IsoPasstau", IsoPasstau, "IsoPasstau[Ntau]/i");


  mytree->Branch("Njet",&Njet, "Njet/I");
  mytree->Branch("Ptjet",Ptjet, "Ptjet[Njet]/F");
  mytree->Branch("Etajet",Etajet, "Etajet[Njet]/F");
  mytree->Branch("Phijet",Phijet, "Phijet[Njet]/F");
  mytree->Branch("Massjet",Massjet, "Massjet[Njet]/F");
  mytree->Branch("IDPassjet", IDPassjet, "IDPassjet[Njet]/i");
  mytree->Branch("bMVAjet",bMVAjet, "bMVAjet[Njet]/F");
  mytree->Branch("DeepCSVjetTag1",DeepCSVjetTag1,"DeepCSVjetTag1[Njet]/F");
  mytree->Branch("DeepCSVjetTag2",DeepCSVjetTag2,"DeepCSVjetTag2[Njet]/F");
  

  mytree->Branch("Nmet",&Nmet, "Nmet/I");
  mytree->Branch("Met", Met, "Met[Nmet]/F");
  mytree->Branch("Phimet",Phimet, "Phimet[Nmet]/F");


 



  if(debug_)      std::cout<<"Here I am : ending constructor "<<std::endl;


}


Validator::~Validator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Validator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   if(debug_) std::cout<<"here starts the event:"<<std::endl;

   Handle<std::vector<reco::Vertex>> vertices;
   iEvent.getByToken(verticesToken_, vertices);

   Handle<std::vector<reco::GenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   
   Handle<std::vector<reco::GenJet>> genJets;
   iEvent.getByToken(genJetsToken_, genJets);

   Handle<std::vector<pat::Photon>> photns;
   iEvent.getByToken(photnsToken_, photns);
   
   Handle<std::vector<pat::Electron>> elecs;
   iEvent.getByToken(elecsToken_, elecs);
   
   Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonsToken_, muons);

   Handle<std::vector<pat::Tau>> taus;
   iEvent.getByToken(tausToken_, taus);
   
   Handle<std::vector<pat::Jet>> jets;
   iEvent.getByToken(jetsToken_, jets);
   
   Handle<std::vector<pat::MET>> met;
   iEvent.getByToken(metToken_, met);

   if(debug_) std::cout<<"Here I am : got handles right "<<std::endl;  
   Nvtx          = 0;
   Ngenjet       = 0;
   Ngenparticle  = 0;
   Nelec         = 0;
   Njet          = 0;
   Nmuon         = 0;
   Nmet          = 0;
   Npho          = 0;
   Ntau          = 0;

   if(debug_) std::cout<<"Here I am : initalised number of particles=0 in the event "<<std::endl;  

   Nevt++;


   /////////////////////////////
   //////vertices info//////
   /////////////////////////////

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




   
   /////////////////////////////
   //////Gen Particle info//////
   /////////////////////////////
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
     	 itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), genParts->at(i).mother());
	 if(itCandidate != vectorCandidate.end())
	   {
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
     Ngenparticle++;
     if(Ngenparticle>kMaxParticle) break;       
   }
   if(debug_)   std::cout<<"Here I am : got genpart infor right "<<std::endl;


   /////////////////////////////
   //GenJet information
   /////////////////////////////
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
   
   /////////////////////////////
   //Photon information         
   /////////////////////////////                                                                                                
   for(size_t ip= 0 ; ip < photns->size(); ip++)
     {
       if(photns->at(ip).pt() < 10.) continue;
       if(fabs(photns->at(ip).eta()) > 3.) continue;
       float mvaValue = photns->at(ip).userFloat("mvaValue");
       bool isEB = photns->at(ip).isEB();
       bool isLoose(0), isMedium(0), isTight(0);

       if( isEB )
       	 {
       	   isLoose  = (mvaValue > 0.00);
 	   isMedium = (mvaValue > 0.2); 
       	   isTight  = (mvaValue > 0.56);
       	 }     
       else
       	 {
       	   isLoose = (mvaValue > 0.20);
 	   isMedium = (mvaValue > 0.4); 
       	   isTight = (mvaValue > 0.68);
       	 }          
       
       Ptpho[Npho]        = photns->at(ip).pt();
       Etapho[Npho]       = photns->at(ip).eta();
       Phipho[Npho]       = photns->at(ip).phi();
       Masspho[Npho]      = photns->at(ip).mass();
       IDVarpho[Npho]     = mvaValue; // MVA
       Isolationpho[Npho] = 1.;
       

       if(isLoose)
	 IDPasspho[Npho] |= 1 << 0;

       if(isMedium)
	 IDPasspho[Npho] |= 1 << 1;
       
       if(isTight)
	 IDPasspho[Npho] |= 1 << 2;
       
       if(Isolationpho[Npho] < 0.1)
	 IsoPasspho[Npho] |= 1 << 0;

       if(Isolationpho[Npho] < 0.2)
	 IsoPasspho[Npho] |= 1 << 1;

       if(Isolationpho[Npho] < 0.3)
	 IsoPasspho[Npho] |= 1 << 2;

       if(Isolationpho[Npho] < 0.4)
	 IsoPasspho[Npho] |= 1 << 3;

       

       Npho++;
       if(Npho>kMaxPhoton) break;
     }
   if(debug_)   std::cout<<"Here I am : got pho infor right "<<std::endl;

   /////////////////////////////
   //Electron information
   /////////////////////////////                       
   for(size_t ie= 0 ; ie < elecs->size(); ie++)  {
     if(elecs->at(ie).pt() < 10.) continue;
     if (fabs(elecs->at(ie).eta()) > 3.) continue;
     float mvaValue = elecs->at(ie).userFloat("mvaValue");
     Ptelec[Nelec]               = elecs->at(ie).pt();
     Etaelec[Nelec]              = elecs->at(ie).eta();
     Phielec[Nelec]              = elecs->at(ie).phi();
     Masselec[Nelec]             = elecs->at(ie).mass();
     Chargeelec[Nelec]           = elecs->at(ie).charge();
     IDVarelec[Nelec]            = mvaValue; //MVA
     bool isEB = elecs->at(ie).isEB();
     if(isEB) 
       Isolationelec[Nelec] = (elecs->at(ie).puppiNoLeptonsChargedHadronIso() + elecs->at(ie).puppiNoLeptonsNeutralHadronIso() + elecs->at(ie).puppiNoLeptonsPhotonIso()) / elecs->at(ie).pt();
     else 
       Isolationelec[Nelec] = (elecs->at(ie).userFloat("hgcElectronID:caloIsoRing1") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing2") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing3") + elecs->at(ie).userFloat("hgcElectronID:caloIsoRing4")) / elecs->at(ie).energy();

     bool isLoose(0), isMedium(0), isTight(0);
   
     if( isEB ) {
       if (elecs->at(ie).pt() < 20.) {
     	 isLoose = (mvaValue  > -0.661);
     	 isMedium = (mvaValue > 0.885);
     	 isTight = (mvaValue  > 0.986);
     	 
     
       }
       else {
     	 isLoose = (mvaValue  > -0.797);
     	 isMedium = (mvaValue > 0.723);
     	 isTight = (mvaValue  > 0.988);
       }
     }
     else {
       if (not (elecs->at(ie).userFloat("hgcElectronID:ecEnergy") > 0)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:sigmaUU") > 0)) continue;
       if (not (elecs->at(ie).fbrem() > -1)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:measuredDepth") < 40)) continue;
       if (not (elecs->at(ie).userFloat("hgcElectronID:nLayers") > 20)) continue;
       if (elecs->at(ie).pt() < 20.) {
     	 isLoose = (mvaValue  > -0.320);
     	 isMedium = (mvaValue > 0.777);
     	 isTight = (mvaValue  > 0.969);
       }
       else {
     	 isLoose = (mvaValue  > -0.919);
     	 isMedium = (mvaValue > 0.591);
     	 isTight = (mvaValue  > 0.983);
       }
     }

     if(isLoose)
       IDPasselec[Nelec] |= 1 << 0;

     if(isMedium)
       IDPasselec[Nelec] |= 1 << 1;

     if(isTight)
       IDPasselec[Nelec] |= 1 << 2;


     if(Isolationelec[Nelec] < 0.1)
       IsoPasselec[Nelec] |= 1 << 0;
     
     if(Isolationelec[Nelec] < 0.2)
       IsoPasselec[Nelec] |= 1 << 1;
     
     if(Isolationelec[Nelec] < 0.3)
       IsoPasselec[Nelec] |= 1 << 2;
     
     if(Isolationelec[Nelec] < 0.4)
       IsoPasselec[Nelec] |= 1 << 3;
          
     Nelec++;
    if(Nelec>kMaxElectron) break;
   }
   
   //std::cout<<Nelec<<std::endl;
   if(debug_)   std::cout<<"Here I am : got elec infor right "<<std::endl;
   
   /////////////////////////////
   // Muon information
   /////////////////////////////
   for(size_t im= 0 ; im < muons->size(); im++)
     {
    if (muons->at(im).pt() < 2.) continue;
    if (fabs(muons->at(im).eta()) > 2.8) continue;

	   Ptmuon[Nmuon]               = muons->at(im).pt();
	   Etamuon[Nmuon]              = muons->at(im).eta();
	   Phimuon[Nmuon]              = muons->at(im).phi();
	   Massmuon[Nmuon]             = muons->at(im).mass();
	   Chargemuon[Nmuon]           = muons->at(im).charge();
	   Isolationmuon[Nmuon]        = muons->at(im).trackIso()/muons->at(im).pt();
	   IsoPassmuon[Nmuon]          = muons->at(im).mass();
	   IDVarmuon[Nmuon]            = 1.0;
	   double dPhiCut = std::min(std::max(1.2/muons->at(im).p(),1.2/100),0.056);
	   double dPhiBendCut = std::min(std::max(0.2/muons->at(im).p(),0.2/100),0.0096);

	   int isLoose = (int) (fabs(muons->at(im).eta()) < 2.4 && muon::isLooseMuon(muons->at(im))) || (fabs(muons->at(im).eta()) > 2.4 && isME0MuonSelNew(muons->at(im), 0.077, dPhiCut, dPhiBendCut, iSetup));
	   if(isLoose)
	     IDPassmuon[Nmuon] |= 1 << 0;

	   //// isLoose = isMedium
	   int isMedium = isLoose;
	   if(isMedium)
	     IDPassmuon[Nmuon] |= 1 << 1;
	   
	   	   
	   bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
	   if (muons->at(im).innerTrack().isNonnull()){
	     ipxy = std::abs(muons->at(im).muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
	     ipz = std::abs(muons->at(im).muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
	     validPxlHit = muons->at(im).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
	     highPurity = muons->at(im).innerTrack()->quality(reco::Track::highPurity);
	   }      
	   dPhiCut = std::min(std::max(1.2/muons->at(im).p(),1.2/100),0.032);
	   dPhiBendCut = std::min(std::max(0.2/muons->at(im).p(),0.2/100),0.0041);
	   int isTight = (int) (fabs(muons->at(im).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(im),vertices->at(prVtx))) || (fabs(muons->at(im).eta()) > 2.4 && isME0MuonSelNew(muons->at(im), 0.048, dPhiCut, dPhiBendCut, iSetup) && ipxy && ipz && validPxlHit && highPurity);
	   if(isTight)
	     IDPassmuon[Nmuon] |= 1 << 2;
	   
	   
	   
	   //CutBasedIdLoosemuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdLoose);
	   //CutBasedIdMediummuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdMedium);
	   //CutBasedIdTightmuon[Nmuon] = muons->at(im).passed(reco::Muon::CutBasedIdTight);
	   //PFIsoLoosemuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoLoose);
	   //PFIsoMediummuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoMedium);
	   //
	   //PFIsoTightmuon[Nmuon] = muons->at(im).passed(reco::Muon::PFIsoTight);



	   if(Isolationmuon[Nmuon] < 0.1)
	     IsoPassmuon[Nmuon] |= 1 << 0;
	   
	   if(Isolationmuon[Nmuon] < 0.2)
	     IsoPassmuon[Nmuon] |= 1 << 1;
	   
	   if(Isolationmuon[Nmuon] < 0.3)
	     IsoPassmuon[Nmuon] |= 1 << 2;
	   
	   if(Isolationmuon[Nmuon] < 0.4)
	     IsoPassmuon[Nmuon] |= 1 << 3;



	   Nmuon++;
	   if(Nmuon>kMaxMuonLoose) break;
     }
     
   if(debug_)   std::cout<<"Here I am : got muon infor right "<<std::endl;
   //std::cout<<"Muon"<<std::endl;

   /////////////////////////////
   // Taus info
   /////////////////////////////
   for (size_t it = 0; it < taus->size(); it++) {
     if (taus->at(it).pt()<15.) continue; 
     if (fabs(taus->at(it).eta()) > 3.0) continue;
     if (taus->at(it).tauID("decayModeFinding")<0) continue;     
     

     Pttau[Ntau]         = taus->at(it).pt();
     Etatau[Ntau]        = taus->at(it).eta();
     Phitau[Ntau]        = taus->at(it).phi();
     Masstau[Ntau]       = taus->at(it).mass();
     Chargetau[Ntau]     = taus->at(it).charge();
     DMtau[Ntau]         = taus->at(it).decayMode();
     float ChargedIsoTau = taus->at(it).tauID("chargedIsoPtSumdR03");
     float NeutralIsoTau = taus->at(it).tauID("neutralIsoPtSumdR03");
     Isolationtau[Ntau]  = ChargedIsoTau + 0.2*max(0.,NeutralIsoTau - 5.) / taus->at(it).pt();
     IsoFunctiontau[Ntau]= calculate_demetraIsolation(taus->at(it));
     
     if(Isolationtau[Ntau] < 0.1)
       IsoPasstau[Ntau] |= 1 << 0;
     
     if(Isolationtau[Ntau] < 0.2)
       IsoPasstau[Ntau] |= 1 << 1;

     if(Isolationtau[Ntau] < 0.3)
	 IsoPasstau[Ntau] |= 1 << 2;
     
     if(Isolationtau[Ntau] < 0.4)
       IsoPasstau[Ntau] |= 1 << 3;
     

     Ntau++;
     if(Ntau>kMaxTau) break; 
   }
   
   if(debug_)   std::cout<<"Here I am : got tau infor right "<<std::endl;
   
   
   /////////////////////////////
   //Jet information
   /////////////////////////////
   for(size_t ij= 0 ; ij < jets->size(); ij++)
     {
       //if (jets->at(ij).pt() < 20.) continue;
       //if (fabs(jets->at(ij).eta()) > 5) continue;
       //
       Ptjet[Njet]  = jets->at(ij).pt();
       Etajet[Njet] = jets->at(ij).eta();
       Phijet[Njet] = jets->at(ij).phi();
       Massjet[Njet]= jets->at(ij).mass();

       bool isLoose(0), isMedium(0), isTight(0);
       if( (jets->at(ij).numberOfDaughters() > 1 ) && ( jets->at(ij).neutralEmEnergyFraction()< 0.99 ) && ( jets->at(ij).neutralHadronEnergyFraction() < 0.99 ))
       	 {
       	   isLoose = 1;
       	 }
       isMedium = isLoose;
       if( (jets->at(ij).numberOfDaughters() > 1 ) && ( jets->at(ij).neutralEmEnergyFraction()< 0.9 ) && ( jets->at(ij).neutralHadronEnergyFraction() < 0.9 ))
       	 {
	   isTight = 1;
	 }

       if(isLoose)
	 IDPassjet[Njet] |= 1 << 0;
	   
       if(isMedium)
	 IDPassjet[Njet] |= 1 << 1;

       if(isTight)
	 IDPassjet[Njet] |= 1 << 2;
       
       bMVAjet[Njet]         = (float) jets->at(ij).bDiscriminator("pfCombinedMVAV2BJetTags");
       DeepCSVjetTag1[Njet]  = (float) jets->at(ij).bDiscriminator("pfDeepCSVJetTags:probb");
       DeepCSVjetTag2[Njet]  = (float) jets->at(ij).bDiscriminator("pfDeepCSVJetTags:probbb");
          
       if(debug_)
	 {
	   std::cout<<"Btaggers::::::::pfCombinedMVAV2BJetTags/pfDeepCSVJetTags:probb/pfDeepCSVJetTags:probbb:::::::::::::"<<std::endl;
	   std::cout<<bMVAjet[Njet]<<" / "<< DeepCSVjetTag1[Njet] <<" / "<< DeepCSVjetTag2[Njet] <<std::endl;
	 }
       Njet++;
       if(Njet>kMaxJet) break;
     }
   if(debug_)   std::cout<<"Here I am : got jet infor right "<<std::endl;
   

   /////////////////////////////
   //Met information
   /////////////////////////////
   for(size_t imet= 0 ; imet < met->size(); imet++)
     {
       Met[Nmet]    = met->at(imet).pt();
       Phimet[Nmet] = met->at(imet).phi();
       Nmet++;
      }
   
   
   if(debug_)   std::cout<<"Here I am : got met infor right "<<std::endl;
   if(debug_) std::cout<<"beforeend"<<std::endl;
 
   mytree->Fill();
   if(debug_) std::cout<<"end"<<std::endl;
}


// ------------ method to improve ME0 muon ID ----------------
bool 
Validator::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut, edm::EventSetup const& iSetup)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;


    if(!ME0Geometry_)
      //throw std::runtime_error("Validator::isME0MuonSelNew: muon geometry not loaded");
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


////Jan's code//
float Validator::calculate_demetraIsolation(const pat::Tau& tau)const{
  /*  unsigned int tau_vertex_idxpf=-1;
  pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
  tau_vertex_idxpf = packedLeadTauCand->vertexRef().key();
  
  float isoDR03pt08dz015=0;
  float gamma_DR03sum=0;
  
  
  for(const auto& IsoCand: tau.isolationChargedHadrCands()){
    pat::PackedCandidate const* cand = dynamic_cast<pat::PackedCandidate const*>(IsoCand.get());
    if (! cand->charge() )continue;
    //WATCH OUT WHICH VERTICES THESE ARE
    /*if(!vertices())continue;
    const auto& tau_vertex = (*vertices())[tau_vertex_idxpf];
    
    if ((cand->pt()<=0.8) || (fabs(cand->dxy(tau_vertex.position()))>=0.05))continue;
    if (cand->hasTrackDetails()){
      const auto &tt = cand->pseudoTrack();
      if (tt.normalizedChi2()>=100. || cand->numberOfHits()<3)continue;
    }
    
    if (reco::deltaR2(tau,*cand)<0.3*0.3
	&& fabs(cand->dz(tau_vertex.position()))<0.15){
	   
	   isoDR03pt08dz015+=cand->pt();
    }
  }
  for(const auto&  IsoCand: tau.isolationGammaCands()){
    pat::PackedCandidate const* cand = dynamic_cast<pat::PackedCandidate const*>(IsoCand.get());
    if ( cand->pt() < 0.5 ) continue;
    if (reco::deltaR2(tau,*cand)<0.3*0.3 && cand->pt()>1.){
      gamma_DR03sum+=cand->pt();
      }
  }
  */
  //return (isoDR03pt08dz015 + 0.2 * std::max(0., gamma_DR03sum - 5.))/tau.pt();
  return 1.;
    
  
}





// ------------ method called once each job just before starting event loop  ------------
void 
Validator::beginJob()
{

}


void
Validator::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
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
  Validator::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }


// ------------ method called once each job just after ending the event loop  ------------
void 
Validator::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Validator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Validator);

