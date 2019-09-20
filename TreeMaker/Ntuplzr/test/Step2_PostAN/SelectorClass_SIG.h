//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 28 05:36:43 2017 by ROOT version 6.10/05
// from TTree mytree/TestTree
// found on file: file.root
//////////////////////////////////////////////////////////

#ifndef SelectorClass_h
#define SelectorClass_h
#include "lepInfo.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TProofOutputFile.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TBranch.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include <TMap.h>
#include "TObject.h"
#include <vector>

using namespace std;
float Pi = 3.1415 ;
// Headers needed by this particular selector


class SelectorClass : public TSelector {
 public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
  TProofOutputFile *fProofFile; // For optimized merging of the ntuple 
  TFile *fileName;   
  string weight,OUTFILENAME;
  float wgt;
  
  std::vector<int> jet_index;
  std::vector<int> matchedjet_index;
  //std::vector<int> unmatchedjet_index;

  std::vector<int> elec_index;
  std::vector<int> matchedelec_index;
  std::vector<int> matchedLooseelec_index;
  std::vector<int> matchedMediumelec_index;
  std::vector<int> matchedTightelec_index;
  std::vector<int> matchedGoodelec_index;

  std::vector<int> muon_index;
  std::vector<int> matchedmuon_index;
  std::vector<int> matchedLoosemuon_index;
  std::vector<int> matchedTightmuon_index;
  std::vector<int> matchedGoodmuon_index;

  TH1F *nvtxAll;  
  TH1F *rel_isoAllElec_EE,*rel_isoAllElec_EB;  
  TH1F *ptAllGenElec,*etaAllGenElec, *phiAllGenElec;
  
  TH1F *ptAllMatchedElec, *etaAllMatchedElec, *phiAllMatchedElec;
  TH1F *nvtxAllMatchedElec;    
  TH1F *rel_isoAllMatchedElec_EE,*rel_isoAllMatchedElec_EB;  
  TH1F *ptAllMatchedElec_EB,*ptAllMatchedElec_EE;    

  TH1F *ptTightMatchedElec, *etaTightMatchedElec, *phiTightMatchedElec;
  TH1F *nvtxTightMatchedElec; 
  TH1F *ptTightMatchedElec_EB,*ptTightMatchedElec_EE;
  TH1F *ptLooseMatchedElec, *etaLooseMatchedElec, *phiLooseMatchedElec;
  TH1F *nvtxLooseMatchedElec;       
  TH1F *ptLooseMatchedElec_EB,*ptLooseMatchedElec_EE;
  TH1F *rel_isoLooseMatchedElec_EE,*rel_isoLooseMatchedElec_EB;  
  TH1F *ptMediumMatchedElec, *etaMediumMatchedElec, *phiMediumMatchedElec;
  TH1F *nvtxMediumMatchedElec;       
  TH1F *ptMediumMatchedElec_EB,*ptMediumMatchedElec_EE;
  TH1F *rel_isoMediumMatchedElec_EE,*rel_isoMediumMatchedElec_EB;  
  TH1F *rel_isoTightMatchedElec_EE,*rel_isoTightMatchedElec_EB;
  TH1F *ptGoodMatchedElec, *etaGoodMatchedElec, *phiGoodMatchedElec;
  TH1F *nvtxGoodMatchedElec;    
  TH1F *ptGoodMatchedElec_EB,*ptGoodMatchedElec_EE;
  TH1F *h_dRElec, *rel_isoAllGenElec;


  TH1F *ptAllGenMuon,*etaAllGenMuon, *phiAllGenMuon;
  TH1F *ptAllMatchedMuon, *etaAllMatchedMuon, *phiAllMatchedMuon;
  TH1F *nvtxAllMatchedMuon;    
  TH1F *rel_isoAllMatchedMuon_EE,*rel_isoAllMatchedMuon_EB;  
  TH1F *rel_isoAllMatchedMuon;  
  TH1F *ptAllMatchedMuon_EB,*ptAllMatchedMuon_EE;    

  TH1F *ptTightMatchedMuon, *etaTightMatchedMuon, *phiTightMatchedMuon;
  TH1F *nvtxTightMatchedMuon; 
  TH1F *ptTightMatchedMuon_EB,*ptTightMatchedMuon_EE;
  TH1F *rel_isoTightMatchedMuon;  
  TH1F *rel_isoTightMatchedMuon_EE,*rel_isoTightMatchedMuon_EB;  
  TH1F *ptLooseMatchedMuon, *etaLooseMatchedMuon, *phiLooseMatchedMuon;
  TH1F *nvtxLooseMatchedMuon;       
  TH1F *ptLooseMatchedMuon_EB,*ptLooseMatchedMuon_EE;
  TH1F *rel_isoLooseMatchedMuon;  
  TH1F *rel_isoLooseMatchedMuon_EE,*rel_isoLooseMatchedMuon_EB;  
  TH1F *ptGoodMatchedMuon, *etaGoodMatchedMuon, *phiGoodMatchedMuon;
  TH1F *nvtxGoodMatchedMuon;    
  TH1F *ptGoodMatchedMuon_EB,*ptGoodMatchedMuon_EE;
  TH1F *h_dRMuon, *rel_isoAllGenMuon;

  
  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> Nvtx = {fReader, "Nvtx"};
  TTreeReaderArray<Float_t> Vtx_pt2 = {fReader, "Vtx_pt2"};
  TTreeReaderValue<Int_t> Nevt = {fReader, "Nevt"};
  TTreeReaderArray<Float_t> genWeight = {fReader, "genWeight"};
   TTreeReaderValue<Int_t> Ngenparticle = {fReader, "Ngenparticle"};
   TTreeReaderArray<Int_t> PIdgenparticle = {fReader, "PIdgenparticle"};
   TTreeReaderArray<Int_t> Statusgenparticle = {fReader, "Statusgenparticle"};
   TTreeReaderArray<Float_t> Ptgenparticle = {fReader, "Ptgenparticle"};
   TTreeReaderArray<Float_t> Etagenparticle = {fReader, "Etagenparticle"};
   TTreeReaderArray<Float_t> Phigenparticle = {fReader, "Phigenparticle"};
   TTreeReaderArray<Float_t> Massgenparticle = {fReader, "Massgenparticle"};
   TTreeReaderArray<Int_t> M1genparticle = {fReader, "M1genparticle"};
   TTreeReaderArray<Int_t> M2genparticle = {fReader, "M2genparticle"};
   TTreeReaderArray<Int_t> D1genparticle = {fReader, "D1genparticle"};
   TTreeReaderArray<Int_t> D2genparticle = {fReader, "D2genparticle"};

  TTreeReaderValue<Int_t> Ngenlepton = {fReader, "Ngenlepton"};
  TTreeReaderArray<Int_t> PIdgenlepton = {fReader, "PIdgenlepton"};
  TTreeReaderArray<Int_t> Chargegenlepton = {fReader, "Chargegenlepton"};
  TTreeReaderArray<Int_t> Statusgenlepton = {fReader, "Statusgenlepton"};
   TTreeReaderArray<Float_t> Pxgenlepton = {fReader, "Pxgenlepton"};
   TTreeReaderArray<Float_t> Pygenlepton = {fReader, "Pygenlepton"};
   TTreeReaderArray<Float_t> Pzgenlepton = {fReader, "Pzgenlepton"};
   TTreeReaderArray<Float_t> Egenlepton = {fReader, "Egenlepton"};
  TTreeReaderArray<Float_t> Ptgenlepton = {fReader, "Ptgenlepton"};
  TTreeReaderArray<Float_t> Etagenlepton = {fReader, "Etagenlepton"};
  TTreeReaderArray<Float_t> Phigenlepton = {fReader, "Phigenlepton"};
   TTreeReaderArray<Float_t> Massgenlepton = {fReader, "Massgenlepton"};
  TTreeReaderArray<Float_t> IsolationVargenlepton = {fReader, "IsolationVargenlepton"};
  TTreeReaderValue<Int_t> Ngenphoton = {fReader, "Ngenphoton"};
  TTreeReaderArray<Int_t> Statusgenphoton = {fReader, "Statusgenphoton"};
  TTreeReaderArray<Float_t> Ptgenphoton = {fReader, "Ptgenphoton"};
  TTreeReaderArray<Float_t> Etagenphoton = {fReader, "Etagenphoton"};
  TTreeReaderArray<Float_t> Phigenphoton = {fReader, "Phigenphoton"};
  TTreeReaderValue<Int_t> Ngenjet = {fReader, "Ngenjet"};
  TTreeReaderArray<Float_t> Ptgenjet = {fReader, "Ptgenjet"};
  TTreeReaderArray<Float_t> Etagenjet = {fReader, "Etagenjet"};
  TTreeReaderArray<Float_t> Phigenjet = {fReader, "Phigenjet"};
  TTreeReaderValue<Int_t> Npho = {fReader, "Npho"};
  TTreeReaderArray<Float_t> Ptpho = {fReader, "Ptpho"};
  TTreeReaderArray<Float_t> Etapho = {fReader, "Etapho"};
  TTreeReaderArray<Float_t> Phipho = {fReader, "Phipho"};
  TTreeReaderArray<Float_t> MVApho = {fReader, "MVApho"}; 
  TTreeReaderArray<Float_t> IsolationVarpho = {fReader, "IsolationVarpho"};
  TTreeReaderArray<Int_t> isLpho = {fReader, "isLpho"};
  TTreeReaderArray<Int_t> isTpho = {fReader, "isTpho"};
   TTreeReaderValue<Int_t> Nelec = {fReader, "Nelec"};
   TTreeReaderArray<Float_t> Ptelec = {fReader, "Ptelec"};
   TTreeReaderArray<Float_t> Etaelec = {fReader, "Etaelec"};
   TTreeReaderArray<Float_t> Phielec = {fReader, "Phielec"};
   TTreeReaderArray<Int_t> Chargeelec = {fReader, "Chargeelec"};
   TTreeReaderArray<Float_t> IsolationVarelec = {fReader, "IsolationVarelec"};
   TTreeReaderArray<Float_t> OldIsolationVarelec = {fReader, "OldIsolationVarelec"};
   TTreeReaderArray<Int_t> MissingHitselec = {fReader, "MissingHitselec"};
   TTreeReaderArray<Int_t> ConvVetoElec = {fReader, "ConvVetoElec"};
   TTreeReaderArray<Int_t> ChargeConsistentelec = {fReader, "ChargeConsistentelec"};
   TTreeReaderArray<Int_t> ChargeConsistent1elec = {fReader, "ChargeConsistent1elec"};
   TTreeReaderArray<Int_t> ChargeConsistent2elec = {fReader, "ChargeConsistent2elec"};
   TTreeReaderArray<Float_t> Pxelec = {fReader, "Pxelec"};
   TTreeReaderArray<Float_t> Pyelec = {fReader, "Pyelec"};
   TTreeReaderArray<Float_t> Pzelec = {fReader, "Pzelec"};
   TTreeReaderArray<Float_t> Eelec = {fReader, "Eelec"};
   TTreeReaderArray<Float_t> Masselec = {fReader, "Masselec"};
   TTreeReaderArray<Int_t> isLE = {fReader, "isLE"};
   TTreeReaderArray<Int_t> isME = {fReader, "isME"};
   TTreeReaderArray<Int_t> isTE = {fReader, "isTE"};
   TTreeReaderArray<Float_t> MVAelec = {fReader, "MVAelec"};
   TTreeReaderArray<Int_t> isEBelec = {fReader, "isEBelec"};
   TTreeReaderArray<Float_t> fbremelec = {fReader, "fbremelec"};
   TTreeReaderArray<Float_t> kfhitselec = {fReader, "kfhitselec"};
   TTreeReaderArray<Float_t> kfchi2elec = {fReader, "kfchi2elec"};
   TTreeReaderArray<Float_t> gsfhitselec = {fReader, "gsfhitselec"};
   TTreeReaderArray<Float_t> gsfchi2elec = {fReader, "gsfchi2elec"};
   TTreeReaderArray<Float_t> eelepoutelec = {fReader, "eelepoutelec"};
   TTreeReaderArray<Float_t> deltaetaelec = {fReader, "deltaetaelec"};
   TTreeReaderArray<Float_t> deltaphielec = {fReader, "deltaphielec"};
   TTreeReaderArray<Float_t> oldsigmaietaietaelec = {fReader, "oldsigmaietaietaelec"};
   TTreeReaderArray<Float_t> oldsigmaiphiiphielec = {fReader, "oldsigmaiphiiphielec"};
   TTreeReaderArray<Float_t> oldcircularityelec = {fReader, "oldcircularityelec"};
   TTreeReaderArray<Float_t> oldr9elec = {fReader, "oldr9elec"};
   TTreeReaderArray<Float_t> scletawidthelec = {fReader, "scletawidthelec"};
   TTreeReaderArray<Float_t> sclphiwidthelec = {fReader, "sclphiwidthelec"};
   TTreeReaderArray<Float_t> heelec = {fReader, "heelec"};
   TTreeReaderArray<Float_t> epelec = {fReader, "epelec"};
   TTreeReaderArray<Float_t> eseedpoutelec = {fReader, "eseedpoutelec"};
   TTreeReaderArray<Float_t> deltaetaseedelec = {fReader, "deltaetaseedelec"};
   TTreeReaderArray<Float_t> deltaphiseedelec = {fReader, "deltaphiseedelec"};
   TTreeReaderValue<Int_t> Nmuon = {fReader, "Nmuon"};
   TTreeReaderArray<Float_t> Ptmuon = {fReader, "Ptmuon"};
   TTreeReaderArray<Float_t> Etamuon = {fReader, "Etamuon"};
   TTreeReaderArray<Float_t> Phimuon = {fReader, "Phimuon"};
   TTreeReaderArray<Int_t> Chargemuon = {fReader, "Chargemuon"};
   TTreeReaderArray<Float_t> IsolationVarmuon = {fReader, "IsolationVarmuon"};
   TTreeReaderArray<Float_t> CutBasedIdLoosemuon = {fReader, "CutBasedIdLoosemuon"};
   TTreeReaderArray<Float_t> CutBasedIdMediummuon = {fReader, "CutBasedIdMediummuon"};
   TTreeReaderArray<Float_t> CutBasedIdTightmuon = {fReader, "CutBasedIdTightmuon"};
   TTreeReaderArray<Float_t> PFIsoLoosemuon = {fReader, "PFIsoLoosemuon"};
   TTreeReaderArray<Float_t> PFIsoMediummuon = {fReader, "PFIsoMediummuon"};
   TTreeReaderArray<Float_t> PFIsoTightmuon = {fReader, "PFIsoTightmuon"};
   TTreeReaderArray<Float_t> Emuon = {fReader, "Emuon"};
   TTreeReaderArray<Float_t> Massmuon = {fReader, "Massmuon"};
   TTreeReaderArray<Int_t> isLM = {fReader, "isLM"};
   TTreeReaderArray<Int_t> isTM = {fReader, "isTM"};
   TTreeReaderValue<Int_t> Njet = {fReader, "Njet"};
   TTreeReaderArray<Float_t> Ptjet = {fReader, "Ptjet"};
   TTreeReaderArray<Float_t> Etajet = {fReader, "Etajet"};
   TTreeReaderArray<Float_t> Phijet = {fReader, "Phijet"};
   TTreeReaderArray<Int_t> Lj = {fReader, "Lj"};
   TTreeReaderArray<Int_t> Tj = {fReader, "Tj"};
   TTreeReaderArray<Float_t> Deepjet = {fReader, "Deepjet"};
   TTreeReaderArray<Float_t> Massjet = {fReader, "Massjet"};
   TTreeReaderArray<Int_t> NoDjet = {fReader, "NoDjet"};
   TTreeReaderArray<Float_t> NHEFjet = {fReader, "NHEFjet"};
   TTreeReaderArray<Float_t> NEEFjet = {fReader, "NEEFjet"};
   TTreeReaderArray<Float_t> CHEFjet = {fReader, "CHEFjet"};
   TTreeReaderArray<Float_t> CEMEFjet = {fReader, "CEMEFjet"};
   TTreeReaderArray<Float_t> JECFjet = {fReader, "JECFjet"};
   TTreeReaderValue<Int_t> Nmet = {fReader, "Nmet"};
   TTreeReaderArray<Float_t> Met = {fReader, "Met"};
   TTreeReaderArray<Float_t> Phimet = {fReader, "Phimet"};
   TTreeReaderValue<Int_t> Ntau = {fReader, "Ntau"};
   TTreeReaderArray<Float_t> Pttau = {fReader, "Pttau"};
   TTreeReaderArray<Float_t> Etatau = {fReader, "Etatau"};
   TTreeReaderArray<Float_t> Phitau = {fReader, "Phitau"};
   TTreeReaderArray<Int_t> Chargetau = {fReader, "Chargetau"};
   TTreeReaderArray<Float_t> Masstau = {fReader, "Masstau"};
   TTreeReaderArray<Float_t> DMtau = {fReader, "DMtau"};
   TTreeReaderArray<Float_t> ChargedIsotau = {fReader, "ChargedIsotau"};
   TTreeReaderArray<Float_t> NeutralIsotau = {fReader, "NeutralIsotau"};
   TTreeReaderArray<Float_t> CombinedIsotau = {fReader, "CombinedIsotau"
};



 SelectorClass(TTree * /*tree*/ =0) : fChain(0) {fProofFile =0;} 
  virtual ~SelectorClass() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  virtual float    dphi(float phi1, float phi2);
  virtual float    dR(float eta1, float phi1, float eta2, float phi2);
  virtual float    Minv(float pt0, float eta0, float phi0, float pt1, float eta1, float phi1);
  
  TLorentzVector FourVector(float pt0, float eta0, float phi0);
  ClassDef(SelectorClass,0);

};

#endif

#ifdef SelectorClass_cxx
void SelectorClass::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  fReader.SetTree(tree);
}

Bool_t SelectorClass::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


float SelectorClass::dphi(float phi1, float phi2)
{
  float deltaphi = fabs(phi1-phi2);
  if(deltaphi > Pi) deltaphi = fabs(deltaphi - 2*Pi);
  return deltaphi;
}

float SelectorClass::dR(float eta1, float phi1, float eta2, float phi2)
{
  float dphi = fabs(phi1-phi2);
  float deta = fabs(eta1-eta2);
  if(dphi > Pi) dphi = fabs(dphi - 2*Pi);
  float deltaR   = sqrt(pow(deta,2) + pow(dphi,2));
  return deltaR;
}

float SelectorClass::Minv(float pt0, float eta0, float phi0, float pt1, float eta1, float phi1){
  float E1   = pt1*cosh(eta1);
  float E0   = pt0*cosh(eta0);
  
  float Pz1  = pt1*sinh(eta1); 
  float Pz0  = pt0*sinh(eta0); 
  
  float Px1  = pt1*cos(phi1); 
  float Px0  = pt0*cos(phi0); 

  float Py1  = pt1*sin(phi1); 
  float Py0  = pt0*sin(phi0); 
	 
  float invMass = sqrt(pow(E1 + E0,2)- pow(Px1 + Px0,2)-pow(Py1+Py0,2)-pow(Pz1+Pz0,2));

  return invMass;
}

TLorentzVector SelectorClass::FourVector(float pt0, float eta0, float phi0){
  float E0   = pt0*cosh(eta0);
  float Pz0  = pt0*sinh(eta0); 
  float Px0  = pt0*cos(phi0); 
  float Py0  = pt0*sin(phi0); 

  TLorentzVector p(Px0, Py0, Pz0, E0);
  return p;
}


#endif // #ifdef SelectorClass_cxx
