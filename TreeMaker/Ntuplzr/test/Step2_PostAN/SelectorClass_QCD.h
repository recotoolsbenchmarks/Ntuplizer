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
 

  TH1F *nvtxAll;  
  TH1F *ptAllGenElec,*etaAllGenElec, *phiAllGenElec;
  TH1F *ptAllElec, *etaAllElec, *phiAllElec;
  TH1F *nvtxAllElec;    
  TH1F *rel_isoAllElec_EE,*rel_isoAllElec_EB;  
  TH1F *ptAllElec_EB,*ptAllElec_EE;    

  TH1F *ptTightElec, *etaTightElec, *phiTightElec;
  TH1F *nvtxTightElec; 
  TH1F *ptTightElec_EB,*ptTightElec_EE;
  TH1F *rel_isoTightElec_EE,*rel_isoTightElec_EB;  
  TH1F *ptLooseElec, *etaLooseElec, *phiLooseElec;
  TH1F *nvtxLooseElec;       
  TH1F *ptLooseElec_EB,*ptLooseElec_EE;
  TH1F *rel_isoLooseElec_EE,*rel_isoLooseElec_EB;  
  TH1F *ptMediumElec, *etaMediumElec, *phiMediumElec;
  TH1F *nvtxMediumElec;       
  TH1F *ptMediumElec_EB,*ptMediumElec_EE;
  TH1F *rel_isoMediumElec_EE,*rel_isoMediumElec_EB;  
  TH1F *ptGoodElec, *etaGoodElec, *phiGoodElec;
  TH1F *nvtxGoodElec;    
  TH1F *ptGoodElec_EB,*ptGoodElec_EE;
  TH1F *h_dRElec;
  //*rel_isoAllGenElec;
 
  TH1F *ptAllGenMuon,*etaAllGenMuon, *phiAllGenMuon;
  TH1F *ptAllMuon, *etaAllMuon, *phiAllMuon;
  TH1F *nvtxAllMuon;    
  TH1F *rel_isoAllMuon;  
  TH1F *rel_isoAllMuon_EE,*rel_isoAllMuon_EB;  
  TH1F *ptAllMuon_EB,*ptAllMuon_EE;    

  TH1F *ptTightMuon, *etaTightMuon, *phiTightMuon;
  TH1F *nvtxTightMuon; 
  TH1F *ptTightMuon_EB,*ptTightMuon_EE;
  TH1F *rel_isoTightMuon;  
  TH1F *rel_isoTightMuon_EE,*rel_isoTightMuon_EB;  
  TH1F *ptLooseMuon, *etaLooseMuon, *phiLooseMuon;
  TH1F *nvtxLooseMuon;       
  TH1F *rel_isoLooseMuon;  
  TH1F *rel_isoLooseMuon_EE,*rel_isoLooseMuon_EB;  
  TH1F *ptLooseMuon_EB,*ptLooseMuon_EE;
  TH1F *ptGoodMuon, *etaGoodMuon, *phiGoodMuon;
  TH1F *nvtxGoodMuon;    
  TH1F *ptGoodMuon_EB,*ptGoodMuon_EE;
  TH1F *h_dRMuon;
    //*rel_isoAllGenMuon;
  
  
  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> Nvtx = {fReader, "Nvtx"};
  TTreeReaderArray<Float_t> Vtx_pt2 = {fReader, "Vtx_pt2"};
  TTreeReaderValue<Int_t> Nevt = {fReader, "Nevt"};
  TTreeReaderArray<Float_t> genWeight = {fReader, "genWeight"};
  TTreeReaderValue<Int_t> Ngenlepton = {fReader, "Ngenlepton"};
  TTreeReaderArray<Int_t> PIdgenlepton = {fReader, "PIdgenlepton"};
  TTreeReaderArray<Int_t> Chargegenlepton = {fReader, "Chargegenlepton"};
  TTreeReaderArray<Int_t> Statusgenlepton = {fReader, "Statusgenlepton"};
  TTreeReaderArray<Float_t> Pgenlepton = {fReader, "Pgenlepton"};
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
  TTreeReaderArray<Float_t> Pgenphoton = {fReader, "Pgenphoton"};
  TTreeReaderArray<Float_t> Pxgenphoton = {fReader, "Pxgenphoton"};
  TTreeReaderArray<Float_t> Pygenphoton = {fReader, "Pygenphoton"};
  TTreeReaderArray<Float_t> Pzgenphoton = {fReader, "Pzgenphoton"};
  TTreeReaderArray<Float_t> Egenphoton = {fReader, "Egenphoton"};
  TTreeReaderArray<Float_t> Ptgenphoton = {fReader, "Ptgenphoton"};
  TTreeReaderArray<Float_t> Etagenphoton = {fReader, "Etagenphoton"};
  TTreeReaderArray<Float_t> Phigenphoton = {fReader, "Phigenphoton"};
  TTreeReaderValue<Int_t> Ngenjet = {fReader, "Ngenjet"};
  TTreeReaderArray<Float_t> Ptgenjet = {fReader, "Ptgenjet"};
  TTreeReaderArray<Float_t> Etagenjet = {fReader, "Etagenjet"};
  TTreeReaderArray<Float_t> Phigenjet = {fReader, "Phigenjet"};
  TTreeReaderArray<Float_t> Massgenjet = {fReader, "Massgenjet"};
  TTreeReaderValue<Int_t> Npho = {fReader, "Npho"};
  TTreeReaderArray<Float_t> Ptpho = {fReader, "Ptpho"};
  TTreeReaderArray<Float_t> Etapho = {fReader, "Etapho"};
  TTreeReaderArray<Float_t> Phipho = {fReader, "Phipho"};
  TTreeReaderArray<Float_t> IsolationVarpho = {fReader, "IsolationVarpho"};
  TTreeReaderArray<Float_t> Pxpho = {fReader, "Pxpho"};
  TTreeReaderArray<Float_t> Pypho = {fReader, "Pypho"};
  TTreeReaderArray<Float_t> Pzpho = {fReader, "Pzpho"};
  TTreeReaderArray<Float_t> Epho = {fReader, "Epho"};
  TTreeReaderArray<Int_t> isLpho = {fReader, "isLpho"};
  TTreeReaderArray<Int_t> isTpho = {fReader, "isTpho"};
  TTreeReaderValue<Int_t> Nelec = {fReader, "Nelec"};
  TTreeReaderArray<Float_t> Ptelec = {fReader, "Ptelec"};
  TTreeReaderArray<Float_t> Etaelec = {fReader, "Etaelec"};
  TTreeReaderArray<Float_t> Phielec = {fReader, "Phielec"};
  TTreeReaderArray<Int_t> Chargeelec = {fReader, "Chargeelec"};
  TTreeReaderArray<Float_t> IsolationVarelec = {fReader, "IsolationVarelec"};
  TTreeReaderArray<Float_t> OldIsolationVarelec = {fReader, "OldIsolationVarelec"};
  TTreeReaderArray<Float_t> Pxelec = {fReader, "Pxelec"};
  TTreeReaderArray<Float_t> Pyelec = {fReader, "Pyelec"};
  TTreeReaderArray<Float_t> Pzelec = {fReader, "Pzelec"};
  TTreeReaderArray<Float_t> Eelec = {fReader, "Eelec"};
  TTreeReaderArray<Int_t> isLE = {fReader, "isLE"};
  TTreeReaderArray<Int_t> isME = {fReader, "isME"};
  TTreeReaderArray<Int_t> isTE = {fReader, "isTE"};
  TTreeReaderArray<Float_t> BDTScore = {fReader, "BDTScore"};
  TTreeReaderValue<Int_t> Nmuon = {fReader, "Nmuon"};
  TTreeReaderArray<Float_t> Ptmuon = {fReader, "Ptmuon"};
  TTreeReaderArray<Float_t> Etamuon = {fReader, "Etamuon"};
  TTreeReaderArray<Float_t> Phimuon = {fReader, "Phimuon"};
  TTreeReaderArray<Int_t> Chargemuon = {fReader, "Chargemuon"};
  TTreeReaderArray<Float_t> IsolationVarmuon = {fReader, "IsolationVarmuon"};
  TTreeReaderArray<Float_t> Emuon = {fReader, "Emuon"};
  TTreeReaderArray<Int_t> isLM = {fReader, "isLM"};
  TTreeReaderArray<Int_t> isTM = {fReader, "isTM"};
  TTreeReaderValue<Int_t> Njet = {fReader, "Njet"};
  TTreeReaderArray<Float_t> Ptjet = {fReader, "Ptjet"};
  TTreeReaderArray<Float_t> Etajet = {fReader, "Etajet"};
  TTreeReaderArray<Float_t> Phijet = {fReader, "Phijet"};
  TTreeReaderArray<Int_t> ovEl_ak = {fReader, "ovEl_ak"};
  TTreeReaderArray<Int_t> ovMu_ak = {fReader, "ovMu_ak"};
  TTreeReaderArray<Int_t> Lj = {fReader, "Lj"};
  TTreeReaderArray<Int_t> Tj = {fReader, "Tj"};
  TTreeReaderArray<Int_t> nc = {fReader, "nc"};
  TTreeReaderArray<Float_t> nhef = {fReader, "nhef"};
  TTreeReaderArray<Float_t> neef = {fReader, "neef"};
  TTreeReaderArray<Float_t> chef = {fReader, "chef"};
  TTreeReaderArray<Float_t> ceef = {fReader, "ceef"};
  TTreeReaderValue<Int_t> Nmet = {fReader, "Nmet"};
  TTreeReaderArray<Float_t> Met = {fReader, "Met"};
  TTreeReaderArray<Float_t> Phimet = {fReader, "Phimet"};


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
