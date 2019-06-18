#define SelectorClass_cxx
#include "SelectorClass.h"
#include <TH2.h>
#include <TStyle.h>
class PtSortCriterium{
public:
  bool operator() (lepInfo p1,lepInfo p2 ){
    return p1.pt > p2.pt;
  }
};


void SelectorClass::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void SelectorClass::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
   TObjArray *args = (TObjArray*)option.Tokenize(" ");
   if (args->GetSize()>1)
     {
       OUTFILENAME = (string)((TObjString*)args->At(0))->GetString();
       weight = (string)((TObjString*)args->At(1))->GetString();
     }
   std::cout<<"outfile name and weight : "<<OUTFILENAME<<" "<<weight<<std::endl;
   wgt = std::stof(weight);
   std::cout<<"weight in float : "<<wgt<<std::endl;
   
   fProofFile = new TProofOutputFile(OUTFILENAME.c_str());
   fOutput->Add(fProofFile);
   
   fileName = fProofFile->OpenFile("RECREATE");
   fileName->cd();
  
   char name[100];

   sprintf(name,"nvtxAll");
   nvtxAll = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxAll->Sumw2();
   nvtxAll->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptAllGenElec");
   ptAllGenElec = new TH1F (name,"" ,  100, 0, 1000);
   ptAllGenElec->Sumw2();
   ptAllGenElec->GetXaxis()->SetTitle("P_{t}(GeV)");
   //ptAllGenElec->GetNbinsX();
   //ptAllGenElec->GetYaxis()->SetTitle("#Events/GeV");
     
   sprintf(name,"etaAllGenElec");
   etaAllGenElec = new TH1F (name,"" , 100, -5, 5);
   etaAllGenElec->Sumw2();
   etaAllGenElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllGenElec");
   phiAllGenElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllGenElec->Sumw2();
   phiAllGenElec->GetXaxis()->SetTitle("#phi");

   //sprintf(name,"rel_isoAllGenElec");
   //rel_isoAllGenElec = new TH1F (name,"" , 100, 0, 1);
   //rel_isoAllGenElec->Sumw2();
   //rel_isoAllGenElec->GetXaxis()->SetTitle("relIso_{e}");

   sprintf(name,"ptAllElec");
   ptAllElec = new TH1F (name,"" ,  100, 0, 1000);
   ptAllElec->Sumw2();
   ptAllElec->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"etaAllElec");
   etaAllElec = new TH1F (name,"" , 100, -5, 5);
   etaAllElec->Sumw2();
   etaAllElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllElec");
   phiAllElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllElec->Sumw2();
   phiAllElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptAllElec_EB");
   ptAllElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptAllElec_EB->Sumw2();
   ptAllElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllElec_EE");
   ptAllElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptAllElec_EE->Sumw2();
   ptAllElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoAllElec_EB");
   rel_isoAllElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllElec_EB->Sumw2();
   rel_isoAllElec_EB->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoAllElec_EE");
   rel_isoAllElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllElec_EE->Sumw2();
   rel_isoAllElec_EE->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"nvtxAllElec");
   nvtxAllElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxAllElec->Sumw2();
   nvtxAllElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptLooseElec");
   ptLooseElec = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseElec->Sumw2();
   ptLooseElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaLooseElec");
   etaLooseElec = new TH1F (name,"" , 100, -5, 5);
   etaLooseElec->Sumw2();
   etaLooseElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiLooseElec");
   phiLooseElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiLooseElec->Sumw2();
   phiLooseElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptLooseElec_EB");
   ptLooseElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseElec_EB->Sumw2();
   ptLooseElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptLooseElec_EE");
   ptLooseElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseElec_EE->Sumw2();
   ptLooseElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoLooseElec_EB");
   rel_isoLooseElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseElec_EB->Sumw2();
   rel_isoLooseElec_EB->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoLooseElec_EE");
   rel_isoLooseElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseElec_EE->Sumw2();
   rel_isoLooseElec_EE->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"nvtxLooseElec");
   nvtxLooseElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxLooseElec->Sumw2();
   nvtxLooseElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptMediumElec");
   ptMediumElec = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumElec->Sumw2();
   ptMediumElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaMediumElec");
   etaMediumElec = new TH1F (name,"" , 100, -5, 5);
   etaMediumElec->Sumw2();
   etaMediumElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiMediumElec");
   phiMediumElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiMediumElec->Sumw2();
   phiMediumElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptMediumElec_EB");
   ptMediumElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumElec_EB->Sumw2();
   ptMediumElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptMediumElec_EE");
   ptMediumElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumElec_EE->Sumw2();
   ptMediumElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoMediumElec_EB");
   rel_isoMediumElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoMediumElec_EB->Sumw2();
   rel_isoMediumElec_EB->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoMediumElec_EE");
   rel_isoMediumElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoMediumElec_EE->Sumw2();
   rel_isoMediumElec_EE->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"nvtxMediumElec");
   nvtxMediumElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxMediumElec->Sumw2();
   nvtxMediumElec->GetXaxis()->SetTitle("#vertices");
   
   sprintf(name,"ptTightElec");
   ptTightElec = new TH1F (name,"" ,  100, 0, 1000);
   ptTightElec->Sumw2();
   ptTightElec->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"etaTightElec");
   etaTightElec = new TH1F (name,"" , 100, -5, 5);
   etaTightElec->Sumw2();
   etaTightElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiTightElec");
   phiTightElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiTightElec->Sumw2();
   phiTightElec->GetXaxis()->SetTitle("#phi");

   
   sprintf(name,"ptTightElec_EB");
   ptTightElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptTightElec_EB->Sumw2();
   ptTightElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptTightElec_EE");
   ptTightElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptTightElec_EE->Sumw2();
   ptTightElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoTightElec_EB");
   rel_isoTightElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightElec_EB->Sumw2();
   rel_isoTightElec_EB->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoTightElec_EE");
   rel_isoTightElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightElec_EE->Sumw2();
   rel_isoTightElec_EE->GetXaxis()->SetTitle("relIso_{e^{-1}}");
     
   sprintf(name,"nvtxTightElec");
   nvtxTightElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxTightElec->Sumw2();
   nvtxTightElec->GetXaxis()->SetTitle("#vertices");


   sprintf(name,"ptGoodElec");
   ptGoodElec = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodElec->Sumw2();
   ptGoodElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaGoodElec");
   etaGoodElec = new TH1F (name,"" , 100, -5, 5);
   etaGoodElec->Sumw2();
   etaGoodElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiGoodElec");
   phiGoodElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiGoodElec->Sumw2();
   phiGoodElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxGoodElec");
   nvtxGoodElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxGoodElec->Sumw2();
   nvtxGoodElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptGoodElec_EB");
   ptGoodElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodElec_EB->Sumw2();
   ptGoodElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptGoodElec_EE");
   ptGoodElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodElec_EE->Sumw2();
   ptGoodElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllGenMuon");
   ptAllGenMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptAllGenMuon->Sumw2();
   ptAllGenMuon->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"etaAllGenMuon");
   etaAllGenMuon = new TH1F (name,"" , 100, -5, 5);
   etaAllGenMuon->Sumw2();
   etaAllGenMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllGenMuon");
   phiAllGenMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllGenMuon->Sumw2();
   phiAllGenMuon->GetXaxis()->SetTitle("#phi");

   //sprintf(name,"rel_isoAllGenMuon");
   //rel_isoAllGenMuon = new TH1F (name,"" , 100, 0, 1);
   //rel_isoAllGenMuon->Sumw2();
   //rel_isoAllGenMuon->GetXaxis()->SetTitle("relIso");

   sprintf(name,"ptAllMuon");
   ptAllMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMuon->Sumw2();
   ptAllMuon->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"etaAllMuon");
   etaAllMuon = new TH1F (name,"" , 100, -5, 5);
   etaAllMuon->Sumw2();
   etaAllMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllMuon");
   phiAllMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllMuon->Sumw2();
   phiAllMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptAllMuon_EB");
   ptAllMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMuon_EB->Sumw2();
   ptAllMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllMuon_EE");
   ptAllMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMuon_EE->Sumw2();
   ptAllMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoAllMuon");
   rel_isoAllMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMuon ->Sumw2();
   rel_isoAllMuon ->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoAllMuon_EB");
   rel_isoAllMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMuon_EB->Sumw2();
   rel_isoAllMuon_EB->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"rel_isoAllMuon_EE");
   rel_isoAllMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMuon_EE->Sumw2();
   rel_isoAllMuon_EE->GetXaxis()->SetTitle("relIso_{e^{-1}}");

   sprintf(name,"nvtxAllMuon");
   nvtxAllMuon = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxAllMuon->Sumw2();
   nvtxAllMuon->GetXaxis()->SetTitle("#vertices");


   sprintf(name,"ptLooseMuon");
   ptLooseMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMuon->Sumw2();
   ptLooseMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaLooseMuon");
   etaLooseMuon = new TH1F (name,"" , 100, -5, 5);
   etaLooseMuon->Sumw2();
   etaLooseMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiLooseMuon");
   phiLooseMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiLooseMuon->Sumw2();
   phiLooseMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxLooseMuon");
   nvtxLooseMuon = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxLooseMuon->Sumw2();
   nvtxLooseMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptLooseMuon_EB");
   ptLooseMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMuon_EB->Sumw2();
   ptLooseMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptLooseMuon_EE");
   ptLooseMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMuon_EE->Sumw2();
   ptLooseMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoLooseMuon");
   rel_isoLooseMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMuon ->Sumw2();
   rel_isoLooseMuon ->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoLooseMuon_EB");
   rel_isoLooseMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMuon_EB->Sumw2();
   rel_isoLooseMuon_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoLooseMuon_EE");
   rel_isoLooseMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMuon_EE->Sumw2();
   rel_isoLooseMuon_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"ptTightMuon");
   ptTightMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMuon->Sumw2();
   ptTightMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaTightMuon");
   etaTightMuon = new TH1F (name,"" , 100, -5, 5);
   etaTightMuon->Sumw2();
   etaTightMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiTightMuon");
   phiTightMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiTightMuon->Sumw2();
   phiTightMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxTightMuon");
   nvtxTightMuon = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxTightMuon->Sumw2();
   nvtxTightMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptTightMuon_EB");
   ptTightMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMuon_EB->Sumw2();
   ptTightMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptTightMuon_EE");
   ptTightMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMuon_EE->Sumw2();
   ptTightMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoTightMuon");
   rel_isoTightMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMuon ->Sumw2();
   rel_isoTightMuon ->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoTightMuon_EB");
   rel_isoTightMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMuon_EB->Sumw2();
   rel_isoTightMuon_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoTightMuon_EE");
   rel_isoTightMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMuon_EE->Sumw2();
   rel_isoTightMuon_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"ptGoodMuon");
   ptGoodMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMuon->Sumw2();
   ptGoodMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaGoodMuon");
   etaGoodMuon = new TH1F (name,"" , 100, -5, 5);
   etaGoodMuon->Sumw2();
   etaGoodMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiGoodMuon");
   phiGoodMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiGoodMuon->Sumw2();
   phiGoodMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxGoodMuon");
   nvtxGoodMuon = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxGoodMuon->Sumw2();
   nvtxGoodMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptGoodMuon_EB");
   ptGoodMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMuon_EB->Sumw2();
   ptGoodMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptGoodMuon_EE");
   ptGoodMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMuon_EE->Sumw2();
   ptGoodMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");
  
}

Bool_t SelectorClass::Process(Long64_t entry)
{
  fReader.SetEntry(entry);

  nvtxAll -> Fill((float) *Nvtx , wgt);
  /*
  std::vector<lepInfo> lepinfo_container;
  lepinfo_container.clear();

  std::vector<lepInfo> phoinfo_container;
  phoinfo_container.clear();
  lepInfo lepton;


  for(Int_t ip= 0 ; ip < *Npho; ip++) {
    if(Ptpho[ip] < 20. ) continue;
    if(fabs(Etapho[ip]) > 3.0 ) continue;
    if(fabs(Etapho[ip])> 1.444 && fabs(Etapho[ip])< 1.566 ) continue;
    //if(fabs(Etapho[ip]) < 1.444 )
    //  if(MVApho[ip] < 0.75) continue;
    //if(fabs(Etapho[ie]) > 1.566 )
    //  if(MVApho[ip] < 0.8) continue;
    */


  for(Int_t ig= 0 ; ig < *Ngenlepton; ig++) {
    if(abs(PIdgenlepton[ig])==11)
      {
	if (Ptgenlepton[ig] < 20. ) continue;
	if(fabs(Etagenlepton[ig]) > 3.0 ) continue;
	if(fabs(Etagenlepton[ig]) > 1.444 && fabs(Etagenlepton[ig]) < 1.566 ) continue;
	//rel_isoAllGenElec -> Fill(IsolationVargenlepton[ig] , wgt);
	//if (IsolationVargenlepton[ig] > 0.02 ) continue;
	ptAllGenElec  -> Fill(Ptgenlepton[ig]  , wgt);
	etaAllGenElec -> Fill(Etagenlepton[ig] , wgt);
	phiAllGenElec -> Fill(Phigenlepton[ig] , wgt);
      }
  }

  for(Int_t ie= 0 ; ie < *Nelec; ie++) {
    if(Ptelec[ie] < 20. ) continue;
    if(fabs(Etaelec[ie]) > 3.0 ) continue;
    if(fabs(Etaelec[ie])> 1.444 && fabs(Etaelec[ie])< 1.566 ) continue;
    ptAllElec   -> Fill(Ptelec[ie]  , wgt);
    etaAllElec  -> Fill(Etaelec[ie] , wgt);
    phiAllElec  -> Fill(Phielec[ie] , wgt);
    nvtxAllElec -> Fill((float) *Nvtx , wgt);
    if(fabs(Etaelec[ie]) < 1.444)
      {
	rel_isoAllElec_EB -> Fill(IsolationVarelec[ie], wgt);
	ptAllElec_EB  -> Fill(Ptelec[ie]  , wgt);
      }
    if(fabs(Etaelec[ie]) > 1.566)
      {
	rel_isoAllElec_EE -> Fill(IsolationVarelec[ie], wgt);
	ptAllElec_EE  -> Fill(Ptelec[ie]  , wgt);
      }


    if(!isLE[ie]) continue;

    ptLooseElec  -> Fill(Ptelec[ie]  , wgt);
    etaLooseElec -> Fill(Etaelec[ie] , wgt);
    phiLooseElec -> Fill(Phielec[ie] , wgt);
    nvtxLooseElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[ie]) < 1.444)
      {
	ptLooseElec_EB  -> Fill(Ptelec[ie]  , wgt);
	rel_isoLooseElec_EB -> Fill(IsolationVarelec[ie], wgt);
      }
    if(fabs(Etaelec[ie]) > 1.566)
      {
	ptLooseElec_EE  -> Fill(Ptelec[ie]  , wgt);
	rel_isoLooseElec_EE -> Fill(IsolationVarelec[ie], wgt);
      }

    if(!isME[ie]) continue;

    ptMediumElec  -> Fill(Ptelec[ie]  , wgt);
    etaMediumElec -> Fill(Etaelec[ie] , wgt);
    phiMediumElec -> Fill(Phielec[ie] , wgt);
    nvtxMediumElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[ie]) < 1.444)
      {
	ptMediumElec_EB  -> Fill(Ptelec[ie]  , wgt);
	rel_isoMediumElec_EB -> Fill(IsolationVarelec[ie], wgt);
      }
    if(fabs(Etaelec[ie]) > 1.566)
      {
	ptMediumElec_EE  -> Fill(Ptelec[ie]  , wgt);
	rel_isoMediumElec_EE -> Fill(IsolationVarelec[ie], wgt);
      }

    
    if((fabs(Etaelec[ie]) < 1.444 && IsolationVarelec[ie] <  0.15)
       || (fabs(Etaelec[ie]) > 1.566 && IsolationVarelec[ie] < 0.21)) {  
    ptGoodElec  ->Fill(Ptelec[ie] , wgt);
    etaGoodElec ->Fill(Etaelec[ie], wgt);
    phiGoodElec ->Fill(Phielec[ie], wgt);
    nvtxGoodElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[ie]) < 1.444)
	ptGoodElec_EB  -> Fill(Ptelec[ie]  , wgt);
    if(fabs(Etaelec[ie]) > 1.566)
	ptGoodElec_EE  -> Fill(Ptelec[ie]  , wgt);

    }

    if(!isTE[ie]) continue;
    ptTightElec  -> Fill(Ptelec[ie]  , wgt);
    etaTightElec -> Fill(Etaelec[ie] , wgt);
    phiTightElec -> Fill(Phielec[ie] , wgt);
    nvtxTightElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[ie]) < 1.444)
      {
	ptTightElec_EB  -> Fill(Ptelec[ie]  , wgt);
	rel_isoTightElec_EB -> Fill(IsolationVarelec[ie], wgt);
      }
    if(fabs(Etaelec[ie]) > 1.566)
      {
	ptTightElec_EE  -> Fill(Ptelec[ie]  , wgt);
	rel_isoTightElec_EE -> Fill(IsolationVarelec[ie], wgt);
      }

  }

  for(Int_t ig= 0 ; ig < *Ngenlepton; ig++) {
    if(abs(PIdgenlepton[ig])==13)
      {
	if (Ptgenlepton[ig] < 20. ) continue;
	if(fabs(Etagenlepton[ig]) > 3.0 ) continue;
	//if(fabs(Etagenlepton[ig]) > 1.444 && fabs(Etagenlepton[ig]) < 1.566 ) continue;
	//rel_isoAllGenMuon -> Fill(IsolationVargenlepton[ig] , wgt);
	//if (IsolationVargenlepton[ig] > 0.02 ) continue;

	ptAllGenMuon  -> Fill(Ptgenlepton[ig]  , wgt);
	etaAllGenMuon -> Fill(Etagenlepton[ig] , wgt);
	phiAllGenMuon -> Fill(Phigenlepton[ig] , wgt);
      }
  }



  for(Int_t ie= 0 ; ie < *Nmuon; ie++) {
    if(Ptmuon[ie] < 20. ) continue;
    if(fabs(Etamuon[ie]) > 3.0 ) continue;
    
    ptAllMuon -> Fill(Ptmuon[ie]  , wgt);
    etaAllMuon -> Fill(Etamuon[ie] , wgt);
    phiAllMuon -> Fill(Phimuon[ie] , wgt);
    nvtxAllMuon -> Fill((float) *Nvtx , wgt);
    rel_isoAllMuon -> Fill(IsolationVarmuon[ie], wgt);
    if(fabs(Etamuon[ie]) < 1.444)
      {
	rel_isoAllMuon_EB -> Fill(IsolationVarmuon[ie], wgt);
	ptAllMuon_EB  -> Fill(Ptmuon[ie]  , wgt);
      }
    else
      {
	rel_isoAllMuon_EE -> Fill(IsolationVarmuon[ie], wgt);
	ptAllMuon_EE  -> Fill(Ptmuon[ie]  , wgt);
      }


    if(!isLM[ie]) continue;
    ptLooseMuon  -> Fill(Ptmuon[ie]  , wgt);
    etaLooseMuon -> Fill(Etamuon[ie] , wgt);
    phiLooseMuon -> Fill(Phimuon[ie] , wgt);
    nvtxLooseMuon -> Fill((float)*Nvtx , wgt);
    rel_isoLooseMuon -> Fill(IsolationVarmuon[ie], wgt);
    if(fabs(Etamuon[ie]) < 1.444)
      {
	rel_isoLooseMuon_EB -> Fill(IsolationVarmuon[ie], wgt);
	ptLooseMuon_EB  -> Fill(Ptmuon[ie]  , wgt);
      }
    else
      {
	rel_isoLooseMuon_EE -> Fill(IsolationVarmuon[ie], wgt);
	ptLooseMuon_EE  -> Fill(Ptmuon[ie]  , wgt);
      }

    if(!isTM[ie]) continue;


    ptTightMuon  -> Fill(Ptmuon[ie]  , wgt);
    etaTightMuon -> Fill(Etamuon[ie] , wgt);
    phiTightMuon -> Fill(Phimuon[ie] , wgt);
    nvtxTightMuon -> Fill((float)*Nvtx , wgt);
    rel_isoTightMuon -> Fill(IsolationVarmuon[ie], wgt);
    if(fabs(Etamuon[ie]) < 1.444)
      {
	rel_isoTightMuon_EB -> Fill(IsolationVarmuon[ie], wgt);
	ptTightMuon_EB  -> Fill(Ptmuon[ie]  , wgt);
      }
    else
      {
	rel_isoTightMuon_EE -> Fill(IsolationVarmuon[ie], wgt);
	ptTightMuon_EE  -> Fill(Ptmuon[ie]  , wgt);
      }


    if(IsolationVarmuon[ie] < 0.15 ){
    ptGoodMuon  ->Fill(Ptmuon[ie] , wgt);
    etaGoodMuon ->Fill(Etamuon[ie], wgt);
    phiGoodMuon ->Fill(Phimuon[ie], wgt);
    nvtxGoodMuon -> Fill((float)*Nvtx , wgt);
    if(fabs(Etamuon[ie]) < 1.444)
	ptGoodMuon_EB  -> Fill(Ptmuon[ie]  , wgt);
    else
	ptGoodMuon_EE  -> Fill(Ptmuon[ie]  , wgt);
    }
  }

  /*

  jet_index.clear();

  for(Int_t ij= 0 ; ij < *Njet; ij++) {
    if(Ptjet[ij] < 30. ) continue;
    if(fabs(Etajet[ij]) > 4.7 ) continue;
    
    bool overlap=false;
    
    for(int ie=0;ie<*Nelec;ie++)
      {
	float dR_ejt   = dR(Etaelec[ie], Phielec[ie], Etajet[ij], Phijet[ij]);
	float ptdiff   = fabs(Ptjet[ij] - Ptelec[ie]);
	dR_jetlep->Fill(dR_ejt, wgt);
        //pt_ptdiff->Fill(Ptelec[ie], ptdiff);
        //relpt    -> Fill(ptdiff/Ptelec[ie]);
        if(dR_ejt<0.3) {overlap=true; break;}
      }
    if(overlap) continue;
    for(int ie=0;ie<*Nmuon;ie++)
      {
        float dR_ejt   = dR(Etamuon[ie], Phimuon[ie], Etajet[ij], Phijet[ij]);
        float ptdiff   = fabs(Ptjet[ij] - Ptmuon[ie]);
	dR_jetlep->Fill(dR_ejt);
        //pt_ptdiff->Fill(Ptmuon[ie], ptdiff);
        relpt    -> Fill(ptdiff/Ptmuon[ie]);
        if(dR_ejt<0.3) {overlap=true; break;}

      }

    if(overlap) continue;

    ptAllJet  -> Fill(Ptjet[ij]  , wgt);
    etaAllJet -> Fill(Etajet[ij] , wgt);
    phiAllJet -> Fill(Phijet[ij] , wgt);

    if(fabs(Etajet[ij]) < 3.0)
      {
      ptAllJet_HB -> Fill(Ptjet[ij]  , wgt);
      etaAllJet_HB -> Fill(Etajet[ij] , wgt);
      phiAllJet_HB -> Fill(Phijet[ij] , wgt);
      }
    else
      {
      ptAllJet_HE -> Fill(Ptjet[ij]  , wgt);
      etaAllJet_HE -> Fill(Etajet[ij] , wgt);
      phiAllJet_HE -> Fill(Phijet[ij] , wgt);
      }

    if(!Lj[ij]) continue;
    ptLooseJet  -> Fill(Ptjet[ij]  , wgt);
    etaLooseJet -> Fill(Etajet[ij] , wgt);
    phiLooseJet -> Fill(Phijet[ij] , wgt);

    if(fabs(Etajet[ij]) < 3.0)
      {
      ptLooseJet_HB  -> Fill(Ptjet[ij]  , wgt);
      etaLooseJet_HB -> Fill(Etajet[ij] , wgt);
      phiLooseJet_HB -> Fill(Phijet[ij] , wgt);
      }
    else
      {
      ptLooseJet_HE -> Fill(Ptjet[ij]  , wgt);
      etaLooseJet_HE -> Fill(Etajet[ij] , wgt);
      phiLooseJet_HE -> Fill(Phijet[ij] , wgt);
      }
    
    if(!Tj[ij]) continue;
    ptTightJet  -> Fill(Ptjet[ij]  , wgt);
    etaTightJet -> Fill(Etajet[ij] , wgt);
    phiTightJet -> Fill(Phijet[ij] , wgt);
    
    if(fabs(Etajet[ij]) < 3.0)
      {
      ptTightJet_HB  -> Fill(Ptjet[ij]  , wgt);
      etaTightJet_HB -> Fill(Etajet[ij] , wgt);
      phiTightJet_HB -> Fill(Phijet[ij] , wgt);
      }
    else
      {
      ptTightJet_HE -> Fill(Ptjet[ij]  , wgt);
      etaTightJet_HE -> Fill(Etajet[ij] , wgt);
      phiTightJet_HE -> Fill(Phijet[ij] , wgt);
      }    
    jet_index.push_back(ij);
  }

  std::cout<<"got jet info"<<endl;

  matchedjet_index.clear();
  //unmatchedjet_index.clear();   
  //get genjet information
  std::cout<<"here1"<<endl;
  count =0;
  for(Int_t igj= 0 ; igj < *Ngenjet; igj++)    {
    if(Ptgenjet[igj] < 30.) continue;
    if(fabs(Etagenjet[igj]) > 4.7) continue;
    ptAllGenJet  -> Fill(Ptgenjet[igj]  , wgt);
    etaAllGenJet -> Fill(Etagenjet[igj] , wgt);
    phiAllGenJet -> Fill(Phigenjet[igj] , wgt);
    // GenJet Matching
    float min_dR(100.), min_pt(1000);
    Int_t index;
    std::cout<<"here2"<<endl;
    // if matching with jet_index[ij] which stores information of tight jets here.
    //for(Int_t ij= 0 ; ij < jet_index.size(); ij++) {
    for(Int_t ij= 0 ; ij < *Njet; ij++) {
      if(Ptjet[ij] < 30. ) continue;
      if(fabs(Etajet[ij]) > 4.7 ) continue;
      bool overlap=false;
      for(int ie=0;ie<*Nelec;ie++)
      {
        float dR_ejt   = dR(Etaelec[ie], Phielec[ie], Etajet[ij], Phijet[ij]);
	  float ptdiff   = fabs(Ptjet[ij] - Ptelec[ie]);
	    if(dR_ejt<0.3) {overlap=true; break;}
	    }
      if(overlap) continue;
      for(int ie=0;ie<*Nmuon;ie++)
      {
        float dR_ejt   = dR(Etamuon[ie], Phimuon[ie], Etajet[ij], Phijet[ij]);
	  float ptdiff   = fabs(Ptjet[ij] - Ptmuon[ie]);
	    if(dR_ejt<0.3) {overlap=true; break;}
	    }
      if(overlap) continue;
    
      //if(dR(Etajet[jet_index[ij]], Phijet[jet_index[ij]], Etagenjet[igj], Phigenjet[igj])< min_dR )

      if(dR(Etajet[ij], Phijet[ij], Etagenjet[igj], Phigenjet[igj])< min_dR )
      {
        min_dR = dR(Etajet[ij], Phijet[ij], Etagenjet[igj], Phigenjet[igj]);
	  min_pt = fabs(Ptgenjet[igj]-Ptjet[ij])/(Ptgenjet[igj]);
	    //index  = jet_index[ij];
	      index  = ij;
	      }
    }
    std::cout<<"here3"<<endl;
    if(min_dR < 0.2 && min_pt < 0.5){
      matchedjet_index.push_back(index);
      if()
      break;
    }
    count++;
  }
  //if(matchedjet_index.size()>count)
  //  std::cout<< "Here is the problem"<< std::endl;

  for(unsigned int ij = 0 ; ij < matchedjet_index.size(); ij++){
  if(matchedjet_index[ij] > *Njet || matchedjet_index[ij] < 0) cout<<"Issue with matchedjet_index"<<endl;
    ptGoodMatchedJet  -> Fill(Ptjet [matchedjet_index[ij]] , wgt);
    etaGoodMatchedJet -> Fill(Etajet[matchedjet_index[ij]] , wgt);
    phiGoodMatchedJet -> Fill(Phijet[matchedjet_index[ij]] , wgt);
    if(matchedjet_index.size()<=3)
    {
      //filling histograms
        ptMatchedJet[ij] ->Fill(Ptjet [matchedjet_index[ij]], wgt);
	  etaMatchedJet[ij]->Fill(Etajet[matchedjet_index[ij]], wgt);
	    phiMatchedJet[ij]->Fill(Phijet[matchedjet_index[ij]], wgt);
	    }
    }
    std::cout<<"got allgenjet info"<<endl;


  */  

  return kTRUE;
}

void SelectorClass::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  fileName->cd();
  fileName->Write();
  fileName->Close();

  //delete tmvatree;
  //delete ptAllGenMuon,etaAllGenMuon, phiAllGenMuon;
  //delete ptAllGenElec,etaAllGenElec, phiAllGenElec;
  //delete fileName;
  //delete fProofFile;
}

void SelectorClass::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
