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

   sprintf(name,"rel_isoAllElec_EB");
   rel_isoAllElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllElec_EB->Sumw2();
   rel_isoAllElec_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoAllElec_EE");
   rel_isoAllElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllElec_EE->Sumw2();
   rel_isoAllElec_EE->GetXaxis()->SetTitle("relIso");

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

   h_dRElec  = new TH1F ("h_dRElec","", 100,0,1);
   h_dRElec->Sumw2();
   h_dRElec->GetXaxis()->SetTitle("#Delta R");

   sprintf(name,"ptAllMatchedElec");
   ptAllMatchedElec = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedElec->Sumw2();
   ptAllMatchedElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaAllMatchedElec");
   etaAllMatchedElec = new TH1F (name,"" , 100, -5, 5);
   etaAllMatchedElec->Sumw2();
   etaAllMatchedElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllMatchedElec");
   phiAllMatchedElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllMatchedElec->Sumw2();
   phiAllMatchedElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptAllMatchedElec_EB");
   ptAllMatchedElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedElec_EB->Sumw2();
   ptAllMatchedElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllMatchedElec_EE");
   ptAllMatchedElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedElec_EE->Sumw2();
   ptAllMatchedElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoAllMatchedElec_EB");
   rel_isoAllMatchedElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMatchedElec_EB->Sumw2();
   rel_isoAllMatchedElec_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoAllMatchedElec_EE");
   rel_isoAllMatchedElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMatchedElec_EE->Sumw2();
   rel_isoAllMatchedElec_EE->GetXaxis()->SetTitle("relIso");


   sprintf(name,"nvtxAllMatchedElec");
   nvtxAllMatchedElec = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxAllMatchedElec->Sumw2();
   nvtxAllMatchedElec->GetXaxis()->SetTitle("#vertices");


   sprintf(name,"ptLooseMatchedElec");
   ptLooseMatchedElec = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedElec->Sumw2();
   ptLooseMatchedElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaLooseMatchedElec");
   etaLooseMatchedElec = new TH1F (name,"" , 100, -5, 5);
   etaLooseMatchedElec->Sumw2();
   etaLooseMatchedElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiLooseMatchedElec");
   phiLooseMatchedElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiLooseMatchedElec->Sumw2();
   phiLooseMatchedElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptLooseMatchedElec_EB");
   ptLooseMatchedElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedElec_EB->Sumw2();
   ptLooseMatchedElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptLooseMatchedElec_EE");
   ptLooseMatchedElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedElec_EE->Sumw2();
   ptLooseMatchedElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoLooseMatchedElec_EB");
   rel_isoLooseMatchedElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMatchedElec_EB->Sumw2();
   rel_isoLooseMatchedElec_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoLooseMatchedElec_EE");
   rel_isoLooseMatchedElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMatchedElec_EE->Sumw2();
   rel_isoLooseMatchedElec_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"nvtxLooseMatchedElec");
   nvtxLooseMatchedElec = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxLooseMatchedElec->Sumw2();
   nvtxLooseMatchedElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptMediumMatchedElec");
   ptMediumMatchedElec = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumMatchedElec->Sumw2();
   ptMediumMatchedElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaMediumMatchedElec");
   etaMediumMatchedElec = new TH1F (name,"" , 100, -5, 5);
   etaMediumMatchedElec->Sumw2();
   etaMediumMatchedElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiMediumMatchedElec");
   phiMediumMatchedElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiMediumMatchedElec->Sumw2();
   phiMediumMatchedElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptMediumMatchedElec_EB");
   ptMediumMatchedElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumMatchedElec_EB->Sumw2();
   ptMediumMatchedElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptMediumMatchedElec_EE");
   ptMediumMatchedElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptMediumMatchedElec_EE->Sumw2();
   ptMediumMatchedElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoMediumMatchedElec_EB");
   rel_isoMediumMatchedElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoMediumMatchedElec_EB->Sumw2();
   rel_isoMediumMatchedElec_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoMediumMatchedElec_EE");
   rel_isoMediumMatchedElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoMediumMatchedElec_EE->Sumw2();
   rel_isoMediumMatchedElec_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"nvtxMediumMatchedElec");
   nvtxMediumMatchedElec = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxMediumMatchedElec->Sumw2();
   nvtxMediumMatchedElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptTightMatchedElec");
   ptTightMatchedElec = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedElec->Sumw2();
   ptTightMatchedElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"rel_isoTightMatchedElec_EB");
   rel_isoTightMatchedElec_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMatchedElec_EB->Sumw2();
   rel_isoTightMatchedElec_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoTightMatchedElec_EE");
   rel_isoTightMatchedElec_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMatchedElec_EE->Sumw2();
   rel_isoTightMatchedElec_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"etaTightMatchedElec");
   etaTightMatchedElec = new TH1F (name,"" , 100, -5, 5);
   etaTightMatchedElec->Sumw2();
   etaTightMatchedElec->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiTightMatchedElec");
   phiTightMatchedElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiTightMatchedElec->Sumw2();
   phiTightMatchedElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxTightMatchedElec");
   nvtxTightMatchedElec = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxTightMatchedElec->Sumw2();
   nvtxTightMatchedElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptTightMatchedElec_EB");
   ptTightMatchedElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedElec_EB->Sumw2();
   ptTightMatchedElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptTightMatchedElec_EE");
   ptTightMatchedElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedElec_EE->Sumw2();
   ptTightMatchedElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptGoodMatchedElec");
   ptGoodMatchedElec = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedElec->Sumw2();
   ptGoodMatchedElec->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaGoodMatchedElec");
   etaGoodMatchedElec = new TH1F (name,"" , 100, -5, 5);
   etaGoodMatchedElec->Sumw2();
   etaGoodMatchedElec->GetXaxis()->SetTitle("#eta");

   sprintf(name,"phiGoodMatchedElec");
   phiGoodMatchedElec = new TH1F (name,"" , 100, -3.5, 3.5);
   phiGoodMatchedElec->Sumw2(); 
   phiGoodMatchedElec->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxGoodMatchedElec");
   nvtxGoodMatchedElec = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxGoodMatchedElec->Sumw2();
   nvtxGoodMatchedElec->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptGoodMatchedElec_EB");
   ptGoodMatchedElec_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedElec_EB->Sumw2();
   ptGoodMatchedElec_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptGoodMatchedElec_EE");
   ptGoodMatchedElec_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedElec_EE->Sumw2();
   ptGoodMatchedElec_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllGenMuon");
   ptAllGenMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptAllGenMuon->Sumw2();
   ptAllGenMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
   //ptAllGenMuon->GetNbinsX();
   //ptAllGenMuon->GetYaxis()->SetTitle("#Events/GeV");
     
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

   sprintf(name,"ptAllMatchedMuon");
   ptAllMatchedMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedMuon->Sumw2();
   ptAllMatchedMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaAllMatchedMuon");
   etaAllMatchedMuon = new TH1F (name,"" , 100, -5, 5);
   etaAllMatchedMuon->Sumw2();
   etaAllMatchedMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiAllMatchedMuon");
   phiAllMatchedMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiAllMatchedMuon->Sumw2();
   phiAllMatchedMuon->GetXaxis()->SetTitle("#phi");

   h_dRMuon  = new TH1F ("h_dRMuon","", 100,0,1);
   h_dRMuon->Sumw2();
   h_dRMuon->GetXaxis()->SetTitle("#Delta R");


   sprintf(name,"ptAllMatchedMuon_EB");
   ptAllMatchedMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedMuon_EB->Sumw2();
   ptAllMatchedMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptAllMatchedMuon_EE");
   ptAllMatchedMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptAllMatchedMuon_EE->Sumw2();
   ptAllMatchedMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoAllMatchedMuon");
   rel_isoAllMatchedMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMatchedMuon ->Sumw2();
   rel_isoAllMatchedMuon ->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoAllMatchedMuon_EB");
   rel_isoAllMatchedMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMatchedMuon_EB->Sumw2();
   rel_isoAllMatchedMuon_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoAllMatchedMuon_EE");
   rel_isoAllMatchedMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoAllMatchedMuon_EE->Sumw2();
   rel_isoAllMatchedMuon_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"nvtxAllMatchedMuon");
   nvtxAllMatchedMuon = new TH1F(name,"" , 75, 50.5, 200.5);
   nvtxAllMatchedMuon->Sumw2();
   nvtxAllMatchedMuon->GetXaxis()->SetTitle("#vertices");


   sprintf(name,"ptLooseMatchedMuon");
   ptLooseMatchedMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedMuon->Sumw2();
   ptLooseMatchedMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaLooseMatchedMuon");
   etaLooseMatchedMuon = new TH1F (name,"" , 100, -5, 5);
   etaLooseMatchedMuon->Sumw2();
   etaLooseMatchedMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiLooseMatchedMuon");
   phiLooseMatchedMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiLooseMatchedMuon->Sumw2();
   phiLooseMatchedMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"ptLooseMatchedMuon_EB");
   ptLooseMatchedMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedMuon_EB->Sumw2();
   ptLooseMatchedMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptLooseMatchedMuon_EE");
   ptLooseMatchedMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptLooseMatchedMuon_EE->Sumw2();
   ptLooseMatchedMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"rel_isoLooseMatchedMuon");
   rel_isoLooseMatchedMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMatchedMuon ->Sumw2();
   rel_isoLooseMatchedMuon ->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoLooseMatchedMuon_EB");
   rel_isoLooseMatchedMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMatchedMuon_EB->Sumw2();
   rel_isoLooseMatchedMuon_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoLooseMatchedMuon_EE");
   rel_isoLooseMatchedMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoLooseMatchedMuon_EE->Sumw2();
   rel_isoLooseMatchedMuon_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"nvtxLooseMatchedMuon");
   nvtxLooseMatchedMuon = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxLooseMatchedMuon->Sumw2();
   nvtxLooseMatchedMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptTightMatchedMuon");
   ptTightMatchedMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedMuon->Sumw2();
   ptTightMatchedMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaTightMatchedMuon");
   etaTightMatchedMuon = new TH1F (name,"" , 100, -5, 5);
   etaTightMatchedMuon->Sumw2();
   etaTightMatchedMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiTightMatchedMuon");
   phiTightMatchedMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiTightMatchedMuon->Sumw2();
   phiTightMatchedMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxTightMatchedMuon");
   nvtxTightMatchedMuon = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxTightMatchedMuon->Sumw2();
   nvtxTightMatchedMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptTightMatchedMuon_EB");
   ptTightMatchedMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedMuon_EB->Sumw2();
   ptTightMatchedMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptTightMatchedMuon_EE");
   ptTightMatchedMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptTightMatchedMuon_EE->Sumw2();
   ptTightMatchedMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");


   sprintf(name,"rel_isoTightMatchedMuon");
   rel_isoTightMatchedMuon  = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMatchedMuon ->Sumw2();
   rel_isoTightMatchedMuon ->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoTightMatchedMuon_EB");
   rel_isoTightMatchedMuon_EB = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMatchedMuon_EB->Sumw2();
   rel_isoTightMatchedMuon_EB->GetXaxis()->SetTitle("relIso");

   sprintf(name,"rel_isoTightMatchedMuon_EE");
   rel_isoTightMatchedMuon_EE = new TH1F (name,"" , 100, 0, 1);
   rel_isoTightMatchedMuon_EE->Sumw2();
   rel_isoTightMatchedMuon_EE->GetXaxis()->SetTitle("relIso");

   sprintf(name,"ptGoodMatchedMuon");
   ptGoodMatchedMuon = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedMuon->Sumw2();
   ptGoodMatchedMuon->GetXaxis()->SetTitle("P_{t}(GeV)");
     
   sprintf(name,"etaGoodMatchedMuon");
   etaGoodMatchedMuon = new TH1F (name,"" , 100, -5, 5);
   etaGoodMatchedMuon->Sumw2();
   etaGoodMatchedMuon->GetXaxis()->SetTitle("#eta");
   
   sprintf(name,"phiGoodMatchedMuon");
   phiGoodMatchedMuon = new TH1F (name,"" , 100, -3.5, 3.5);
   phiGoodMatchedMuon->Sumw2(); 
   phiGoodMatchedMuon->GetXaxis()->SetTitle("#phi");

   sprintf(name,"nvtxGoodMatchedMuon");
   nvtxGoodMatchedMuon = new TH1F (name,"" , 75, 50.5, 200.5);
   nvtxGoodMatchedMuon->Sumw2();
   nvtxGoodMatchedMuon->GetXaxis()->SetTitle("#vertices");

   sprintf(name,"ptGoodMatchedMuon_EB");
   ptGoodMatchedMuon_EB = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedMuon_EB->Sumw2();
   ptGoodMatchedMuon_EB->GetXaxis()->SetTitle("P_{t}(GeV)");

   sprintf(name,"ptGoodMatchedMuon_EE");
   ptGoodMatchedMuon_EE = new TH1F (name,"" ,  100, 0, 1000);
   ptGoodMatchedMuon_EE->Sumw2();
   ptGoodMatchedMuon_EE->GetXaxis()->SetTitle("P_{t}(GeV)");

}

Bool_t SelectorClass::Process(Long64_t entry)
{
  fReader.SetEntry(entry);
  std::cout<<"Event:"<<entry<<std::endl;
  nvtxAll -> Fill((float) *Nvtx , wgt);

  matchedpho_index.clear();
  matchedLoosepho_index.clear();
  matchedMediumpho_index.clear();
  matchedTightpho_index.clear();
  matchedGoodpho_index.clear();
  
  
  
  
  
  for(Int_t ip= 0 ; ip < *Npho; ip++) {
    if(Ptpho[ip] < 20. ) continue;
    if(fabs(Etapho[ip]) > 3.0 ) continue;
    if(fabs(Etapho[ip])> 1.444 && fabs(Etapho[ip])< 1.566 ) continue;
    if(fabs(Etapho[ip]) < 1.444)
      rel_isoAllPho_EB -> Fill(IsolationVarpho[ip], wgt);
    if(fabs(Etapho[ip]) > 1.566)
      rel_isoAllPho_EE -> Fill(IsolationVarpho[ip], wgt);
  }


  for(Int_t ig= 0 ; ig < *Ngenphoton; ig++)
    {
      if (Ptgenphoton[ig] < 20. ) continue;
      if(fabs(Etagenphoton[ig]) > 3.0 ) continue;
      if(fabs(Etagenphoton[ig]) > 1.444 && fabs(Etagenphoton[ig]) < 1.566 ) continue;
      //rel_isoAllGenPho -> Fill(IsolationVargenphoton[ig] , wgt);
      if (IsolationVargenphoton[ig] > 0.02 ) continue;
      
      ptAllGenPho  -> Fill(Ptgenphoton[ig]  , wgt);
      etaAllGenPho -> Fill(Etagenphoton[ig] , wgt);
      phiAllGenPho -> Fill(Phigenphoton[ig] , wgt);
      //Gen Pho matching
      float min_dR(100.), min_pt(1000);
      for(Int_t ij= 0 ; ij < *Npho; ij++) 
	{
	  if(Ptpho[ij] < 20. ) continue;
	  if(fabs(Etapho[ij]) > 3.0 ) continue;
	  if(fabs(Etapho[ij])> 1.444 && fabs(Etapho[ij])< 1.566 ) continue;
	  if(dR(Etapho[ij], Phipho[ij], Etagenphoton[ig], Phigenphoton[ig])< min_dR )
	    {
	      min_dR = dR(Etapho[ij], Phipho[ij], Etagenphoton[ig], Phigenphoton[ig]);
	      min_pt = fabs(Ptgenphoton[ig]-Ptpho[ij])/(Ptgenphoton[ig]);
	    }
	  
	  h_dRPho->Fill(min_dR);
	  if(min_dR < 0.1 && min_pt < 0.5)
	    {
	      matchedpho_index.push_back(ij);
	      
	      if(isLpho[ij])
		matchedLoosepho_index.push_back(ij);
	      
	      if(isTpho[ij])
		matchedTightpho_index.push_back(ij);
	      break;
	    }
	}
    }
  
  for(unsigned int ij = 0 ; ij < matchedpho_index.size(); ij++){
    ptAllMatchedPho  -> Fill(Ptpho [matchedpho_index[ij]] , wgt);
    etaAllMatchedPho -> Fill(Etapho[matchedpho_index[ij]] , wgt);
    phiAllMatchedPho -> Fill(Phipho[matchedpho_index[ij]] , wgt);
    nvtxAllMatchedPho -> Fill((float) *Nvtx , wgt);
    if(fabs(Etapho[matchedpho_index[ij]]) < 1.444)
      { 
	ptAllMatchedPho_EB       -> Fill(Ptpho [matchedpho_index[ij]] , wgt);
	rel_isoAllMatchedPho_EB  -> Fill(IsolationVarpho[matchedpho_index[ij]] , wgt);
      }
    else
      {
	ptAllMatchedPho_EE  -> Fill(Ptpho [matchedpho_index[ij]] , wgt);
	rel_isoAllMatchedPho_EE  -> Fill(IsolationVarpho[matchedpho_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedLoosepho_index.size(); ij++){
    ptLooseMatchedPho  -> Fill(Ptpho [matchedLoosepho_index[ij]] , wgt);
    etaLooseMatchedPho -> Fill(Etapho[matchedLoosepho_index[ij]] , wgt);
    phiLooseMatchedPho -> Fill(Phipho[matchedLoosepho_index[ij]] , wgt);
    nvtxLooseMatchedPho -> Fill((float)*Nvtx , wgt);
    if(fabs(Etapho[matchedLoosepho_index[ij]]) < 1.444)
      {
	ptLooseMatchedPho_EB  -> Fill(Ptpho [matchedLoosepho_index[ij]] , wgt);
	rel_isoLooseMatchedPho_EB  -> Fill(IsolationVarpho[matchedLoosepho_index[ij]] , wgt);
      }
    else
      {
	ptLooseMatchedPho_EE  -> Fill(Ptpho [matchedLoosepho_index[ij]] , wgt);
	rel_isoLooseMatchedPho_EE  -> Fill(IsolationVarpho[matchedLoosepho_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedTightpho_index.size(); ij++){
    ptTightMatchedPho  -> Fill(Ptpho [matchedTightpho_index[ij]] , wgt);
    etaTightMatchedPho -> Fill(Etapho[matchedTightpho_index[ij]] , wgt);
    phiTightMatchedPho -> Fill(Phipho[matchedTightpho_index[ij]] , wgt);
    nvtxTightMatchedPho -> Fill((float)*Nvtx , wgt);
    if(fabs(Etapho[matchedTightpho_index[ij]]) < 1.444)
      {
	ptTightMatchedPho_EB  -> Fill(Ptpho [matchedTightpho_index[ij]] , wgt);
	rel_isoTightMatchedPho_EB  -> Fill(IsolationVarpho[matchedTightpho_index[ij]] , wgt);
      }
    else
      {
	ptTightMatchedPho_EE  -> Fill(Ptpho [matchedTightpho_index[ij]] , wgt);
	rel_isoTightMatchedPho_EE  -> Fill(IsolationVarpho[matchedTightpho_index[ij]] , wgt);
      }
  }


  //// For Electrons
  matchedelec_index.clear();
  matchedLooseelec_index.clear();
  matchedMediumelec_index.clear();
  matchedTightelec_index.clear();
  matchedGoodelec_index.clear();

  //nvtxAll -> Fill((float) *Nvtx , wgt);
  for(Int_t ie= 0 ; ie < *Nelec; ie++) {
    if(Ptelec[ie] < 20. ) continue;
    if(fabs(Etaelec[ie]) > 3.0 ) continue;
    if(fabs(Etaelec[ie])> 1.444 && fabs(Etaelec[ie])< 1.566 ) continue;
    if(fabs(Etaelec[ie]) < 1.444)
      rel_isoAllElec_EB -> Fill(IsolationVarelec[ie], wgt);
    if(fabs(Etaelec[ie]) > 1.566)
      rel_isoAllElec_EE -> Fill(IsolationVarelec[ie], wgt);
  }


  for(Int_t ig= 0 ; ig < *Ngenlepton; ig++) {
    if(abs(PIdgenlepton[ig])==11)
      {
	if (Ptgenlepton[ig] < 20. ) continue;
	if(fabs(Etagenlepton[ig]) > 3.0 ) continue;
	if(fabs(Etagenlepton[ig]) > 1.444 && fabs(Etagenlepton[ig]) < 1.566 ) continue;
	//rel_isoAllGenElec -> Fill(IsolationVargenlepton[ig] , wgt);
	if (IsolationVargenlepton[ig] > 0.02 ) continue;

	ptAllGenElec  -> Fill(Ptgenlepton[ig]  , wgt);
	etaAllGenElec -> Fill(Etagenlepton[ig] , wgt);
	phiAllGenElec -> Fill(Phigenlepton[ig] , wgt);
	//Gen Elec matching
	float min_dR(100.), min_pt(1000);
	for(Int_t ij= 0 ; ij < *Nelec; ij++) {
	  if(Ptelec[ij] < 20. ) continue;
	  if(fabs(Etaelec[ij]) > 3.0 ) continue;
	  if(fabs(Etaelec[ij])> 1.444 && fabs(Etaelec[ij])< 1.566 ) continue;
	  if(dR(Etaelec[ij], Phielec[ij], Etagenlepton[ig], Phigenlepton[ig])< min_dR )
	    {
	      min_dR = dR(Etaelec[ij], Phielec[ij], Etagenlepton[ig], Phigenlepton[ig]);
	      min_pt = fabs(Ptgenlepton[ig]-Ptelec[ij])/(Ptgenlepton[ig]);
	    }

	  h_dRElec->Fill(min_dR);
	  if(min_dR < 0.1 && min_pt < 0.5) {
	    matchedelec_index.push_back(ij);
	    
	    if(isLE[ij])
	      matchedLooseelec_index.push_back(ij);
	    
	    if(isME[ij])
	      matchedMediumelec_index.push_back(ij);
	    
	    if(isTE[ij])
	      matchedTightelec_index.push_back(ij);
	    
	    if(isME[ij])
	      { 
		if((fabs(Etaelec[ij]) < 1.444 && IsolationVarelec[ij] <  0.15)
		   ||(fabs(Etaelec[ij]) > 1.566 && IsolationVarelec[ij] < 0.21))
		  matchedGoodelec_index.push_back(ij);
	      }
	    break;
	  }
	}
      }
  }

  for(unsigned int ij = 0 ; ij < matchedelec_index.size(); ij++){
    ptAllMatchedElec  -> Fill(Ptelec [matchedelec_index[ij]] , wgt);
    etaAllMatchedElec -> Fill(Etaelec[matchedelec_index[ij]] , wgt);
    phiAllMatchedElec -> Fill(Phielec[matchedelec_index[ij]] , wgt);
    nvtxAllMatchedElec -> Fill((float) *Nvtx , wgt);
    if(fabs(Etaelec[matchedelec_index[ij]]) < 1.444)
      { 
	ptAllMatchedElec_EB       -> Fill(Ptelec [matchedelec_index[ij]] , wgt);
	rel_isoAllMatchedElec_EB  -> Fill(IsolationVarelec[matchedelec_index[ij]] , wgt);
	//etaAllMatchedElec_EB -> Fill(Etaelec[matchedelec_index[ij]] , wgt);
	//phiAllMatchedElec_EB  -> Fill(Ptelec [matchedelec_index[ij]] , wgt);
	//nvtxAllMatchedElec_EB -> Fill(Etaelec[matchedelec_index[ij]] , wgt);
      }
    else
      {
	ptAllMatchedElec_EE  -> Fill(Ptelec [matchedelec_index[ij]] , wgt);
	rel_isoAllMatchedElec_EE  -> Fill(IsolationVarelec[matchedelec_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedLooseelec_index.size(); ij++){
    ptLooseMatchedElec  -> Fill(Ptelec [matchedLooseelec_index[ij]] , wgt);
    etaLooseMatchedElec -> Fill(Etaelec[matchedLooseelec_index[ij]] , wgt);
    phiLooseMatchedElec -> Fill(Phielec[matchedLooseelec_index[ij]] , wgt);
    nvtxLooseMatchedElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[matchedLooseelec_index[ij]]) < 1.444)
      {
	ptLooseMatchedElec_EB  -> Fill(Ptelec [matchedLooseelec_index[ij]] , wgt);
	rel_isoLooseMatchedElec_EB  -> Fill(IsolationVarelec[matchedLooseelec_index[ij]] , wgt);
      }
    else
      {
	ptLooseMatchedElec_EE  -> Fill(Ptelec [matchedLooseelec_index[ij]] , wgt);
	rel_isoLooseMatchedElec_EE  -> Fill(IsolationVarelec[matchedLooseelec_index[ij]] , wgt);
      }
  }

 for(unsigned int ij = 0 ; ij < matchedMediumelec_index.size(); ij++){
    ptMediumMatchedElec  -> Fill(Ptelec [matchedMediumelec_index[ij]] , wgt);
    etaMediumMatchedElec -> Fill(Etaelec[matchedMediumelec_index[ij]] , wgt);
    phiMediumMatchedElec -> Fill(Phielec[matchedMediumelec_index[ij]] , wgt);
    nvtxMediumMatchedElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[matchedMediumelec_index[ij]]) < 1.444)
      {
	ptMediumMatchedElec_EB  -> Fill(Ptelec [matchedMediumelec_index[ij]] , wgt);
	rel_isoMediumMatchedElec_EB  -> Fill(IsolationVarelec[matchedMediumelec_index[ij]] , wgt);
      }
    else
      {
	ptMediumMatchedElec_EE  -> Fill(Ptelec [matchedMediumelec_index[ij]] , wgt);
	rel_isoMediumMatchedElec_EE  -> Fill(IsolationVarelec[matchedMediumelec_index[ij]] , wgt);
      }
  }


  for(unsigned int ij = 0 ; ij < matchedTightelec_index.size(); ij++){
    ptTightMatchedElec  -> Fill(Ptelec [matchedTightelec_index[ij]] , wgt);
    etaTightMatchedElec -> Fill(Etaelec[matchedTightelec_index[ij]] , wgt);
    phiTightMatchedElec -> Fill(Phielec[matchedTightelec_index[ij]] , wgt);
    nvtxTightMatchedElec -> Fill((float)*Nvtx , wgt);
    if(fabs(Etaelec[matchedTightelec_index[ij]]) < 1.444)
      {
	ptTightMatchedElec_EB  -> Fill(Ptelec [matchedTightelec_index[ij]] , wgt);
	rel_isoTightMatchedElec_EB  -> Fill(IsolationVarelec[matchedTightelec_index[ij]] , wgt);
      }
    else
      {
	ptTightMatchedElec_EE  -> Fill(Ptelec [matchedTightelec_index[ij]] , wgt);
	rel_isoTightMatchedElec_EE  -> Fill(IsolationVarelec[matchedTightelec_index[ij]] , wgt);
      }
  }
  for(unsigned int ij = 0 ; ij < matchedGoodelec_index.size(); ij++){
    ptGoodMatchedElec  -> Fill(Ptelec [matchedGoodelec_index[ij]] , wgt);
    etaGoodMatchedElec -> Fill(Etaelec[matchedGoodelec_index[ij]] , wgt);
    phiGoodMatchedElec -> Fill(Phielec[matchedGoodelec_index[ij]] , wgt);
    nvtxGoodMatchedElec -> Fill((float) *Nvtx , wgt);
    if(fabs(Etaelec[matchedGoodelec_index[ij]]) < 1.444)
      ptGoodMatchedElec_EB  -> Fill(Ptelec [matchedGoodelec_index[ij]] , wgt);
    else
      ptGoodMatchedElec_EE  -> Fill(Ptelec [matchedGoodelec_index[ij]] , wgt);
  }


  matchedmuon_index.clear();
  matchedLoosemuon_index.clear();
  matchedTightmuon_index.clear();
  matchedGoodmuon_index.clear();

  //for(Int_t ie= 0 ; ie < *Nmuon; ie++) {
  //  if(Ptmuon[ie] < 30. ) continue;
  //  if(fabs(Etamuon[ie]) > 2.8 ) continue;
  //  //if(fabs(Etamuon[ie])> 1.444 && fabs(Etamuon[ie])< 1.566 ) continue;
  //  //if(fabs(Etamuon[ie]) < 1.444)
  //  //  rel_isoAllMuon_EB -> Fill(IsolationVarmuon[ie], wgt);
  //  //if(fabs(Etamuon[ie]) > 1.566)
  //  //  rel_isoAllMuon_EE -> Fill(IsolationVarmuon[ie], wgt);
  //}

  for(Int_t ig= 0 ; ig < *Ngenlepton; ig++) {
    if(abs(PIdgenlepton[ig])==13)
      {
	if (Ptgenlepton[ig] < 20. ) continue;
	if(fabs(Etagenlepton[ig]) > 3.0 ) continue;
	//if(fabs(Etagenlepton[ig]) > 1.444 && fabs(Etagenlepton[ig]) < 1.566 ) continue;
	//rel_isoAllGenMuon -> Fill(IsolationVargenlepton[ig] , wgt);
	if (IsolationVargenlepton[ig] > 0.02 ) continue;

	ptAllGenMuon  -> Fill(Ptgenlepton[ig]  , wgt);
	etaAllGenMuon -> Fill(Etagenlepton[ig] , wgt);
	phiAllGenMuon -> Fill(Phigenlepton[ig] , wgt);
	//Gen Muon matching
	float min_dR(100.), min_pt(1000);
	for(Int_t ij= 0 ; ij < *Nmuon; ij++) {
	  if(Ptmuon[ij] < 20. ) continue;
	  if(fabs(Etamuon[ij]) > 3.0 ) continue;
	  //if(fabs(Etamuon[ij])> 1.444 && fabs(Etamuon[ij])< 1.566 ) continue;
	  if(dR(Etamuon[ij], Phimuon[ij], Etagenlepton[ig], Phigenlepton[ig])< min_dR )
	    {
	      min_dR = dR(Etamuon[ij], Phimuon[ij], Etagenlepton[ig], Phigenlepton[ig]);
	      min_pt = fabs(Ptgenlepton[ig]-Ptmuon[ij])/(Ptgenlepton[ig]);
	    }

	  h_dRMuon->Fill(min_dR);
	  if(min_dR < 0.1 && min_pt < 0.5) {
	    matchedmuon_index.push_back(ij);
	    
	    if(isLM[ij])
	      matchedLoosemuon_index.push_back(ij);

	    if(isTM[ij])
	      matchedTightmuon_index.push_back(ij);

	    if(isTM[ij]&&IsolationVarmuon[ij] < 0.15 )
	      matchedGoodmuon_index.push_back(ij);
	    

	    
	    break;
	  }
	}
      }
  }

  for(unsigned int ij = 0 ; ij < matchedmuon_index.size(); ij++){
    ptAllMatchedMuon  -> Fill(Ptmuon [matchedmuon_index[ij]] , wgt);
    etaAllMatchedMuon -> Fill(Etamuon[matchedmuon_index[ij]] , wgt);
    phiAllMatchedMuon -> Fill(Phimuon[matchedmuon_index[ij]] , wgt);
    nvtxAllMatchedMuon -> Fill((float) *Nvtx , wgt);
    rel_isoAllMatchedMuon -> Fill(IsolationVarmuon[matchedmuon_index[ij]] , wgt);
    if(fabs(Etamuon[matchedmuon_index[ij]]) < 1.444)
      { 
	ptAllMatchedMuon_EB       -> Fill(Ptmuon [matchedmuon_index[ij]] , wgt);
	rel_isoAllMatchedMuon_EB  -> Fill(IsolationVarmuon[matchedmuon_index[ij]] , wgt);
	//etaAllMatchedMuon_EB -> Fill(Etamuon[matchedmuon_index[ij]] , wgt);
	//phiAllMatchedMuon_EB  -> Fill(Ptmuon [matchedmuon_index[ij]] , wgt);
	//nvtxAllMatchedMuon_EB -> Fill(Etamuon[matchedmuon_index[ij]] , wgt);
      }
    else
      {
	ptAllMatchedMuon_EE  -> Fill(Ptmuon [matchedmuon_index[ij]] , wgt);
	rel_isoAllMatchedMuon_EE  -> Fill(IsolationVarmuon[matchedmuon_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedLoosemuon_index.size(); ij++){
    ptLooseMatchedMuon  -> Fill(Ptmuon [matchedLoosemuon_index[ij]] , wgt);
    etaLooseMatchedMuon -> Fill(Etamuon[matchedLoosemuon_index[ij]] , wgt);
    phiLooseMatchedMuon -> Fill(Phimuon[matchedLoosemuon_index[ij]] , wgt);
    nvtxLooseMatchedMuon -> Fill((float)*Nvtx , wgt);
    rel_isoLooseMatchedMuon  -> Fill(IsolationVarmuon[matchedmuon_index[ij]] , wgt);
    if(fabs(Etamuon[matchedLoosemuon_index[ij]]) < 1.444)
      {
	ptLooseMatchedMuon_EB  -> Fill(Ptmuon [matchedLoosemuon_index[ij]] , wgt);
	rel_isoLooseMatchedMuon_EB  -> Fill(IsolationVarmuon[matchedLoosemuon_index[ij]] , wgt);
      }
    else
      {
	ptLooseMatchedMuon_EE  -> Fill(Ptmuon [matchedLoosemuon_index[ij]] , wgt);
	rel_isoLooseMatchedMuon_EE  -> Fill(IsolationVarmuon[matchedLoosemuon_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedTightmuon_index.size(); ij++){
    ptTightMatchedMuon  -> Fill(Ptmuon [matchedTightmuon_index[ij]] , wgt);
    etaTightMatchedMuon -> Fill(Etamuon[matchedTightmuon_index[ij]] , wgt);
    phiTightMatchedMuon -> Fill(Phimuon[matchedTightmuon_index[ij]] , wgt);
    nvtxTightMatchedMuon -> Fill((float)*Nvtx , wgt);
    rel_isoTightMatchedMuon  -> Fill(IsolationVarmuon[matchedTightmuon_index[ij]] , wgt);
    if(fabs(Etamuon[matchedTightmuon_index[ij]]) < 1.444)
      {
	ptTightMatchedMuon_EB  -> Fill(Ptmuon [matchedTightmuon_index[ij]] , wgt);
	rel_isoTightMatchedMuon_EB  -> Fill(IsolationVarmuon[matchedTightmuon_index[ij]] , wgt);
      }
    else
      {
	ptTightMatchedMuon_EE  -> Fill(Ptmuon [matchedTightmuon_index[ij]] , wgt);
	rel_isoTightMatchedMuon_EE  -> Fill(IsolationVarmuon[matchedTightmuon_index[ij]] , wgt);
      }
  }

  for(unsigned int ij = 0 ; ij < matchedGoodmuon_index.size(); ij++){
    ptGoodMatchedMuon  -> Fill(Ptmuon [matchedGoodmuon_index[ij]] , wgt);
    etaGoodMatchedMuon -> Fill(Etamuon[matchedGoodmuon_index[ij]] , wgt);
    phiGoodMatchedMuon -> Fill(Phimuon[matchedGoodmuon_index[ij]] , wgt);
    nvtxGoodMatchedMuon -> Fill((float) *Nvtx , wgt);
    if(fabs(Etamuon[matchedGoodmuon_index[ij]]) < 1.444)
      {
	ptGoodMatchedMuon_EB  -> Fill(Ptmuon [matchedGoodmuon_index[ij]] , wgt);
      }
    else
      {
	ptGoodMatchedMuon_EE  -> Fill(Ptmuon [matchedGoodmuon_index[ij]] , wgt);
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
        relpt    -> Fill(ptdiff/Ptelec[ie]);
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
